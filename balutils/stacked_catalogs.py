from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt

# NOTE: Try using generators for Table chunks if files get too large!
# http://docs.astropy.org/en/stable/io/ascii/read.html#reading-large-tables-in-chunks

# NOTE: Could also switch to pandas DF as base catalog type, though I'm not sure
# matters much if we're careful about what we load into memory in the first place
# (i.e. using `cols`)

# NOTE: `Catalog` should relaly be an abstract class, but something in how
# I'm using ABCMeta is causing problems. Can fix later if we care
# class Catalog(ABCMeta):
class Catalog(object):

    def __init__(self, filename, cols=None):
        self.filename = filename
        self.cols = cols

        self._load_catalog()

        return

#     @abstractmethod
    def _load_catalog(self):
        pass

    @staticmethod
    def flux2mag(flux, zp=30.):
        return -2.5 * np.log10(flux) + zp

    def apply_cut(self, cut):
        self._cat = self._cat[cut]
        self.Nobjs = len(self._cat)

        return

    def get_cat(self):
        return self._cat

    # The following are so we can access the catalog
    # values similarly to a dict
    def __getitem__(self, key):
        return self._cat[key]

    def __setitem__(self, key, value):
        self._cat[key] = value

    def __delitem__(self, key):
        del self._cat[key]

    def __contains__(self, key):
        return key in self._cat

    def __len__(self):
        return len(self._cat)

    def __repr__(self):
        return repr(self._cat)

class FitsCatalog(Catalog):
    def _load_catalog(self):
        self._cat = Table(fitsio.read(self.filename, columns=self.cols))
        self.Nobjs = len(self._cat)

        return

class DetectionCatalog(FitsCatalog):

    def __init__(self, filename, cols=None):
        super(DetectionCatalog, self).__init__(filename, cols=cols)

        self._check_for_duplicates()

        return

    def _check_for_duplicates(self):
        '''
        Balrog stack versions 1.4 and below have a small bug that
        seems to duplicate exactly 1 object, so check for these

        NOTE: Only works if there are exactly 1 extra duplicate for
        a given bal_id!
        '''
        unq, unq_idx, unq_cnt = np.unique(self._cat['bal_id'],
                                          return_inverse=True,
                                          return_counts=True)
        Nunq = len(unq)
        if Nunq != self.Nobjs:
            Ndups = self.Nobjs - Nunq
            dup_ids = unq[np.where(unq_cnt > 1)]
            print('Warning: Detection catalog has {} duplicate(s)!'.format(Ndups))
            print('Removing the following duplicates from detection catalog:')
            print(dup_ids)

        Nbefore = self.Nobjs
        for did in dup_ids:
            indx = np.where(self._cat['bal_id']==did)[0]
            print 'indx',indx
            self._cat.remove_row(indx[0])

        self.Nobjs = len(self._cat)
        assert self.Nobjs == (Nbefore - Ndups)

        print('{} duplicates removed, catalog size now {}'.format(Ndups, self.Nobjs))

        return

class H5Catalog(Catalog):

    def __init__(self, filename, basepath, cols=None):
        self.basepath = basepath
        super(H5Catalog, self).__init__(filename, cols=cols)

        return

    def _load_catalog(self):
        self._h5cat = h5py.File(self.filename)
        self._cat = Table()

        if self.cols is not None:
            for col in self.cols:
                path = os.path.join(self.basepath, col)
                self._cat[col] = self._h5cat[path][:]

        self.Nobjs = len(self._cat)

        return

    def add_col(self, col):
        path = os.path.join(self.basepath, col)
        self._cat[col] = self._h5cat[path]

        return

    def delete_col(self, col):
        self._cat.remove_column(col)

        return

    def __delete__(self):
        self._h5cat.close()
        super(Catalog, self).__delete__()

        return

class McalCatalog(H5Catalog):
    pass

class BalrogMcalCatalog(Catalog):

    _gold_cut_cols = ['flags_foreground',
                      'flags_badregions',
                      'flags_footprint',
                      'meas_FLAGS_GOLD'
                     ]

    _shape_cut_cols = ['flags',
                       'T',
                       'psf_T',
                       'snr'
                      ]

    def __init__(self, mcal_file, det_file, mcal_cols=None, det_cols=None,
                 mcal_path='catalog/unsheared', save_all=False, vb=False):

        self.mcal_file = mcal_file
        self.det_file = det_file
        self.mcal_cols = mcal_cols
        self.det_cols = det_cols
        self.mcal_path = mcal_path
        self.save_all = save_all
        self.vb = vb

        self._load_catalog()

        return

    def _load_catalog(self):
        if self.vb is True: print('Loading Mcal catalog...')
        mcal = H5Catalog(self.mcal_file, self.mcal_path, cols=self.mcal_cols)

        if self.vb is True: print('Loading Detection catalog...')
        det  = DetectionCatalog(self.det_file, cols=self.det_cols)

        if self.vb is True: print('Joining catalogs...')
        self._join(mcal.get_cat(), det.get_cat())

        return


    def _join(self, mcal, det):
        self._cat = join(mcal, det, join_type='left')

        if self.save_all is True:
            self.mcal = mcal
            self.det = det
        else:
            self.mcal = None
            self.det = None

        return

    def apply_gold_cuts(self):
        self._check_for_cols(self._gold_cut_cols)
        gold_cuts = np.where( (self._cat['flags_foreground'] == 0) &
                              (self._cat['flags_badregions'] < 2) &
                              (self._cat['flags_footprint'] == 1) &
                              (self._cat['meas_FLAGS_GOLD'] < 2)
                            )
        self.apply_cut(gold_cuts)

        return

    def apply_shape_cuts(self):
        self._check_for_cols(self._shape_cut_cols)
        shape_cuts = np.where( (self._cat['flags'] == 0) &
                              ((self._cat['T']/self._cat['psf_T']) > 0.5) &
                              (self._cat['snr'] > 10) &
                              (self._cat['snr'] < 100)
                            )
        self.apply_cut(shape_cuts)

        return

    def _check_for_cols(self, cols):
        for col in cols:
            if col not in self._cat.colnames:
                raise AttributeError('{} not found in joined '.format(col) +
                                     'catalog but required for requested cuts!')
