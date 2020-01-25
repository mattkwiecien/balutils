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

    # @abstractmethod
    def _load_catalog(self):
        pass

    @staticmethod
    def flux2mag(flux, zp=30., clip_val=0.001):
        return -2.5 * np.log10(flux.clip(clip_val)) + zp

    def apply_cut(self, cut):
        self._cat = self._cat[cut]
        self.Nobjs = len(self._cat)

        return

    def get_cat(self):
        return self._cat

    def fill_cat(self):
        self._cat = self._cat.filled()

    def _check_for_cols(self, cols):
        for col in cols:
            if col not in self._cat.colnames:
                raise AttributeError('{} not found in joined '.format(col) +
                                     'catalog but required for requested cuts!')

        return

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

class GoldCatalog(Catalog):
    _gold_cut_cols_default = ['flags_foreground',
                              'flags_badregions',
                              'flags_footprint',
                              'meas_FLAGS_GOLD'
                              ]
    _gold_cut_cols_mof_only = ['flags_foreground',
                               'flags_badregions',
                               'flags_footprint',
                               'meas_FLAGS_GOLD_MOF_ONLY'
                               ]
    _gold_cut_cols_sof_only = ['flags_foreground',
                               'flags_badregions',
                               'flags_footprint',
                               'meas_FLAGS_GOLD_SOF_ONLY'
                               ]

    _gold_cut_cols = {'default':_gold_cut_cols_default,
                      'mof_only':_gold_cut_cols_mof_only,
                      'sof_only':_gold_cut_cols_sof_only
                      }

    def __init__(self, filename, cols=None, match_type='default'):
        super(GoldCatalog, self).__init__(filename, cols=None)

        self.match_type = match_type
        self._set_gold_colname(match_type)

        return

    def _set_gold_colname(self, match_type):
        if match_type == 'default':
            self.flags_gold_colname = 'meas_FLAGS_GOLD'
        elif match_type == 'mof_only':
            self.flags_gold_colname = 'meas_FLAGS_GOLD_MOF_ONLY'
        elif match_type == 'sof_only':
            self.flags_gold_colname = 'meas_FLAGS_GOLD_SOF_ONLY'
        else:
            raise ValueError('match_type can only be default, mof_only, or sof_only')

        return

    def apply_gold_cuts(self):
        self._check_for_cols(self._gold_cut_cols[self.match_type])
        gold_cuts = np.where( (self._cat['flags_foreground'] == 0) &
                              (self._cat['flags_badregions'] < 2) &
                              (self._cat['flags_footprint'] == 1) &
                              (self._cat[self.flags_gold_colname] < 2)
                            )
        self.apply_cut(gold_cuts)

        return

    pass

class FitsCatalog(Catalog):
    def _load_catalog(self):
        self._cat = Table(fitsio.read(self.filename, columns=self.cols))
        self.Nobjs = len(self._cat)

        return

# TODO: Remove if not useful
class GoldFitsCatalog(FitsCatalog, GoldCatalog):
    def __init__(self, filename, cols=None, match_type='default'):
        super(GoldFitsCatalog, self).__init__(filename, cols=cols, match_type=match_type)

        return

class DetectionCatalog(FitsCatalog, GoldCatalog):

    def __init__(self, filename, cols=None, match_type='default'):
        super(DetectionCatalog, self).__init__(filename, cols=cols, match_type=match_type)

        self._check_for_duplicates()

        return

    def _check_for_duplicates(self):
        '''
        Balrog stack versions 1.4 and below have a small bug that
        seems to duplicate exactly 1 object, so check for these
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

                L = len(indx)
                for i in range(L-1): # keep last one
                    self._cat.remove_row(indx[i])

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

    _match_flag_col = 'match_flag_1.5_asec'

    _shape_cut_cols = ['flags',
                       'T',
                       'psf_T',
                       'snr'
                      ]
    _sompz_cut_cols = ['mag_r',
                       'mag_i',
                       'mag_z',
                       _match_flag_col
                      ]

    _gap_flux_cols = ['e_1',
                      'e_2',
                      'T',
                      'flux_r',
                      'flux_i',
                      'flux_z'
                      ]

    def __init__(self, filename, basepath, cols=None, compute_gap_fluxes=False, vb=False):
        super(McalCatalog, self).__init__(filename, basepath, cols=cols)

        self.calc_mags()

        if compute_gap_fluxes is True:
            self._compute_gap_fluxes(vb)

        return

    def calc_mags(self):
        '''
        Mcal catalogs don't automatically come with magnitudes
        '''

        fluxes = [c for c in self._cat.colnames if 'flux' in c.lower()]
        bands = [f[-1] for f in fluxes]

        for b in bands:
            self._cat['mag_{}'.format(b)]= self.flux2mag(self._cat['flux_{}'.format(b)])

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

    def apply_sompz_cuts(self, use_match_flag=True):
        '''
        This is the same as "selection2"
        '''
        self._check_for_cols(self._sompz_cut_cols)

        sompz_cuts = np.where( (self._cat['mag_i'] >= 18.) &
                               (self._cat['mag_i'] <= 23.5) &
                               (self._cat['mag_r'] >= 15.) &
                               (self._cat['mag_r'] <= 26.) &
                               (self._cat['mag_z'] >= 15.) &
                               (self._cat['mag_z'] <= 26.) &
                               ((self._cat['mag_z'] - self._cat['mag_i']) <= 1.5) &
                               ((self._cat['mag_z'] - self._cat['mag_i']) >= -4.) &
                               ((self._cat['mag_r'] - self._cat['mag_i']) <= 4.) &
                               ((self._cat['mag_r'] - self._cat['mag_i']) >= -1.5)
                              )

        self.apply_cut(sompz_cuts)

        if use_match_flag is True:
            match_flag_cut = np.where(self._cat[self._match_flag_col] < 2)
            self.apply_cut(match_flag_cut)

        return

    def _compute_gap_fluxes(self, vb):
        import ngmix

        self._check_for_cols(self._gap_flux_cols)

        # ngmix profile pars to reconstruct light profile
        # Centroid offset is set to 0
        mcal_pars = np.array([
            len(self)*[0.0],
            len(self)*[0.0],
            self._cat['e_1'],
            self._cat['e_2'],
            self._cat['T'],
            self._cat['flux_r'],
            self._cat['flux_i'],
            self._cat['flux_z']
        ]).T

        if vb is True:
            print('Computing Gaussian Aperture fluxes...')
        gap_flux, gap_flags = ngmix.gaussap.get_gaussap_flux(mcal_pars,
                                                             'gauss',
                                                             4.0,
                                                             verbose=vb)

        self._cat['gap_flux_r'] = gap_flux[:,0]
        self._cat['gap_flux_i'] = gap_flux[:,1]
        self._cat['gap_flux_z'] = gap_flux[:,2]

        self._cat['gap_flags_r'] = gap_flags[:,0]
        self._cat['gap_flags_i'] = gap_flags[:,1]
        self._cat['gap_flags_z'] = gap_flags[:,2]

        return

class BalrogMcalCatalog(GoldCatalog, McalCatalog):

    def __init__(self, mcal_file, det_file, mcal_cols=None, det_cols=None,
                 mcal_path='catalog/unsheared', match_type='default', save_all=False,
                 compute_gap_fluxes=False, vb=False):

        self.mcal_file = mcal_file
        self.det_file = det_file
        self.mcal_cols = mcal_cols
        self.det_cols = det_cols
        self.mcal_path = mcal_path
        self.match_type = match_type
        self.save_all = save_all
        self.compute_gap_fluxes = compute_gap_fluxes
        self.vb = vb

        self._set_gold_colname(match_type)
        self._load_catalog()

        return

    def _load_catalog(self):
        if self.vb is True: print('Loading Mcal catalog...')
        mcal = McalCatalog(self.mcal_file, self.mcal_path, cols=self.mcal_cols,
                           compute_gap_fluxes=self.compute_gap_fluxes, vb=self.vb)
        mcal.calc_mags()

        if self.vb is True: print('Loading Detection catalog...')
        det = DetectionCatalog(self.det_file, cols=self.det_cols)

        if self.vb is True: print('Joining catalogs...')
        self._join(mcal.get_cat(), det.get_cat())

        return


    def _join(self, mcal, det):
        self._cat = join(mcal, det, join_type='inner')

        if self.save_all is True:
            self.mcal = mcal
            self.det = det
        else:
            self.mcal = None
            self.det = None

        return


class BalrogDetectionCatalog(DetectionCatalog):

    def __init__(self, filename, cols=None, match_type='default', profile='bdf',
                 has_mags=True, dereddened=True, real=0):

        super(BalrogDetectionCatalog, self).__init__(filename,
                                                     cols=cols,
                                                     match_type=match_type)

        self.profile = profile
        self.prefix = profile + '_'
        self.has_mags = has_mags
        self.dereddened = dereddened
        self.b_indx = dict(zip('griz', range(4)))

        # Could be optional parameters in the future, but enforced for now
        self.true_prefix = 'true_'
        self.meas_prefix = 'meas_'

        p = self.prefix
        if has_mags is True:
            if dereddened is True:
                self.true_mag_colname  = self.true_prefix + p + 'mag_deredden'
        else:
            self.true_mag_colname  = self.true_prefix + p + 'mag'

        if not isinstance(real, int):
            raise TypeError('real must be an int!')
        self.real = real

        return

    def plot_detection_efficiency(self, bands='griz', xlim=[16.0, 30.0], ylim=[0.0, 1.0],
                                  S=8, title=None, cmap='inferno', dim=2, dx=0.1,
                                  vline=None):
            p = self.prefix
            mcol = self.true_mag_colname

            columns = [mcol, 'detected']
            cat = fitsio.read(self.filename, columns=columns)

            N = 2. / (dx) + 1
            mag_min, mag_max = xlim[0], xlim[1]
            bins = np.linspace(mag_min, mag_max, N)

            IN_BIN = {}
            DET = {}
            EFF = {}
            EFF_ERR = {}
            MAG = {}

            for x in [IN_BIN, DET, EFF, EFF_ERR, MAG]:
                for band in bands:
                    x[band] = []

            for band in bands:
                print('Band {}'.format(band))
                bi = self.b_indx[band]
                for i, bstart in enumerate(bins):
                    if i == len(bins)-1:
                        break
                    else:
                        # print 'i, bins[i], bins[i+1] = ',i, bins[i], bins[i+1]

                        in_bin = len(cat[(bins[i]<=cat[mcol][:,bi]) & (cat[mcol][:,bi]<bins[i+1])])
                        det = len(cat[(bins[i]<=cat[mcol][:,bi]) &
                                (cat[mcol][:,bi]<bins[i+1]) &
                                (cat['detected']==1)])
                    try:
                        eff = 100. * det / in_bin
                    except ZeroDivisionError:
                        eff = 0.0

                    mag = np.mean([bins[i], bins[i+1]])

                    try:
                        eff_err = eff * np.sqrt( (1. / det) + (1. / det) )
                    except ZeroDivisionError:
                        eff_err = 0.0

                    IN_BIN[band].append(in_bin)
                    DET[band].append(det)
                    EFF[band].append(eff)
                    EFF_ERR[band].append(eff_err)
                    MAG[band].append(mag)

            bk = 0
            for band in bands:
                # ax = plt.subplot(2, 2, bi+1)

                plt.errorbar(MAG[band], EFF[band], EFF_ERR[band], fmt='o', ms=8, label=band)
                # plt.gca().axhline(0, ls='--', c='k', lw=4)
                # med = np.median(diff[cuts])
                # ax.axhline(med, ls=':', c='w', lw=3, label='Median={:.3f}'.format(med))
                # cb = plt.colorbar(hb, ax=ax)
                # legend = plt.legend(bbox_to_anchor=(0.6, 0.925), bbox_transform=ax.transAxes, fontsize=18)
                # plt.setp(legend.get_texts(), color='w')

                # label = val
                # if self.use_deredden is True:
                #     label += ' (dereddened)'

                # for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                #             ax.get_xticklabels() + ax.get_yticklabels()):
                #     item.set_fontsize(20)

            if title: plt.suptitle(title)
            ax = plt.gca()
            ax.set_xlabel('True {}_mag_deredden'.format(self.profile))
            ax.set_ylabel('Detection Efficiency')
            plt.legend()

            if vline is not None:
                plt.axvline(vline, linewidth=2, ls='--', color='k')

            plt.gcf().set_size_inches(S, S)

            return

class BalrogMatchedCatalog(FitsCatalog):

    def __init__(self, match_file, match_cols=None, match_type='default',
                 vb=False):

        super(BalrogMatchedCatalog, self).__init__(match_file, match_cols)

        # Would normally connect this to a Gold catalog, but the matched
        # catalogs don't have the positional flags yet
        self.match_type = match_type
        self.vb = vb

        return

    def cut_by_bal_id(self, bal_ids):
        '''
        Some sample cuts depend on columns not present in the matched
        catalog. In those cases, it is easiest to make the cut on
        the detection catalog and then filter by bal_id.
        '''

        self.apply_cut(self._cat['bal_id'] in bal_ids)

        return

