from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb

class BalrogMcalCatalogs(BalrogMcalCatalog):

    _valid_shear_types = ['unsheared',
                          'sheared_1m',
                          'sheared_1p',
                          'sheared_2m',
                          'sheared_2p'
                          ]

    # Some cols are only saved to unsheared - often we want them
    # in all shear types for cuts
    _mcal_flats = {'bal_id',
                   'coadd_object_id',
                   'flags',
                   'mask_frac',
                   'ra',
                   'dec',
                   'psf_T',
                   'mcal_psf_T',
                   'e_1',
                   'e_2',
                   'R11',
                   'R12',
                   'R21',
                   'R22',
                   }

    def __init__(self, mcal_file, det_file, mcal_cols=None, det_cols=None,
                 stypes='all', match_type='default', save_all=False, vb=False):

        if stypes == 'all':
            self.stypes = self._valid_shear_types
        else:
            self.stypes = stypes
        if 'unsheared' in self.stypes:
            self.has_unsheared = True
        else:
            self.has_unsheared = False

        super(BalrogMcalCatalogs, self).__init__(mcal_file, det_file, mcal_cols=mcal_cols,
                                                 det_cols=det_cols, mcal_path=None,
                                                 match_type=match_type, save_all=save_all,
                                                 vb=vb)
        return

    def _load_catalog(self):
        # Different from other classes, we need the underlying cat to be a dict
        self._cat = {}
        self.mcal = {}
        self.det = {}

        if self.vb is True:
            print('Loading Detection catalog...')
        det = DetectionCatalog(self.det_file, cols=self.det_cols)

        flats = {}
        mcal_cols = self.mcal_cols[:]
        for stype in self.stypes:
            if self.vb is True:
                print('Loading {} Mcal catalog...'.format(stype))

            if stype != 'unsheared':
                for c in self._mcal_flats:
                    if c in mcal_cols:
                        mcal_cols.remove(c)

            path = 'catalog/'+stype+'/'
            mcal = McalGoldCatalog(self.mcal_file, path, cols=mcal_cols, match_type=self.match_type)

            if stype == 'unsheared':
                for col in self._mcal_flats:
                    if col in self.mcal_cols:
                        flats[col] = mcal[col]
            else:
                for col in self._mcal_flats:
                    if col in self.mcal_cols:
                        mcal[col] = flats[col]

            if self.vb is True:
                print('Joining catalogs...')
            self._join(stype, mcal, det)

        return

    def _join(self, stype, mcal, det):
        mcal._cat = join(mcal.get_cat(), det.get_cat(), join_type='inner')
        self._cat[stype] = mcal

        if self.save_all is True:
            self.mcal[stype] = mcal
            self.det[stype] = det
        else:
            self.mcal[stype] = None
            self.det[stype] = None

        return

    def apply_gold_cuts(self):
        for stype in self.stypes:
            if self.vb is True:
                print('Applying gold cuts to {}'.format(stype))
            self._cat[stype].apply_gold_cuts()

        return

    def apply_shape_cuts(self):
        for stype in self.stypes:
            if self.vb is True:
                print('Applying shape cuts to {}'.format(stype))
            self._cat[stype].apply_shape_cuts()

        return

    def apply_sompz_cuts(self, use_match_flag=True):
        for stype in self.stypes:
            if self.vb is True:
                print('Applying sompz cuts to {}'.format(stype))
            self._cat[stype].apply_sompz_cuts(use_match_flag=use_match_flag)

        return

    def apply_cut(self, cut, stype='all'):
        if stype == 'all':
            for s in self.stypes:
                if self.vb is True:
                    print('Applying cut to {}'.format(stype))
                self._cat[s].apply_cut(cut)
        else:
            elf._cat[stype].apply_cut(cut)

        return
