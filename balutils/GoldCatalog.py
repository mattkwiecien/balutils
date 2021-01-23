from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb

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

    def __init__(self, filename, cols=None, match_type='default', **kwargs):
        super(GoldCatalog, self).__init__(filename, cols=cols, **kwargs)

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
