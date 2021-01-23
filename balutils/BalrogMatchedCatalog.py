from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb

class BalrogMatchedCatalog(GoldCatalog):

    def __init__(self, match_file, det_file, match_cols=None, det_cols=None,
                 match_type='default', save_all=False, vb=False):

        self.match_file = match_file
        self.det_file = det_file
        self.match_cols = match_cols
        self.det_cols = det_cols
        self.match_type = match_type
        self.save_all = save_all
        self.vb = vb

        self._set_gold_colname(match_type)
        self._load_catalog()

        return

    def _load_catalog(self):
        if self.vb is True:
            print('Loading Matched catalog...')
        match = GoldFitsCatalog(self.match_file,
                                cols=self.match_cols,
                                match_type=self.match_type)

        if self.vb is True:
            print('Loading Detection catalog...')
        det = DetectionCatalog(self.det_file, cols=self.det_cols)

        if self.vb is True:
            print('Joining catalogs...')
        self._join(match.get_cat(), det.get_cat())

        return

    def _join(self, match, det):
        self._cat = join(match, det, join_type='inner')

        if self.save_all is True:
            self.match = match
            self.det = det
        else:
            self.match = None
            self.det = None

        return

    def cut_by_bal_id(self, bal_ids):
        '''
        Some sample cuts depend on columns not present in the matched
        catalog. In those cases, it is easiest to make the cut on
        the detection catalog and then filter by bal_id.
        '''

        self.apply_cut(self._cat['bal_id'] in bal_ids)

        return
