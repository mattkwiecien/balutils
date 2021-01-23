from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb

class BalrogMcalCatalog(GoldCatalog, McalCatalog):

    def __init__(self, mcal_file, det_file, mcal_cols=None, det_cols=None,
                 mcal_path='catalog/unsheared', match_type='default',
                 save_all=False, vb=False):

        self.mcal_file = mcal_file
        self.det_file = det_file
        self.mcal_cols = mcal_cols
        self.det_cols = det_cols
        self.mcal_path = mcal_path
        self.match_type = match_type
        self.save_all = save_all
        self.vb = vb

        self._set_gold_colname(match_type)
        self._load_catalog()

        return

    def _load_catalog(self):
        if self.vb is True:
            print('Loading Mcal catalog...')
        mcal = McalCatalog(self.mcal_file, self.mcal_path, cols=self.mcal_cols)

        if self.vb is True:
            print('Loading Detection catalog...')
        det = DetectionCatalog(self.det_file, cols=self.det_cols)

        if self.vb is True:
            print('Joining catalogs...')
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
