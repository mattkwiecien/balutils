from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb

class Catalog(object):

    def __init__(self, filename, cols=None):
        self.filename = filename
        self.cols = cols
        self._cat = None

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
        if not isinstance(cols, list):
            cols = [cols]
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
