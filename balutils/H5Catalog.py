from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb


class H5Catalog(Catalog):

    def __init__(self, filename, basepath, cols=None, **kwargs):
        self.basepath = basepath
        super(H5Catalog, self).__init__(filename, cols=cols, **kwargs)

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
