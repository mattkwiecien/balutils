from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb
class FitsCatalog(Catalog):
    def _load_catalog(self):
        self._cat = Table(fitsio.read(self.filename, columns=self.cols))
        self.Nobjs = len(self._cat)

        return