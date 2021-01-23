from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb

class GoldFitsCatalog(FitsCatalog, GoldCatalog):
    def __init__(self, filename, cols=None, match_type='default'):
        super(GoldFitsCatalog, self).__init__(filename, cols=cols, match_type=match_type)

        return
