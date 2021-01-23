from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb

class McalGoldCatalog(McalCatalog, GoldCatalog):
    pass
    # def __init__(self, filename, basepath, cols=None):
    #     super(McalGoldCatalog, self).__init__(filename, basepath, cols=cols)

    #     return
