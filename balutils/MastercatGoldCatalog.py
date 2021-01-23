from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb

# TODO: In progress!!!
class MastercatGoldCatalog(Catalog):
    def __init__(filename='/global/project/projectdirs/des/www/y3_cats/Y3_mastercat_12_3_19.h5',
                 cols=None, bands='griz', vb=False):

        super(MastercatGoldCatalog, self).__init__(filename, cols)

        self.bands = bands
        self.vb = vb

        self._cat = Table()

        h5 = h5py.File(filename, mode='r')
        select = h5['index/select'][:]

        for b in bands:
            if self.vb is True:
                print('Starting band {}'.format(b))
            gld_mcal['flux_{}'.format(b)] = gld['catalog/gold...fdk/flux_{}'.format(b)][:][select]
            gld_mcal['flux_err_{}'.format(b)] = gld['catalog/metacal/unsheared/flux_err_{}'.format(b)][:][select]
            gld_mcal['mag_{}'.format(b)] = flux2mag(gld_mcal['flux_{}'.format(b)])

        colnames = ['T']

        for col in colnames:
            gld_mcal[col] = gld['catalog/metacal/unsheared/'+col][:][select]
