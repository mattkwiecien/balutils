from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb

class DetectionCatalog(FitsCatalog, GoldCatalog):

    def __init__(self, filename, cols=None, match_type='default'):
        super(DetectionCatalog, self).__init__(filename, cols=cols, match_type=match_type)

        self._check_for_duplicates()

        return

    def _check_for_duplicates(self):
        '''
        Balrog stack versions 1.4 and below have a small bug that
        seems to duplicate exactly 1 object, so check for these
        '''
        unq, unq_idx, unq_cnt = np.unique(self._cat['bal_id'],
                                          return_inverse=True,
                                          return_counts=True)
        Nunq = len(unq)
        if Nunq != self.Nobjs:
            Ndups = self.Nobjs - Nunq
            dup_ids = unq[np.where(unq_cnt > 1)]
            print('Warning: Detection catalog has {} duplicate(s)!'.format(Ndups))
            print('Removing the following duplicates from detection catalog:')
            print(dup_ids)

            Nbefore = self.Nobjs
            for did in dup_ids:
                indx = np.where(self._cat['bal_id']==did)[0]

                L = len(indx)
                for i in range(L-1): # keep last one
                    self._cat.remove_row(indx[i])

            self.Nobjs = len(self._cat)
            assert self.Nobjs == (Nbefore - Ndups)

            print('{} duplicates removed, catalog size now {}'.format(Ndups, self.Nobjs))

        return
