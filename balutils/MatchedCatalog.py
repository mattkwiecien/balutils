from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb


class MatchedCatalog(FitsCatalog):

    def __init__(self, match_file, match_cols=None, match_type='default',
                 vb=False):

        super(BalrogMatchedCatalog, self).__init__(match_file, match_cols)

        # Would normally connect this to a Gold catalog, but the matched
        # catalogs don't have the positional flags yet
        self.match_type = match_type
        self.vb = vb

        return

    def cut_by_bal_id(self, bal_ids):
        '''
        Some sample cuts depend on columns not present in the matched
        catalog. In those cases, it is easiest to make the cut on
        the detection catalog and then filter by bal_id.
        '''

        self.apply_cut(self._cat['bal_id'] in bal_ids)

        return
