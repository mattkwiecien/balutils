from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb


class McalCatalog(H5Catalog):

    _match_flag_col = 'match_flag_1.5_asec'

    _shape_cut_cols = ['flags',
                       'size_ratio',
                       'snr'
                      ]
    _sompz_cut_cols = ['mag_r',
                       'mag_i',
                       'mag_z',
                       'e_1',
                       'e_2',
                       'T',
                      ]

    _gap_flux_cols = ['e_1',
                      'e_2',
                      'T',
                      'flux_r',
                      'flux_i',
                      'flux_z'
                      ]

    def __init__(self, filename, basepath, cols=None, **kwargs):
        super(McalCatalog, self).__init__(filename, basepath, cols=cols, **kwargs)

        self.calc_mags()

        return

    def calc_mags(self):
        '''
        Mcal catalogs don't automatically come with magnitudes
        '''

        fluxes = [c for c in self._cat.colnames if 'flux' in c.lower()]
        bands = [f[-1] for f in fluxes]

        for b in bands:
            self._cat['mag_{}'.format(b)]= self.flux2mag(self._cat['flux_{}'.format(b)])

        return

    def apply_shape_cuts(self):
        self._check_for_cols(self._shape_cut_cols)
        shape_cuts = np.where( (self._cat['flags'] == 0) &
                               (self._cat['size_ratio'] > 0.5) &
                               (self._cat['snr'] > 10) &
                               (self._cat['snr'] < 1000)
                              )
        self.apply_cut(shape_cuts)

        return

    def apply_sompz_cuts(self, use_match_flag=True):
        '''
        This is the same as "selection2"
        '''
        self._check_for_cols(self._sompz_cut_cols)

        # Mag & color cuts
        color_cuts = np.where( (self._cat['mag_i'] >= 18.) &
                               (self._cat['mag_i'] <= 23.5) &
                               (self._cat['mag_r'] >= 15.) &
                               (self._cat['mag_r'] <= 26.) &
                               (self._cat['mag_z'] >= 15.) &
                               (self._cat['mag_z'] <= 26.) &
                               ((self._cat['mag_z'] - self._cat['mag_i']) <= 1.5) &
                               ((self._cat['mag_z'] - self._cat['mag_i']) >= -4.) &
                               ((self._cat['mag_r'] - self._cat['mag_i']) <= 4.) &
                               ((self._cat['mag_r'] - self._cat['mag_i']) >= -1.5)
                              )

        self.apply_cut(color_cuts)

        # Binary star cut, taken from Alex A.
        highe_cut = np.greater(np.sqrt(np.power(self._cat['e_1'],2.)
                               + np.power(self._cat['e_2'],2)), 0.8)

        c = 22.5
        m = 3.5

        magT_cut = np.log10(self._cat['T']) < (c - self.flux2mag(self._cat['flux_r'])) / m

        binaries = highe_cut * magT_cut

        self.apply_cut(~binaries)

        if use_match_flag is True:
            self._check_for_cols(self._match_flag_col)
            match_flag_cut = np.where(self._cat[self._match_flag_col] < 2)
            self.apply_cut(match_flag_cut)

        return

    def compute_gap_fluxes(self, vb=False):
        import ngmix

        self._check_for_cols(self._gap_flux_cols)

        # ngmix profile pars to reconstruct light profile
        # Centroid offset is set to 0
        mcal_pars = np.array([
            len(self)*[0.0],
            len(self)*[0.0],
            self._cat['e_1'],
            self._cat['e_2'],
            self._cat['T'],
            self._cat['flux_r'],
            self._cat['flux_i'],
            self._cat['flux_z']
        ]).T

        if vb is True:
            print('Computing Gaussian Aperture fluxes...')
        gap_flux, gap_flags = ngmix.gaussap.get_gaussap_flux(mcal_pars,
                                                             'gauss',
                                                             4.0,
                                                             verbose=vb)

        self._cat['gap_flux_r'] = gap_flux[:,0]
        self._cat['gap_flux_i'] = gap_flux[:,1]
        self._cat['gap_flux_z'] = gap_flux[:,2]

        self._cat['gap_flags_r'] = gap_flags[:,0]
        self._cat['gap_flags_i'] = gap_flags[:,1]
        self._cat['gap_flags_z'] = gap_flags[:,2]

        return
