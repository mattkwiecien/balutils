from abc import ABCMeta, abstractmethod
import os
import fitsio
import h5py
import numpy as np
from astropy.table import Table, vstack, join
import matplotlib.pyplot as plt
import pudb


class BalrogDetectionCatalog(DetectionCatalog):

    def __init__(self, filename, cols=None, match_type='default', profile='bdf',
                 has_mags=True, dereddened=True, real=0):

        super(BalrogDetectionCatalog, self).__init__(filename,
                                                     cols=cols,
                                                     match_type=match_type)

        self.profile = profile
        self.prefix = profile + '_'
        self.has_mags = has_mags
        self.dereddened = dereddened
        self.b_indx = dict(zip('griz', range(4)))

        # Could be optional parameters in the future, but enforced for now
        self.true_prefix = 'true_'
        self.meas_prefix = 'meas_'

        p = self.prefix
        if has_mags is True:
            if dereddened is True:
                self.true_mag_colname  = self.true_prefix + p + 'mag_deredden'
        else:
            self.true_mag_colname  = self.true_prefix + p + 'mag'

        if not isinstance(real, int):
            raise TypeError('real must be an int!')
        self.real = real

        return

    def plot_detection_efficiency(self, bands='griz', xlim=[16.0, 30.0], ylim=[0.0, 1.0],
                                  S=8, title=None, cmap='inferno', dim=2, dx=0.1,
                                  vline=None):
            p = self.prefix
            mcol = self.true_mag_colname

            columns = [mcol, 'detected']
            cat = fitsio.read(self.filename, columns=columns)

            N = 2. / (dx) + 1
            mag_min, mag_max = xlim[0], xlim[1]
            bins = np.linspace(mag_min, mag_max, N)

            IN_BIN = {}
            DET = {}
            EFF = {}
            EFF_ERR = {}
            MAG = {}

            for x in [IN_BIN, DET, EFF, EFF_ERR, MAG]:
                for band in bands:
                    x[band] = []

            for band in bands:
                print('Band {}'.format(band))
                bi = self.b_indx[band]
                for i, bstart in enumerate(bins):
                    if i == len(bins)-1:
                        break
                    else:
                        # print 'i, bins[i], bins[i+1] = ',i, bins[i], bins[i+1]

                        in_bin = len(cat[(bins[i]<=cat[mcol][:,bi]) & (cat[mcol][:,bi]<bins[i+1])])
                        det = len(cat[(bins[i]<=cat[mcol][:,bi]) &
                                (cat[mcol][:,bi]<bins[i+1]) &
                                (cat['detected']==1)])
                    try:
                        eff = 100. * det / in_bin
                    except ZeroDivisionError:
                        eff = 0.0

                    mag = np.mean([bins[i], bins[i+1]])

                    try:
                        eff_err = eff * np.sqrt( (1. / det) + (1. / det) )
                    except ZeroDivisionError:
                        eff_err = 0.0

                    IN_BIN[band].append(in_bin)
                    DET[band].append(det)
                    EFF[band].append(eff)
                    EFF_ERR[band].append(eff_err)
                    MAG[band].append(mag)

            bk = 0
            for band in bands:
                # ax = plt.subplot(2, 2, bi+1)

                plt.errorbar(MAG[band], EFF[band], EFF_ERR[band], fmt='o', ms=8, label=band)
                # plt.gca().axhline(0, ls='--', c='k', lw=4)
                # med = np.median(diff[cuts])
                # ax.axhline(med, ls=':', c='w', lw=3, label='Median={:.3f}'.format(med))
                # cb = plt.colorbar(hb, ax=ax)
                # legend = plt.legend(bbox_to_anchor=(0.6, 0.925), bbox_transform=ax.transAxes, fontsize=18)
                # plt.setp(legend.get_texts(), color='w')

                # label = val
                # if self.use_deredden is True:
                #     label += ' (dereddened)'

                # for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                #             ax.get_xticklabels() + ax.get_yticklabels()):
                #     item.set_fontsize(20)

            if title: plt.suptitle(title)
            ax = plt.gca()
            ax.set_xlabel('True {}_mag_deredden'.format(self.profile))
            ax.set_ylabel('Detection Efficiency')
            plt.legend()

            if vline is not None:
                plt.axvline(vline, linewidth=2, ls='--', color='k')

            plt.gcf().set_size_inches(S, S)

            return
