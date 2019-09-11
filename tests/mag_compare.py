import numpy as np
# import fitsio
# import os
import matplotlib.pyplot as plt

import balutils.stacked_catalogs as sc

df_file   = '/data/des80.a/data/yanny/deepmeds4/run-ugriz-mof02.fits'
gold_file = '/data/des81.a/data/severett/tests/bright_objs/Y3_GOLD_2_2_0.02_mags.fits'
det_file  = '/data/des41.b/data/severett/Balrog/prod020419/stacked_catalogs/1.4/sof/balrog_detection_catalog_sof_v1.4.fits'
mcal_file = '/data/des41.b/data/severett/Balrog/prod020419/stacked_catalogs/1.4/mcal/balrog_mcal_stack-y3v02-0-riz-noNB-mcal_v1.4.h5'

def plot_mags(mags, mmin=15, mmax=30, dm=0.25, label='', a=0.5, s=10, show=True):
    bins = np.arange(mmin, mmax+dm, dm)
    bindx = dict(zip('griz', range(4)))

    k = 0
    for b, bi in bindx.items():
        plt.subplot(2,2,bi)
        plt.hist(mags[bi], bins=bins, ec='k', label=label, alpha=a)
        plt.xlabel('{}-band mag'.format(b))

    if show is True:
        plt.legend()
        plt.gcf().set_size_inches(s,s)
        plt.show()

    return

def main():
    # Code below assumes photo-z gold cuts have already been
    # applied to gold cat to save memory

    # DF mag vs Y3 mag
    print('Starting DF mag vs Y3 mag...')
    gold_cols = ['SOF_CM_MAG_G', 'SOF_CM_MAG_R', 'SOF_CM_MAG_I', 'SOF_CM_MAG_Z']
    gold = sc.GoldFitsCatalog(gold_file, cols=gold_cols)

    df_cols = ['bdf_mag']
    df = sc.FitsCatalog(df_file, cols=df_cols)

    plot_mags([gold['SOF_CM_MAG_G'],
               gold['SOF_CM_MAG_R'],
               gold['SOF_CM_MAG_I'],
               gold['SOF_CM_MAG_Z']],
              label='Y3 Gold',
              show=False)
    plot_mags([df[df['bdf_mag'][0]],
               df[df['bdf_mag'][1]],
               df[df['bdf_mag'][2]],
               df[df['bdf_mag'][3]]],
              label='DF mof02')

    ## More plots...

    return

if __name__ == '__main__':
    main()
