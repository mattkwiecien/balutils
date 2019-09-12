import numpy as np
# import fitsio
# import os
import matplotlib.pyplot as plt

import balutils.stacked_catalogs as sc

#df_file   = '/data/des80.a/data/yanny/deepmeds4/run-ugriz-mof02.fits'
df_file   = '/data/des81.a/data/severett/tests/bright_objs/run-ugriz-mof02-tests.fits'
inj_file  = '/data/des80.a/data/yanny/deepmeds4/BALROG_RUN2_DEEP_CAT_FLAG0INVHS1BDFLT254v4.fits'
gold_file = '/data/des81.a/data/severett/tests/bright_objs/Y3_GOLD_2_2_1.0_bobjs_mags.fits'
det_file  = '/data/des41.b/data/severett/Balrog/prod020419/stacked_catalogs/1.4/sof/balrog_detection_catalog_sof_v1.4.fits'
#mcal_file = '/data/des41.b/data/severett/Balrog/prod020419/stacked_catalogs/1.4/mcal/balrog_mcal_stack-y3v02-0-riz-noNB-mcal_v1.4.h5'
mcal_file = '/data/des41.b/data/severett/Balrog/prod020419/stacked_catalogs/1.4/mcal/balrog_mcal_stack-y3v02-0-riz-noNB-mcal_v1.4.h5'
sof_file = '/data/des41.b/data/severett/Balrog/prod020419/stacked_catalogs/1.4/sof/balrog_matched_catalog_sof_v1.4.fits'

def plot_mags(mags, mmin=15, mmax=30, dm=0.25, label='', a=0.5, c=None, s=10, title=None, density=False, show=True):
    bins = np.arange(mmin, mmax+dm, dm)
    bindx = dict(zip('griz', range(4)))

    k = 0
    for b, bi in bindx.items():
        plt.subplot(2,2,bi+1)
        plt.hist(mags[bi], bins=bins, ec='k', label=label, alpha=a, color=c,
                 density=density)
	plt.yscale('log')
        plt.xlabel('{}-band mag'.format(b))
        if density is True:
	    plt.ylabel('Norm. Counts')
	else:
            plt.ylabel('Counts')
	plt.legend()

    if title is not None:
	plt.title(title)

    if show is True:
        plt.gcf().set_size_inches(s,s)
        plt.show()

    return

def main():
    # Code below assumes photo-z gold cuts have already been
    # applied to gold cat to save memory

    # Pick a max where both are complete
    mmax = 22.

    print('Starting DF mag vs Y3 mag...')
    gold_cols = ['SOF_CM_MAG_G', 'SOF_CM_MAG_R', 'SOF_CM_MAG_I', 'SOF_CM_MAG_Z',
                 'SOF_CM_T', 'SOF_CM_T_ERR']
    print('Loading Gold...')
    gold = sc.GoldFitsCatalog(gold_file, cols=gold_cols)
    gcut = np.where(
                     (gold['SOF_CM_MAG_I']<mmax)
                   )
    gold.apply_cut(gcut)

    df_cols = ['bdf_mag']
    print('Loading DF...')
    df = sc.FitsCatalog(df_file, cols=df_cols)
    dcut = np.where(
                     (df['bdf_mag'][:,2]<mmax)
                   )
    df.apply_cut(dcut)

    inj_cols = ['bdf_mag']
    print('Loading Inj...')
    inj = sc.FitsCatalog(inj_file, cols=inj_cols)
    icut = np.where(
                     (inj['bdf_mag'][:,2]<mmax)
                   )
    inj.apply_cut(icut)

    mcal_cols = ['flux_r', 'flux_i', 'flux_z', 'flags', 'psf_T', 'T', 'snr', 'bal_id']
    det_cols  = ['bal_id', 'meas_FLAGS_GOLD', 'flags_footprint', 'flags_badregions', 'flags_foreground']
    bal = sc.BalrogMcalCatalog(mcal_file, det_file, mcal_cols=mcal_cols, det_cols=det_cols, vb=True)
    bal.apply_gold_cuts()
    bcut = np.where(
                     (bal['mag_i']<mmax)
                   )
    bal.apply_cut(bcut)

    #sof_cols = ['meas_cm_mag']
    #print('Loading SOF...')
    #sof = sc.FitsCatalog(sof_file, cols=sof_cols)
    #scut = np.where(
    #                 (sof['cm_mag'][:,2]<mmax)
    #               )
    #sof.apply_cut(scut)

    # DF mag vs Y3 mag
    print('Plotting Gold...')
    plot_mags([gold['SOF_CM_MAG_G'],
               gold['SOF_CM_MAG_R'],
               gold['SOF_CM_MAG_I'],
               gold['SOF_CM_MAG_Z']],
              label='Y3 Gold',
	      c='tab:orange',
              density=True,
              mmax=mmax+1,
              show=False)
    print('Plotting DF...')
    plot_mags([df['bdf_mag'][:,0],
               df['bdf_mag'][:,1],
               df['bdf_mag'][:,2],
               df['bdf_mag'][:,3]],
              label='DF mof02',
	      c='tab:blue',
              density=True,
              mmax=mmax+1,
	      title='DF vs Y3 Gold Sample')

    # Inj vs DF
    print('Plotting Injections...')
    plot_mags([inj['bdf_mag'][:,0],
               inj['bdf_mag'][:,1],
               inj['bdf_mag'][:,2],
               inj['bdf_mag'][:,3]],
              label='Run2 v4 Injections',
              c='tab:green',
              density=True,
              mmax=mmax+1,
              show=False)
    print('Plotting DF...')
    plot_mags([df['bdf_mag'][:,0],
               df['bdf_mag'][:,1],
               df['bdf_mag'][:,2],
               df['bdf_mag'][:,3]],
              label='DF mof02',
              c='tab:blue',
              density=True,
              mmax=mmax+1,
	      title='Injections vs DF')

    # Inj vs WF
    print('Plotting Injections...')
    plot_mags([inj['bdf_mag'][:,0],
               inj['bdf_mag'][:,1],
               inj['bdf_mag'][:,2],
               inj['bdf_mag'][:,3]],
              label='Run2 v4 Injections',
              c='tab:green',
              density=True,
              mmax=mmax+1,
              show=False)
    print('Plotting Gold...')
    plot_mags([gold['SOF_CM_MAG_G'],
               gold['SOF_CM_MAG_R'],
               gold['SOF_CM_MAG_I'],
               gold['SOF_CM_MAG_Z']],
              label='Y3 Gold',
              c='tab:orange',
              density=True,
              mmax=mmax+1,
              title='Injections vs Y3 Gold Sample')

    # Bal vs WF SOF 
    print('Plotting Balrog...')
    plot_mags([np.nan,
               bal['mag_r'],
               bal['mag_i'],
               bal['mag_z']],
              label='Balrog Prerun2',
              c='tab:cyan',
              density=True,
              mmax=mmax+1,
              show=False)
    print('Plotting Gold...')
    plot_mags([gold['SOF_CM_MAG_G'],
               gold['SOF_CM_MAG_R'],
               gold['SOF_CM_MAG_I'],
               gold['SOF_CM_MAG_Z']],
              label='Y3 Gold',
              c='tab:orange',
              density=True,
              mmax=mmax+1,
              title='Balrog Prerun2 (Gold cuts) vs Y3 Gold Sample')
    return

if __name__ == '__main__':
    main()
