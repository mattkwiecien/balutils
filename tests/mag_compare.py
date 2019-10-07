import numpy as np
# import fitsio
# import os
import matplotlib.pyplot as plt

import balutils.stacked_catalogs as sc

# Files for Run1:
#run_name = 'Run1'
#df_file  = '/data/des81.a/data/severett/tests/bright_objs/Y3A2_MISC_MOF_V1_COADD_TRUTH_MAGS.fits'
#inj_file  = '/data/des81.a/data/severett/tests/bright_objs/BALROG_RUN1_DEEP_CAT_FLAG0INVHS1SNGTM3DET254.fits'
#det_file  = '/data/des41.b/data/severett/Balrog/prod020419/stacked_catalogs/1.4/sof/balrog_detection_catalog_sof_v1.4.fits'
#gold_file = '/data/des81.a/data/severett/tests/bright_objs/Y3_GOLD_2_2_1.0_bobjs_mags.fits'
#mcal_file = '/data/des41.b/data/severett/Balrog/prod020419/stacked_catalogs/1.4/mcal/balrog_mcal_stack-y3v02-0-riz-noNB-mcal_v1.4.h5'
#sof_file = '/data/des41.b/data/severett/Balrog/prod020419/stacked_catalogs/1.4/sof/balrog_matched_catalog_sof_v1.4.fits'

# Files for Prerun2:
run_name = 'Prerun2'
df_file   = '/data/des80.a/data/yanny/deepmeds4/run-ugriz-mof02.fits'
df_file   = '/data/des81.a/data/severett/tests/bright_objs/run-ugriz-mof02-tests.fits'
inj_file  = '/data/des80.a/data/yanny/deepmeds4/BALROG_RUN2_DEEP_CAT_FLAG0INVHS1BDFLT254v4.fits'
gold_file = '/data/des81.a/data/severett/tests/bright_objs/Y3_GOLD_2_2_1.0_bobjs_mags.fits'

# det_file  = '/data/des81.a/data/severett/matched_stacks/prerun2/sof/balrog_detection_catalog_sof_v1.0.fits'
# sof_file  = '/data/des81.a/data/severett/matched_stacks/prerun2/sof/balrog_matched_catalog_sof_v1.0.fits'
# mcal_file = '/data/des81.a/data/severett/matched_stacks/prerun2/mcal/balrog_mcal_stack-y3v02-0-riz-mcal_v1.0.h5'

det_file = '/data/des81.a/data/severett/matched_stacks/run2/sof/balrog_detection_catalog_sof_run2_v1.0.fits'
sof_file  = '/data/des81.a/data/severett/matched_stacks/run2/sof/balrog_matched_catalog_sof_v1.0.fits'
mcal_file = '/data/des81.a/data/severett/matched_stacks/run2/mcal/balrog_mcal_stack-y3v02-0-riz-noNB-mcal_v1.0.h5'

if '1' in run_name:
    match_type = 'default'
elif '2' in run_name:
    match_type = 'sof_only'

make_plots = {
	'df_vs_wf':True,
	'inj_vs_df':True,
	'inj_vs_wf':False,
	'sof_vs_wf':False,
	'mcal_vs_wf':False,
	'sof_vs_wf_s2n':False
}

def plot_mags(mags, mmin=15, mmax=30, dm=0.25, label='', a=0.5, c=None, s=10, title=None, density=False, show=True):
    bins = np.arange(mmin, mmax+dm, dm)
    bindx = dict(zip('griz', range(4)))

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

def plot_dist(x, xmin=0, xmax=100, dx=1, xlabel='', label='', a=0.5, c=None, s=10, title=None, density=False, show=True):
    bins = np.arange(xmin, xmax+dx, dx)

    plt.hist(x, bins=bins, ec='k', label=label, alpha=a, color=c, density=density)
    plt.yscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(xlabel)

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

    if '1' in run_name:
        df_cols = ['CM_MAG_G', 'CM_MAG_R', 'CM_MAG_I', 'CM_MAG_Z']
        print('Loading DF...')
        df = sc.FitsCatalog(df_file, cols=df_cols)
        dcut = np.where(
                         (df['CM_MAG_I']<mmax)
                       )
        df.apply_cut(dcut)
    elif '2' in run_name:
        df_cols = ['bdf_mag', 'bdf_mask_flags']
        print('Loading DF...')
        df = sc.FitsCatalog(df_file, cols=df_cols)
        dcut = np.where(
                         (df['bdf_mag'][:,2]<mmax) &
                         (df['bdf_mask_flags']!=0)
                       )
        df.apply_cut(dcut)

    if '1' in run_name:
        inj_cols = ['cm_mag']
        print('Loading Inj...')
        inj = sc.FitsCatalog(inj_file, cols=inj_cols)
        icut = np.where(
                         (inj['cm_mag'][:,2]<mmax)
                       )
        inj.apply_cut(icut)
    elif '2' in run_name:
        inj_cols = ['bdf_mag']
        print('Loading Inj...')
        inj = sc.FitsCatalog(inj_file, cols=inj_cols)
        icut = np.where(
                         (inj['bdf_mag'][:,2]<mmax)
                       )
        inj.apply_cut(icut)

    mcal_cols = ['flux_r', 'flux_i', 'flux_z', 'flags', 'psf_T', 'T', 'snr', 'bal_id']
    det_cols  = ['bal_id', 'meas_FLAGS_GOLD', 'flags_footprint', 'flags_badregions', 'flags_foreground']
    mcal = sc.BalrogMcalCatalog(mcal_file, det_file, mcal_cols=mcal_cols, det_cols=det_cols,
                               match_type=match_type, vb=True)
    mcal.apply_gold_cuts()
    bcut = np.where(
                     (mcal['mag_i']<mmax)
                   )
    mcal.apply_cut(bcut)

    sof_cols = ['bal_id', 'meas_cm_mag', 'meas_cm_s2n_r']
    print('Loading SOF...')
    sof = sc.BalrogMatchedCatalog(sof_file, det_file, match_cols=sof_cols, det_cols=det_cols, match_type=match_type, vb=True)
    scut = np.where(
                     (sof['meas_cm_mag'][:,2]<mmax)
                   )
    sof.apply_cut(scut)

    # DF mag vs Y3 mag
    if make_plots['df_vs_wf'] is True:
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
	if '1' in run_name:
            plot_mags([df['CM_MAG_G'],
                       df['CM_MAG_R'],
                       df['CM_MAG_I'],
                       df['CM_MAG_Z']],
                      label='DF COADD_TRUTH',
                      c='tab:blue',
                      density=True,
                      mmax=mmax+1,
                      title='DF vs Y3 Gold Sample')
	elif '2' in run_name:
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
    if make_plots['inj_vs_df'] is True:
        print('Plotting Injections...')
        if '1' in run_name:
            plot_mags([inj['cm_mag'][:,0],
                       inj['cm_mag'][:,1],
                       inj['cm_mag'][:,2],
                       inj['cm_mag'][:,3]],
                      label='Run1 Injections',
                      c='tab:green',
                      density=True,
                      mmax=mmax+1,
                      show=False)
        elif '2' in run_name:
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
	if '1' in run_name:
            plot_mags([df['CM_MAG_G'],
                       df['CM_MAG_R'],
                       df['CM_MAG_I'],
                       df['CM_MAG_Z']],
                      label='DF COADD_TRUTH',
                      c='tab:blue',
                      density=True,
                      mmax=mmax+1,
                      title='DF vs Y3 Gold Sample')
	elif '2' in run_name:
            plot_mags([df['bdf_mag'][:,0],
                       df['bdf_mag'][:,1],
                       df['bdf_mag'][:,2],
                       df['bdf_mag'][:,3]],
                      label='DF mof02',
                      c='tab:blue',
                      density=True,
                      mmax=mmax+1,
                      title='DF vs Y3 Gold Sample')

    # Inj vs WF
    if make_plots['inj_vs_wf'] is True:
        print('Plotting Injections...')
        if '1' in run_name:
            plot_mags([inj['cm_mag'][:,0],
                       inj['cm_mag'][:,1],
                       inj['cm_mag'][:,2],
                       inj['cm_mag'][:,3]],
                      label='Run1 Injections',
                      c='tab:green',
                      density=True,
                      mmax=mmax+1,
                      show=False)
        elif '2' in run_name:
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

    # Bal SOF (all) vs WF SOF
    if make_plots['sof_vs_wf'] is True:
        print('Plotting Balrog SOF...')
        plot_mags([sof['meas_cm_mag'][:,0],
                   sof['meas_cm_mag'][:,1],
                   sof['meas_cm_mag'][:,2],
                   sof['meas_cm_mag'][:,3]],
                  label='Balrog SOF {}'.format(run_name),
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
                  title='Balrog SOF {} (all) vs Y3 Gold Sample'.format(run_name))

    # Bal SOF (gold) vs WF SOF
    sof.apply_gold_cuts()
    if make_plots['sof_vs_wf'] is True:
        print('Plotting Balrog SOF...')
        plot_mags([sof['meas_cm_mag'][:,0],
                   sof['meas_cm_mag'][:,1],
                   sof['meas_cm_mag'][:,2],
                   sof['meas_cm_mag'][:,3]],
                  label='Balrog SOF {}'.format(run_name),
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
                  title='Balrog SOF {} (Gold cuts) vs Y3 Gold Sample'.format(run_name)
		)

    # Bal SOF (gold) vs WF SOF S2N
    if make_plots['sof_vs_wf_s2n'] is True:
        print('Plotting Balrog SOF...')
        plot_dist(sof['meas_cm_s2n_r'],
                  label='Balrog SOF {}'.format(run_name),
                  xlabel='cm_s2n_r',
                  c='tab:cyan',
                  density=True,
                  xmax=150,
                  show=False)
        print('Plotting Gold...')
        plot_dist([gold['SOF_CM_FLUX_S2N_G'],
                   gold['SOF_CM_FLUX_S2N_R'],
                   gold['SOF_CM_FLUX_S2N_I'],
                   gold['SOF_CM_FLUX_S2N_Z']],
                  label='Y3 Gold',
                  xlabel='cm_s2n_r',
                  c='tab:orange',
                  density=True,
                  mmax=mmax+1,
                  title='Balrog SOF {} (Gold cuts) vs Y3 Gold Sample'.format(run_name))

    # Bal mcal vs WF SOF 
    if make_plots['mcal_vs_wf'] is True:
        print('Plotting Balrog mcal...')
        plot_mags([np.nan,
                   mcal['mag_r'],
                   mcal['mag_i'],
                   mcal['mag_z']],
                  label='Balrog {}'.format(run_name),
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
                  title='Balrog Mcal {} (Gold cuts) vs Y3 Gold Sample'.format(run_name))

    return

if __name__ == '__main__':
    main()
