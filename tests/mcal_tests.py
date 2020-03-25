import stacked_catalogs as sc
import numpy as np
import matplotlib.pyplot as plt
import pudb

def main():

    # det_file = '/home/spencer/research/balrog/outputs/run2/stacks/1.0/sof/balrog_detection_catalog_sof_run2_v1.0.fits'
    # det_cols = ['bal_id', 'meas_FLAGS_GOLD_SOF_ONLY', 'flags_footprint', 'flags_badregions', 'flags_foreground']
    det_file = '/home/spencer/research/balrog/outputs/y3-merged/stacks/1.1/sof/balrog_detection_catalog_sof_y3-merged_v1.1.fits'
    det_cols = ['bal_id', 'meas_FLAGS_GOLD_SOF_ONLY', 'flags_footprint', 'flags_badregions', 'flags_foreground', 'match_flag_1.5_asec']

    bal_mcal_file = '/home/spencer/research/balrog/outputs/run2/stacks/1.0/mcal/balrog_mcal_stack-y3v02-0-riz-noNB-mcal_v1.0.h5'
    bal_mcal_cols = ['flux_i', 'flux_r', 'flux_z', 'flags', 'mcal_psf_T', 'T', 'snr', 'bal_id',
                     'flux_err_r', 'flux_err_i', 'flux_err_z', 'e_1', 'e_2', 'size_ratio',
                     'R11', 'R22']

    print('First test - unsheared only')
    print(sc.McalGoldCatalog.mro())
    bal_mcal = sc.BalrogMcalCatalogs(bal_mcal_file, det_file, mcal_cols=bal_mcal_cols, det_cols=det_cols, match_type='sof_only', stypes=['unsheared'], vb=True)

    print('Applying cuts...')
    bal_mcal.apply_gold_cuts()
    bal_mcal.apply_shape_cuts()
    bal_mcal.apply_sompz_cuts()

    print('Accessing unsheared')
    unsheared = bal_mcal['unsheared']
    print(unsheared.get_cat().colnames)
    print('Accessing unsheared/snr')
    pudb.set_trace()
    snr = bal_mcal['unsheared']['snr']
    plt.hist(snr, bins=50, ec='k')
    plt.xlabel('unsheared/snr')
    plt.show()

    print('Applying singular cut')
    bal_mcal.apply_cut(np.where(snr > 50))

    print('Second test - all shear types')
    bal_mcal = sc.BalrogMcalCatalogs(bal_mcal_file, det_file, mcal_cols=bal_mcal_cols, det_cols=det_cols, match_type='sof_only', vb=True)

    print('Applying cuts...')
    bal_mcal.apply_gold_cuts()
    bal_mcal.apply_shape_cuts()
    bal_mcal.apply_sompz_cuts()

    print('Accessing unsheared')
    unsheared = bal_mcal['unsheared']
    print(unsheared.get_cat().colnames)
    print('Accessing unsheared/snr')
    pudb.set_trace()
    snr = bal_mcal['unsheared']['snr']
    plt.hist(snr, bins=50, ec='k')
    plt.xlabel('unsheared/snr')
    plt.show()

    for stype in bal_mcal.stypes:
        b = np.arange(0, 5+.1, .1)
        plt.hist(bal_mcal[stype]['size_ratio'], label=stype, histtype='step', bins=b)
    plt.xlabel('size_ratio')
    plt.ylabel('counts')
    plt.yscale('log')
    plt.legend()
    plt.gcf().set_size_inches(9,6)
    plt.plot()

    return

if __name__ == '__main__':
    main()
