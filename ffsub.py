'''
Script that utilizes H41a continuum-subracted cubes to estimate the free-free continuum
at B3 and B6, and subtracts this contribution from "cleanest" ALMA-IMF continuum maps.
'''


# ff_tols.py needs to be in folder or in path
from ff_tools import *
import os
from astropy.io import fits

################### USER DEFINED ############################
data_path = '/home/roberto/ALMA_IMF/freefree/H41a/'


detections = {
    'G008.67':['G008.67_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G008.67_B3_uid___A001_X1296_X1c1_continuum_merged_12M_robust0_selfcal5_finaliter.image.tt0.fits',
    'G008.67_B6_uid___A001_X1296_X1b7_continuum_merged_12M_robust0_selfcal5_finaliter.image.tt0.fits',
    [6.8,84.0],[5.0,5.0],[7000.,0.08]],
    'G010.62':['G010.62_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G010.62_B3_uid___A001_X1296_X1e5_continuum_merged_12M_robust0_selfcal9_finaliter.image.tt0.fits',
    'G010.62_B6_uid___A001_X1296_X1db_continuum_merged_12M_robust0_selfcal5_finaliter.image.tt0.fits',
    [-39.1,41.9],[5.0,5.0],[7000.,0.071]],
    'G012.80':['G012.80_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G012.80_B3_uid___A001_X1296_X1fb_continuum_merged_12M_robust0_selfcal7_finaliter.image.tt0.fits',
    'G012.80_B6_uid___A001_X1296_X1ef_continuum_merged_12M_robust0_selfcal6_finaliter.image.tt0.fits',
    [-5.8,84.4],[5.0,5.0],[7000.,0.065]],
    'G327.29':['G327.29_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G327.29_B3_uid___A001_X1296_X17d_continuum_merged_12M_robust0_selfcal2_finaliter.image.tt0.fits',
    'G327.29_B6_uid___A001_X1296_X175_continuum_merged_12M_robust0_selfcal5_finaliter.image.tt0.fits',
    [-63.4,-19.2],[5.0,5.0],[7000.,0.08]],
    'G333.60':['G333.60_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G333.60_B3_uid___A001_X1296_X1a3_continuum_merged_12M_robust0_selfcal6_finaliter.image.tt0.fits',
    'G333.60_B6_uid___A001_X1296_X19b_continuum_merged_12M_robust0_selfcal6_finaliter.image.tt0.fits',
    [-103.3,5.2],[5.0,5.0],[7000.,0.068]],
    'G337.92':['G337.92_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G337.92_B3_uid___A001_X1296_X147_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    'G337.92_B6_uid___A001_X1296_X13b_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    [-52.8,-8.7],[5.0,5.0],[7000.,0.08]],
    'G351.77':['G351.77_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G351.77_B3_uid___A001_X1296_X209_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    'G351.77_B6_uid___A001_X1296_X201_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    [-29.0,20.6],[5.0,5.0],[7000.,0.08]],
    'G353.41':['G353.41_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G353.41_B3_uid___A001_X1296_X1d5_continuum_merged_12M_robust0_selfcal6_finaliter.image.tt0.fits',
    'G353.41_B6_uid___A001_X1296_X1c9_continuum_merged_12M_robust0_selfcal6_finaliter.image.tt0.fits',
    [-46.7,19.5],[5.0,5.0],[7000.,0.028]],
    'W43-MM2':['W43-MM2_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'W43-MM2_B3_uid___A001_X1296_X11b_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    'W43-MM2_B6_uid___A001_X1296_X113_continuum_merged_12M_robust0_selfcal5_finaliter.image.tt0.fits',
    [61.9,122.6],[5.0,5.0],[7000.,0.258]],
    'W43-MM3':['W43-MM3_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'W43-MM3_B3_uid___A001_X1296_X12f_continuum_merged_12M_robust0_selfcal5_finaliter.image.tt0.fits',
    'W43-MM3_B6_uid___A001_X1296_X129_continuum_merged_12M_robust0_selfcal5_finaliter.image.tt0.fits',
    [61.1,120.0],[5.0,5.0],[7000.,0.077]],
    'W51-E':['W51-E_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'W51-E_B3_uid___A001_X1296_X10b_continuum_merged_12M_robust0_selfcal7_finaliter.image.tt0.fits',
    'W51-E_B6_uid___A001_X1296_X213_continuum_merged_12M_robust0_selfcal7_finaliter.image.tt0.fits',
    [22.8,100],[5.0,5.0],[7000.,0.094]],
    'W51-IRS2':['W51-IRS2_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'W51-IRS2_B3_uid___A001_X1296_X18f_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    'W51-IRS2_B6_uid___A001_X1296_X187_continuum_merged_12M_robust0_selfcal9_finaliter.image.tt0.fits',
    [8.7,93.4],[5.0,5.0],[7000.,0.104]]
    }

nondetections = {
    'G328.25':['G328.25_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G328.25_B3_uid___A001_X1296_X16d_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    'G328.25_B6_uid___A001_X1296_X163_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    [-53.0,-33.0],[5.0,5.0],[7000.,0.08]],
    'G338.93':['G338.93_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'G338.93_B3_uid___A001_X1296_X159_continuum_merged_12M_robust0_selfcal3_finaliter.image.tt0.fits',
    'G338.93_B6_uid___A001_X1296_X14f_continuum_merged_12M_robust0_selfcal6_finaliter.image.tt0.fits',
    [-72.0,-52.0],[5.0,5.0],[7000.,0.08]],
    'W43-MM1':['W43-MM1_B3_spw1_12M_h41a.JvM.image.contsub.fits',
    'W43-MM1_B3_uid___A001_X1296_X1af_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    'W43-MM1_B6_uid___A002_X996c88_X87_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.fits',
    [80.3,115.3],[5.0,5.0],[7000.,0.08]],
}


#Single source for tests
#data_files = {
#    'G008.67':['G008.67_B3_spw1_12M_h41a.JvM.image.contsub.fits',
#    'G008.67_B3_uid___A001_X1296_X1c1_continuum_merged_12M_robust0_selfcal5_finaliter.image.tt0.fits',
#    'G008.67_B6_uid___A001_X1296_X1b7_continuum_merged_12M_robust0_selfcal5_finaliter.image.tt0.fits',
#    [6.8,84.0],[5.0,5.0]]}

###########################################################################################################


for field in detections:
    print('\n**********')
    print('Free-free estimation for '+field)
    print('**********')

    # Correct files for RL cube and continuum images. Extract frequency from continuum images.
    input_cube = os.path.join(data_path,detections[field][0])
    cont_image_B3 = os.path.join(data_path,detections[field][1])
    freq_B3 = fits.getheader(cont_image_B3)['CRVAL3']
    cont_image_B6 = os.path.join(data_path,detections[field][2])
    freq_B6 = fits.getheader(cont_image_B6)['CRVAL3']

    # Moment 0 of RL image with custom velocity range
    mom0_image = mom0(input_cube,vel_range=detections[field][3])

    ##################################################
    # For estimating and subtracting free-free at 3mm:
    print('Estimating B3:')
    ##################################################

    # Free-free estimation image from RL mom0, and error image:
    mmff_image_B3, err_image_B3 = ff_rl(mom0_image,mask_thr=detections[field][4][0], freq_scale=True, new_freq=freq_B3/1e9, 
    fill_zeros=True, error_map=True, Te=detections[field][5][0], NHe_NH=detections[field][5][1])

    # Convolution of continuum (smaller beam) to mom0 image beam:
    conv_cont_B3 = convol(cont_image_B3,mom0_image)

    # Regridding image of ff estimation and error image to pixel geometry of convolved continuum image:
    regr_mmff_B3 = regrid_tocont(mmff_image_B3,conv_cont_B3)
    regr_err_B3 = regrid_tocont(err_image_B3,conv_cont_B3)

    # Subtraction of convolved continuum - regridded ff estimate
    ffsub_image_B3 = subtr_ff(conv_cont_B3,regr_mmff_B3)

    # PB correction of ff estimate, error map, and subtracted map
    pb_im_B3 = cont_image_B3.replace('.fits','.pb.fits')
    pbcor_mmff_B3 = pbcor(regr_mmff_B3,pb_im_B3,squeeze=True)
    pbcor_err_B3 = pbcor(regr_err_B3,pb_im_B3,squeeze=True)
    pbcor_ffsub_B3 = pbcor(ffsub_image_B3,pb_im_B3)

    ##################################################
    # For estimating and subtracting free-free at 1mm:
    print('Estimating B6:')
    ##################################################

    # Free-free estimation image from RL mom0, and error image:
    mmff_image_B6, err_image_B6 = ff_rl(mom0_image,mask_thr=detections[field][4][1], freq_scale=True, new_freq=freq_B6/1e9, 
    fill_zeros=True, error_map=True, Te=detections[field][5][0], NHe_NH=detections[field][5][1])

    # Convolution of continuum (smaller beam) to mom0 image beam:
    conv_cont_B6 = convol(cont_image_B6,mom0_image)

    # Regridding image of ff estimation and error image to pixel geometry of convolved continuum image:
    regr_mmff_B6 = regrid_tocont(mmff_image_B6,conv_cont_B6)
    regr_err_B6 = regrid_tocont(err_image_B6,conv_cont_B6)

    # Subtraction of convolved continuum - regridded ff estimate
    ffsub_image_B6 = subtr_ff(conv_cont_B6,regr_mmff_B6)

    # PB correction of ff estimate, error map, and subtracted map
    pb_im_B6 = cont_image_B6.replace('.fits','.pb.fits')
    pbcor_mmff_B6 = pbcor(regr_mmff_B6,pb_im_B6,squeeze=True)
    pbcor_err_B6 = pbcor(regr_err_B6,pb_im_B6,squeeze=True)
    pbcor_ffsub_B6 = pbcor(ffsub_image_B6,pb_im_B6)



for field in nondetections:
    print('\n**********')
    print('Free-free estimation for non-detection in '+field)
    print('**********')

    # Correct files for RL cube and continuum images. 
    input_cube = os.path.join(data_path,nondetections[field][0])

    # Moment 0 of RL image with custom velocity range
    mom0_image = mom0(input_cube,vel_range=nondetections[field][3])

    ##################################################
    # Fill with zeros:
    print('Filling free-free image with zeros:')
    ##################################################

    # Free-free estimation image from RL mom0, and error image:
    mmff_image = zero_image(mom0_image)
