'''
Script that aglomerates tools for free-free estimates from RL cubes, as
well as subtraction of these estimates from ALMA-IMF continuum bands.
'''

import numpy as np

from astropy.io import fits
from astropy import units as u
from astropy.convolution import convolve, convolve_fft
from astropy import wcs
from astropy.stats import median_absolute_deviation, mad_std
from astropy.stats import SigmaClip

from spectral_cube import SpectralCube
from spectral_cube import VaryingResolutionSpectralCube
from radio_beam import Beam
from reproject import reproject_interp




# CREATE MOMENT 0 OF CUBE
def mom0(input_im,vel_range,beam_threshold=0.05,**kwargs):
    #cube = VaryingResolutionSpectralCube.read(input_im)
    cube = SpectralCube.read(input_im)
    cube.beam_threshold = beam_threshold
    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')

    #noise_cube = cube.spectral_slab((vel_range[0]-10)*u.km/u.s ,vel_range[0]*u.km/u.s)
    #noise_data = noise_cube.unmasked_data[:,:,:]
    #rms_cube = mad_std(noise_data, ignore_nan=True) #MAD

    cube = cube.spectral_slab(vel_range[0]*u.km/u.s ,vel_range[1]*u.km/u.s)
    moment_0 = cube.moment(order=0)
    #print('Channels used in mom0 = {0}, noise per channel = {1}'.format(cube.shape[0], rms_cube))
    #rms_mom0 = np.sqrt(cube.shape[0])*rms_cube #Jy/beam km/s

    output_mom0 = input_im.replace('.fits','.m0.fits')
    moment_0.write(output_mom0, overwrite=True)
    print('Created '+output_mom0)
    return output_mom0

# CREATE MOMENT 1 OF CUBE
def mom1and2(input_im,vel_range,thr=5.0,beam_threshold=0.05,**kwargs):
    #cube = VaryingResolutionSpectralCube.read(input_im)
    cube = SpectralCube.read(input_im)
    cube.beam_threshold = beam_threshold
    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')

    noise_cube = cube.spectral_slab((vel_range[0]-10)*u.km/u.s ,vel_range[0]*u.km/u.s)
    noise_data = noise_cube.unmasked_data[:,:,:]
    rms_cube = mad_std(noise_data, ignore_nan=True) #MAD

    cube = cube.spectral_slab(vel_range[0]*u.km/u.s ,vel_range[1]*u.km/u.s)
    masked_cube=cube.with_mask(cube > thr*rms_cube)

    moment_1 = masked_cube.moment(order=1)
    #print('Channels used in mom1 = {0}, noise per channel = {1}'.format(cube.shape[0], rms_cube))
    output_mom1 = input_im.replace('.fits','.m1.fits')
    moment_1.write(output_mom1, overwrite=True)
    print('Created '+output_mom1)

    moment_2 = masked_cube.moment(order=2)
    #print('Channels used in mom1 = {0}, noise per channel = {1}'.format(cube.shape[0], rms_cube))
    output_mom2 = input_im.replace('.fits','.m2.fits')
    moment_2.write(output_mom2, overwrite=True)
    print('Created '+output_mom2)

    return output_mom1, output_mom2



# CONVOLUTION OF HIGH-RES IMAGE TO LOW-RES IMAGE
def convol(input_im,template_im,**kwargs):
    # Read in images and headers
    fh1 = fits.open(input_im)
    header1 = fits.getheader(input_im)
    beam1 = Beam.from_fits_header(header1)
    print(beam1)
    fh2 = fits.open(template_im)
    header2 = fits.getheader(template_im)
    beam2 = Beam.from_fits_header(header2)
    print(beam2)
    # Find and create convolution kernel
    try:
        conv_beam = beam2.deconvolve(beam1)
        print(conv_beam)
        pix_scale = header1['CDELT2']*u.deg
        pix_scale = pix_scale.to(u.arcsec)
        conv_beam_kernel = conv_beam.as_kernel(pix_scale)
        # Convolve in Fourier space. Need to scale to Jy/pix before
        abeam1 = 1.442*(np.pi/4)*header1['BMAJ']*header1['BMIN']*3600**2
        abeam2 = 1.442*(np.pi/4)*header2['BMAJ']*header2['BMIN']*3600**2
        apix  = (header1['CDELT2']*3600)**2
        print('apix/abeam1 = {:.3f}'.format(apix/abeam1))
        print('abeam2/apix = {:.3f}'.format(abeam2/apix))
        convldata = np.copy(fh1[0].data)
        convldata[0,0,:,:]=(abeam2/apix)*convolve_fft((apix/abeam1)
        *fh1[0].data[0,0,:,:],conv_beam_kernel,preserve_nan=True)
        #convldata = (abeam2/apix)*convolve_fft((apix/abeam1)
        #*fh1[0].data, conv_beam_kernel, preserve_nan=True)
        #convldata = np.where(np.abs(convldata) < 1e-5, np.nan, convldata)
    except:
        #raise Exception("Can't convolve")
        convldata = np.copy(fh1[0].data)
    # Update header and save to fits
    #del header1['HISTORY']
    header1['BMAJ'] = header2['BMAJ']
    header1['BMIN'] = header2['BMIN']
    header1['BPA'] = header2['BPA']
    conv_image = input_im.replace('.fits','.conv.fits')
    fits.writeto(conv_image,convldata,header1,overwrite=True)
    fh1.close()
    fh2.close()
    print('Created '+conv_image)
    return conv_image


# REGRIDDING
def regrid_tom0(input_im,template_im,**kwargs):
    # No change of frame, just pixel regridding
    fh1 = fits.open(input_im)
    fh2 = fits.open(template_im)
    w1 = wcs.WCS(fh1[0].header)
    w2 = wcs.WCS(fh2[0].header)
    repr,_ = reproject_interp((fh1[0].data.squeeze(), w1.celestial), w2.celestial, shape_out=fh2[0].data.squeeze().shape)
    outhead = w2.celestial.to_header()
    outhead['BMAJ'] = fh1[0].header['BMAJ']
    outhead['BMIN'] = fh1[0].header['BMIN']
    outhead['BPA'] = fh1[0].header['BPA']
    outhead['BTYPE'] = fh1[0].header['BTYPE']
    outhead['OBJECT'] = fh1[0].header['OBJECT']
    outhead['BUNIT'] = fh1[0].header['BUNIT']
    hdu = fits.PrimaryHDU(data=repr, header=outhead)
    # or
    #hdu = fits.PrimaryHDU(data=repr, header=fh2[0].header)
    regr_image = input_im.replace('.fits','.regr.fits')
    hdu.writeto(regr_image,overwrite=True)
    fh1.close()
    fh2.close()
    print('Created '+regr_image)
    return regr_image

def regrid_tocont(input_im,template_im,**kwargs):
    fh1 = fits.open(input_im)
    fh2 = fits.open(template_im)
    w1 = wcs.WCS(fh1[0].header)
    w2 = wcs.WCS(fh2[0].header)
    repr,_ = reproject_interp((fh1[0].data, w1.celestial), w2.celestial, shape_out=fh2[0].data.squeeze().shape)
    outhead = w2.celestial.to_header()
    outhead['BMAJ'] = fh1[0].header['BMAJ']
    outhead['BMIN'] = fh1[0].header['BMIN']
    outhead['BPA'] = fh1[0].header['BPA']
    outhead['BTYPE'] = fh1[0].header['BTYPE']
    outhead['OBJECT'] = fh1[0].header['OBJECT']
    outhead['BUNIT'] = fh1[0].header['BUNIT']
    outhead['SIGCLIP'] = fh1[0].header['SIGCLIP']
    hdu = fits.PrimaryHDU(data=repr, header=outhead)
    regr_image = input_im.replace('.fits','.regr.fits')
    hdu.writeto(regr_image,overwrite=True)
    fh1.close()
    fh2.close()
    print('Created '+regr_image)
    return regr_image

# Free-free estimation from RL
# input_im = RL mom0 image
# mask_thr = sigma clipping levelm default 5.
# Te = electron temperature in K, default 8000.
# NHe_NH = N(He+)/N(H+), default 0.08
def ff_rl(input_im,mask_thr=5.0,Te=8000.0,NHe_NH=0.08,freq_scale=False,new_freq=230.,
    fill_zeros=False,error_map=True,fsigma_Te=0.2,**kwargs):

    fh1 = fits.open(input_im)
    header1 = fh1[0].header
    sigma = mad_std(fh1[0].data, ignore_nan=True) #MAD
    print('MAD sigma in mom0 of RL cube = {} Jy beam^-1 km s^-1'.format(sigma))
    data1 = np.nan_to_num(fh1[0].data) #convert nans to zeroes
    fh1.close()

    #mask pixels above n*sigma
    if fill_zeros == True:
        data2 = np.where(data1 > mask_thr*sigma, data1, 0.0) #clipped data
    else:
        data2 = np.where(data1 > mask_thr*sigma, data1, np.nan) #clipped data

    #del header1['HISTORY']
    nu0 = header1['RESTFRQ']/1e9 #GHz
    header1['BUNIT'] = 'Jy/beam'
    header1['SIGCLIP'] = mask_thr

    #free-free estimation
    data3 = data2*(1/6.985e3)*(nu0**-1.1)*(Te**1.15)*(1+NHe_NH)
    str_freq = "{:.1f}".format(nu0)

    if freq_scale == True:
        header1['RESTFRQ'] = new_freq*1e9 #Hz for header
        data3 = data3*(new_freq/nu0)**(-0.1)
        str_freq = "{:.1f}".format(new_freq)

    mmff_image = input_im.replace('.fits','.ff'+str_freq+'GHz_'+str(mask_thr)+'sigma.fits')
    err_mmff_im = None
    hdu = fits.PrimaryHDU(data=data3, header=header1)
    hdu.writeto(mmff_image, overwrite=True)
    print('Created '+mmff_image)

    if error_map == True:
        m0_data = np.where(data1 > mask_thr*sigma, data1, np.nan)
        error_factor = np.sqrt((sigma/m0_data)**2 + fsigma_Te**2)
        error_data = data3*error_factor
        err_mmff_im = input_im.replace('.fits','.ff'+str_freq+'GHz_'+str(mask_thr)+'sigma.err.fits')
        hdu_err = fits.PrimaryHDU(data=error_data, header=header1)
        hdu_err.writeto(err_mmff_im, overwrite=True)
        print('Created '+err_mmff_im)

    return mmff_image, err_mmff_im

def subtr_ff(cont_im,mmff_im,**kwargs):
    cont_hdu = fits.open(cont_im)
    mmff_hdu = fits.open(mmff_im)
    cont_header = cont_hdu[0].header
    cont_data = cont_hdu[0].data
    sigma_clip = mmff_hdu[0].header['SIGCLIP']
    cont_header['SIGCLIP'] = sigma_clip
    mmff_data = mmff_hdu[0].data
    cont_hdu.close()
    mmff_hdu.close()
    subtr = cont_data - mmff_data
    hdu_out = fits.PrimaryHDU(data=subtr, header=cont_header)
    ffsub_image = cont_im.replace('.fits','.ffsub_'+str(sigma_clip)+'sigma.fits')
    hdu_out.writeto(ffsub_image,overwrite=True)
    print('Created '+ffsub_image)
    return ffsub_image

def pbcor(input_im, pb_im, squeeze=False, **kwargs):
    input_hdu = fits.open(input_im)
    pb_hdu = fits.open(pb_im)
    input_data = input_hdu[0].data
    input_hd = input_hdu[0].header
    if squeeze == True:
        pb_data = pb_hdu[0].data.squeeze()
    else:
        pb_data = pb_hdu[0].data
    input_hdu.close()
    pb_hdu.close()
    pbcor_data = input_data/pb_data
    pbcor_im = input_im.replace('.fits','.pbcor.fits')
    hdu_out = fits.PrimaryHDU(data=pbcor_data, header=input_hd)
    hdu_out.writeto(pbcor_im, overwrite=True)
    print('Created '+pbcor_im)
    return pbcor_im

def zero_image(input_im, squeeze=False, **kwargs):
    input_hdu = fits.open(input_im)
    input_data = input_hdu[0].data
    input_hd = input_hdu[0].header
    if squeeze == True:
        data = np.zeros(input_data.squeeze().shape)
    else:
        data = np.zeros(input_data.shape)
    input_hdu.close()
    output_im = input_im.replace('.fits','.zero.fits')
    hdu_out = fits.PrimaryHDU(data=data, header=input_hd)
    hdu_out.writeto(output_im, overwrite=True)
    print('Created '+output_im)
    return output_im

