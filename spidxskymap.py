#!/usr/bin/python

import os, sys, glob
import numpy as np
import astropy.io.fits as pyfits
from astropy import wcs
from lofar import bdsm
from scipy.ndimage.measurements import label
from scipy.ndimage.measurements import center_of_mass
from linearfit import twopoint_spidx_bootstrap

# options
wdir = '/home/fdg/data/spidxskymap/'

def clip_gaus(img, rmsimg, nsigma=1):
    """
    Find the mean rms from rmsimg and clip fitsimage
    to nsigma times this value. Then return the mask.
    """
    with pyfits.open(rmsimg) as fits:
        rms = np.mean(fits[0].data)
    print "Avg rms for ",rmsimg,":", rms

    fitsimg = pyfits.open(img)
    with pyfits.open(img) as fits:
        data = fits[0].data
        mask = (data > nsigma*rms).astype(int)
        fits[0].data = mask
        fits.writeto(img, clobber=True)
    
    return img

# assume NVSS/TGSS files have same name apart from initial part
images_nvss = glob.glob(wdir+'NVSS/*fits')

for image_nvss in images_nvss:
    image_tgss = image_nvss.replace('NVSS','TGSS')
    if not os.path.exists(image_tgss):
        print "TGSS file: ", image_tgss, "does not exist, continue."
        continue

    print "Work on:", image_nvss, image_tgss

    # source finder
    image_rms_nvss = image_nvss.replace('.fits','-rms.fits').replace('NVSS','NVSS/rms',1)
    image_gaus_nvss = image_nvss.replace('.fits','-gaus.fits').replace('NVSS','NVSS/gaus',1)
    if not os.path.exists(image_rms_nvss) or not os.path.exists(image_gaus_nvss):
        c = bdsm.process_image(image_nvss, frequency=1400e6, group_by_isl=True, mean_map='zero')
        c.export_image(outfile=image_rms_nvss, img_type='rms', clobber=True)
        c.export_image(outfile=image_gaus_nvss, img_type='gaus_model', clobber=True)
        # clip the gaus images and create island masks
        # we use gaus images because they catch better the extended emission than isl_mask
        clip_gaus(image_gaus_nvss, image_rms_nvss)

    image_rms_tgss = image_tgss.replace('.fits','-rms.fits').replace('TGSS','TGSS/rms',1)
    image_gaus_tgss = image_tgss.replace('.fits','-gaus.fits').replace('TGSS','TGSS/gaus',1)
    if not os.path.exists(image_rms_tgss) or not os.path.exists(image_gaus_tgss):
        c = bdsm.process_image(image_tgss, frequency=147e6, group_by_isl=True, mean_map='zero')
        c.export_image(outfile=image_rms_tgss, img_type='rms', clobber=True)
        c.export_image(outfile=image_gaus_tgss, img_type='gaus_model', clobber=True)
        clip_gaus(image_gaus_tgss, image_rms_tgss)

    # do an OR between the two island masks
    # create a mask with 1 where at least one survey has a detection, 0 otherwise
    print "Making mask..."
    image_mask = wdir+'masks/'+os.path.basename(image_nvss).replace('NVSS_','mask_')
    catalog = wdir+'catalog.txt'
    if not os.path.exists(image_mask):
        with pyfits.open(image_gaus_nvss) as fits_nvss:
            with pyfits.open(image_gaus_tgss) as fits_tgss:
                fits_nvss[0].data = ((fits_nvss[0].data == 1) | (fits_tgss[0].data == 1)).astype(float)
                fits_nvss.writeto(image_mask, clobber=True)

        print "Making catalogue..."
        with pyfits.open(image_mask) as fits_mask:
            blobs, number_of_blobs = label(fits_mask[0].data.astype(int))
        print "# of source found:", number_of_blobs
        data_nvss =  np.array( pyfits.getdata(image_nvss, 0) )
        data_rms_nvss =  np.array( pyfits.getdata(image_rms_nvss, 0) )
        data_tgss =  np.array( pyfits.getdata(image_tgss, 0) )
        data_rms_tgss =  np.array( pyfits.getdata(image_rms_tgss, 0) )
        area = 10.1978092553 # beam area in pixels - with beam: 45"x45" and pixels size: 15"x15"
        if not os.path.exists('spidx_cat.txt'):
            with open('spidx_cat.txt', 'a') as cat:
                cat.write('#ra dec S_nvss S_err-nvss S_tgss S_err-tgss spidx spidx_err status\n')
        with open('spidx_cat.txt', 'a') as cat:
            w = wcs.WCS(pyfits.open(image_nvss)[0].header)
            for s in xrange(1,number_of_blobs+1):
                if s % 10 == 0: print '.',
                status = 'g'
                idx = np.where(blobs == s)
                flux_nvss = data_nvss[idx].sum()/area
                flux_tgss = data_tgss[idx].sum()/area
                # len(idx[0]) gives the number of pixels of the source
                rms_nvss = data_rms_nvss[idx].mean()*np.sqrt(len(idx[0])/area)
                rms_tgss = data_rms_tgss[idx].mean()*np.sqrt(len(idx[0])/area)
                # upper limit
                if flux_nvss<3*rms_nvss:
                    status = 'u'
                    flux_nvss = 3*rms_nvss
                # lower limit
                if flux_tgss<3*rms_tgss:
                    status = 'l'
                    flux_tgss = 3*rms_tgss
                if flux_nvss<3*rms_nvss and flux_tgss<3*rms_tgss: 
                    print "ERROR: source undetected in both surveys?!"
                    sys.exit()

                # get coords
                _, _, x, y = center_of_mass((blobs == s).astype(int))
                ra, dec, _, _ = w.all_pix2world([[y,x,0,0]], 0)[0]

                spidx, spidx_err = twopoint_spidx_bootstrap([147.,1400.], [flux_tgss,flux_nvss], [rms_tgss, rms_nvss], niter=10000)
                cat.write('%.4f %.4f %.5f %.5f %.5f %.5f %.2f %.2f %s\n' % (ra, dec, flux_nvss, rms_nvss, flux_tgss, rms_tgss, spidx, spidx_err, status))
                #print 'Flux NVSS:', flux_nvss, '+/-',rms_nvss,' - TGSS:', flux_tgss, '+/-', rms_nvss,\
                #        ' - spidx: ', spidx, '+/-', spidx_err, '- status='+status
            print "done."

    #######################
    # spectral index map
    print "Makign spidx map..."
    image_spidx = wdir+'spidx/'+os.path.basename(image_nvss).replace('NVSS_','spidx_')
    image_spidx_ul = wdir+'spidx/'+os.path.basename(image_nvss).replace('NVSS_','spidxUL_')
    image_spidx_err = wdir+'spidx/'+os.path.basename(image_nvss).replace('NVSS_','spidx_').replace('.fits','-err.fits')
    if not os.path.exists(image_spidx) or not os.path.exists(image_spidx_err) or not os.path.exists(image_spidx_ul):
        idx = (pyfits.getdata(image_nvss, 0) > 3 * pyfits.getdata(image_rms_nvss, 0)) | (pyfits.getdata(image_tgss, 0) > 3 * pyfits.getdata(image_rms_tgss, 0))
        idx_g = (pyfits.getdata(image_nvss, 0) > 3 * pyfits.getdata(image_rms_nvss, 0)) & (pyfits.getdata(image_tgss, 0) > 3 * pyfits.getdata(image_rms_tgss, 0))
        idx_u = (pyfits.getdata(image_nvss, 0) <= 3 * pyfits.getdata(image_rms_nvss, 0)) & (pyfits.getdata(image_tgss, 0) > 3 * pyfits.getdata(image_rms_tgss, 0))
        idx_l = (pyfits.getdata(image_nvss, 0) > 3 * pyfits.getdata(image_rms_nvss, 0)) & (pyfits.getdata(image_tgss, 0) <= 3 * pyfits.getdata(image_rms_tgss, 0))

        data_nvss = np.array( pyfits.getdata(image_nvss, 0) )
        data_rms_nvss =  np.array( pyfits.getdata(image_rms_nvss, 0) )
        data_nvss[idx_u] = 3*data_rms_nvss[idx_u] # set val for upper limits

        data_tgss = np.array( pyfits.getdata(image_tgss, 0) )
        data_rms_tgss =  np.array( pyfits.getdata(image_rms_tgss, 0) )
        data_tgss[idx_l] = 3*data_rms_nvss[idx_l] # set val for lower limits

        data_spidx_g = np.zeros(shape=data_nvss[idx_g].shape)
        data_spidx_g_err = np.zeros(shape=data_nvss[idx_g].shape)

#        for d in xrange(len(data_nvss[idx_g])):
#            if d % 10 == 0: print '.',
        data_spidx_g, data_spidx_g_err = twopoint_spidx_bootstrap([147.,1400.], [data_tgss[idx_g],data_nvss[idx_g]], \
                    [data_rms_tgss[idx_g],data_rms_nvss[idx_g]], niter=10000)
        # upper limits
        data_spidx_u = (data_tgss[idx_u]-data_nvss[idx_u])/np.log10(147./1400.)
        # lower limits
        data_spidx_l = (data_tgss[idx_l]-data_nvss[idx_l])/np.log10(147./1400.)

        # write spidx, spidx+upper/lower limits and spidx error map
        with pyfits.open(image_nvss) as fits_nvss:
            fits_nvss[0].data[:] = np.nan
            fits_nvss[0].data[idx_g] = data_spidx_g
            fits_nvss.writeto(image_spidx, clobber=True)

            fits_nvss[0].data[:] = np.nan
            fits_nvss[0].data[idx_g] = data_spidx_g
            fits_nvss[0].data[idx_u] = data_spidx_u
            fits_nvss[0].data[idx_l] = data_spidx_l
            fits_nvss.writeto(image_spidx_ul, clobber=True)

            fits_nvss[0].data[:] = np.nan
            fits_nvss[0].data[idx_g] = data_spidx_g_err
            fits_nvss[0].data[idx_u] = -1
            fits_nvss[0].data[idx_l] = -2
            fits_nvss.writeto(image_spidx_err, clobber=True)
        print "done."

#    c = bdsm.process_image(image_nvss, frequency=1400e6, group_by_isl=True, detection_image=image_mask, beam=(1.25000e-02, 1.25000e-02, 0.0), mean_map='zero', rms_map=False, rms_value=0.1)
#    cat = image_nvss.replace('.fits','.cat')
#    if c.sources == []:
#        open(cat, 'a').close()
#    else:
#        c.write_catalog(outfile=cat, format='ascii', catalog_type='srl', clobber=True)
#
#    c = bdsm.process_image(image_tgss, frequency=147e6, group_by_isl=True, detection_image=image_mask, beam=(1.25000e-02, 1.25000e-02, 0.0), mean_map='zero', rms_map=False, rms_value=0.1)
#    cat = image_tgss.replace('.fits','.cat')
#    if c.sources == []:
#        open(cat, 'a').close()
#    else:
#        c.write_catalog(outfile=cat, format='ascii', catalog_type='srl', clobber=True)
