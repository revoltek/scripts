#!/usr/bin/python

import os, sys, glob
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
import lofar.bdsm as bdsm
from scipy.ndimage.measurements import label
from linearfit import linear_fit_bootstrap

wdir = '/home/fdg/data/spidxskymap/'

def clip_gaus(img, rmsimg, nsigma=1):
    """
    Find the mean rms from rmsimg and clip fitsimage
    to nsigma times this value. Then return the mask.
    """
    # TODO: use with?
    rms = np.mean(pyfits.open(rmsimg)[0].data)
    print "Avg rms for ",rmsimg,":", rms

    fitsimg = pyfits.open(img)
    data = fitsimg[0].data
    mask = (data > nsigma*rms).astype(int)
    fitsimg[0].data = mask
    fitsimg.writeto(img, clobber=True)
    fitsimg.close()
    return img

# assume file have same name apart from initial
images_nvss = glob.glob(wdir+'NVSS/*fits')

for image_nvss in images_nvss:
    image_tgss = image_nvss.replace('NVSS','TGSS')
    if not os.path.exists(image_tgss):
        print "TGSS file: ", image_tgss, "does not exist, continue."
        continue
    print "Work on", image_nvss, image_tgss

    # source finder
    image_rms_nvss = image_nvss.replace('.fits','-rms.fits').replace('NVSS','NVSS/rms',1)
    image_gaus_nvss = image_nvss.replace('.fits','-gaus.fits').replace('NVSS','NVSS/gaus',1)
    if not os.path.exists(image_rms_nvss) or not os.path.exists(image_gaus_nvss):
        c = bdsm.process_image(image_nvss, frequency=1400e6, group_by_isl=True, mean_map='zero')
        c.export_image(outfile=image_rms_nvss, img_type='rms', clobber=True)
        c.export_image(outfile=image_gaus_nvss, img_type='gaus_model', clobber=True)

    image_rms_tgss = image_tgss.replace('.fits','-rms.fits').replace('TGSS','TGSS/rms',1)
    image_gaus_tgss = image_tgss.replace('.fits','-gaus.fits').replace('TGSS','TGSS/gaus',1)
    if not os.path.exists(image_rms_tgss) or not os.path.exists(image_gaus_tgss):
        c = bdsm.process_image(image_tgss, frequency=147e6, group_by_isl=True, mean_map='zero')
        c.export_image(outfile=image_rms_tgss, img_type='rms', clobber=True)
        c.export_image(outfile=image_gaus_tgss, img_type='gaus_model', clobber=True)

    # clip the gaus images and create island masks
    # we use gaus images because they catch better the extended emission than isl_mask
    clip_gaus(image_gaus_nvss, image_rms_nvss)
    clip_gaus(image_gaus_tgss, image_rms_tgss)

    # do an OR between the two island masks
    print "Making mask..."
    image_mask = wdir+'masks/'+os.path.basename(image_nvss).replace('NVSS_','mask_')
    catalog = wdir+'catalog.txt'
    if not os.path.exists(image_mask):
        fits_nvss = pyfits.open(image_gaus_nvss)
        fits_tgss = pyfits.open(image_gaus_tgss)
        data = ((fits_nvss[0].data == 1) | (fits_tgss[0].data == 1)).astype(float)
#        fits_nvss[0].data[:] = np.nan
#        fits_nvss[0].data[(fits_nvss[0].data == 1) | (fits_tgss[0].data == 1)] = 1
        fits_nvss.writeto(image_mask, clobber=True)
        fits_tgss.close()
        fits_nvss.close()

    print "Making catalogue..."
    fits_mask = pyfits.open(image_mask)
    blobs, number_of_blobs = label(fits_mask[0].data.astype(int))
    print "# of source found:", number_of_blobs
    data_nvss =  np.array( pyfits.getdata(image_nvss, 0) )
    data_rms_nvss =  np.array( pyfits.getdata(image_rms_nvss, 0) )
    data_tgss =  np.array( pyfits.getdata(image_tgss, 0) )
    data_rms_tgss =  np.array( pyfits.getdata(image_rms_tgss, 0) )
    area = 10.1978092553 # beam area in pixels - with beam: 45"x45" and pixels size: 15"x15"
    with open('spidx_cat.txt', 'a') as cat:
        for s in xrange(1,number_of_blobs+1):
            status = 'g'
            idx = np.where(blobs == s)
            flux_nvss = data_nvss[idx].sum()/area
            flux_tgss = data_tgss[idx].sum()/area
            f = np.sqrt(len(idx[0])/area) # len(idx[0]) gives the number of pixels of the source
            rms_nvss = data_rms_nvss[idx].mean()*f
            rms_tgss = data_rms_tgss[idx].mean()*f
            if flux_nvss<3*rms_nvss:
                status = 'u'
                flux_nvss = 3*rms_nvss
            if flux_tgss<3*rms_tgss:
                status = 'l'
                flux_tgss = 3*rms_tgss
            if flux_nvss<3*rms_nvss and flux_tgss<3*rms_tgss: 
                print "ERROR"
                sys.exit()
            spidx, _, spidx_err, _ = linear_fit_bootstrap(np.log10([147.,1400.]), \
                    np.log10([flux_tgss,flux_nvss]), [0.434*rms_tgss/flux_tgss, 0.434*rms_nvss/flux_nvss], niter=100)
            cat.write('%f %f %f %f %f %f %s\n' % (flux_nvss, rms_nvss, flux_tgss, rms_tgss, spidx, spidx_err, status))
            print 'Flux NVSS:', flux_nvss, '+/-',rms_nvss,' - TGSS:', flux_tgss, '+/-', rms_nvss,' - spidx: ', spidx, '+/-', spidx_err, '- status='+status

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

        # to be quick update only those pixels that will eventually be used 
        data_rms_nvss[idx] = 0.434*data_rms_nvss[idx]/data_nvss[idx]
        data_rms_tgss[idx] = 0.434*data_rms_tgss[idx]/data_tgss[idx]
        data_nvss[idx] = np.log10( data_nvss[idx] )
        data_tgss[idx] = np.log10( data_tgss[idx] )

        data_spidx_g = np.zeros(shape=data_nvss[idx_g].shape)
        data_spidx_g_err = np.zeros(shape=data_nvss[idx_g].shape)

        for d in xrange(len(data_nvss[idx_g])):
            if d % 10 == 0: print '.',
            data_spidx_g[d], _, data_spidx_g_err[d], _ = linear_fit_bootstrap(np.log10([147.,1400.]), \
                    [data_tgss[idx_g][d],data_nvss[idx_g][d]], [data_rms_tgss[idx_g][d],data_rms_nvss[idx_g][d]], niter=100)
        # upper limits
        data_spidx_u = (data_tgss[idx_u]-data_nvss[idx_u])/np.log10(147./1400.)
        # lower limits
        data_spidx_l = (data_tgss[idx_l]-data_nvss[idx_l])/np.log10(147./1400.)

        # write spidx and error map
        fits_nvss = pyfits.open(image_nvss)
        fits_nvss[0].data[:] = np.nan
        fits_nvss[0].data[idx_g] = data_spidx_g
        fits_nvss.writeto(image_spidx, clobber=True)
        fits_nvss = pyfits.open(image_nvss)
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
        fits_nvss.close()

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
