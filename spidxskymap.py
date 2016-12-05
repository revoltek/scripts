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

def remove_duplicates(file_cat='spidx_cat.txt'):
    """
    Remove duplicates from overlapping regions.
    For each source check if the closest file-center is the one from which is extracted.
    If not it means the same source is closer in another file, delete it.
    """
    from astropy.table import Table
    from astropy.coordinates import match_coordinates_sky
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    # get all file centers
    print "Collecting centers..."
    centers = Table([[],[],[]], names=('name','ra','dec'), dtype=['S100',float,float])
    centers['ra'].unit = 'deg'
    centers['dec'].unit = 'deg'
    for i, mask_file in enumerate(glob.glob('masks/mask*.fits')):
        with pyfits.open(mask_file) as fits:
            head = fits[0].header
            ra = head['CRVAL1']
            dec = head['CRVAL2']
            centers.add_row([os.path.basename(mask_file), ra, dec])

    print "Removing duplicates..."
    sources = Table.read('spidx_cat.txt', format='ascii')
    #print centers,sources
    for i, source in enumerate(sources):
        if i % 10 == 0: print i
        idx, _, _ = match_coordinates_sky(SkyCoord(source['ra']*u.deg, source['dec']*u.deg),\
                                    SkyCoord(centers['ra'], centers['dec']))
        # check if closest
        if source['mask'] != centers[int(idx)]['name']:
            sources.remove_rows(idx)
            print "Removing source ", i
    sources.write('spidx_cat-nodup.txt', format='ascii')

def clip_gaus(img, rmsimg, nsigma=1):
    """
    Find the mean rms from rmsimg and clip fitsimage
    to nsigma times this value. Then return the mask.
    """
    with pyfits.open(rmsimg) as fits:
        rms = np.nanmean(fits[0].data)
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
        c = bdsm.process_image(image_nvss, frequency=1400e6, group_by_isl=True, rms_box=(102,34))
        c.export_image(outfile=image_rms_nvss, img_type='rms', clobber=True)
        c.export_image(outfile=image_gaus_nvss, img_type='gaus_model', clobber=True)
        # clip the gaus images and create island masks
        # we use gaus images because they catch better the extended emission than isl_mask
        clip_gaus(image_gaus_nvss, image_rms_nvss)

    image_rms_tgss = image_tgss.replace('.fits','-rms.fits').replace('TGSS','TGSS/rms',1)
    image_gaus_tgss = image_tgss.replace('.fits','-gaus.fits').replace('TGSS','TGSS/gaus',1)
    if not os.path.exists(image_rms_tgss) or not os.path.exists(image_gaus_tgss):
        c = bdsm.process_image(image_tgss, frequency=147e6, group_by_isl=True, rms_box=(102,34))
        c.export_image(outfile=image_rms_tgss, img_type='rms', clobber=True)
        c.export_image(outfile=image_gaus_tgss, img_type='gaus_model', clobber=True)
        clip_gaus(image_gaus_tgss, image_rms_tgss)

    # do an OR between the two island masks
    # create a mask with 1 where at least one survey has a detection, 0 otherwise
    image_mask = wdir+'masks/'+os.path.basename(image_nvss).replace('NVSS_','mask_')
    if not os.path.exists(image_mask):
        print "Making mask..."
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
                # deg deg Jy Jy Jy Jy - - pixels -
                cat.write('#ra dec S_nvss S_nvss_err S_tgss S_tgss_err spidx spidx_err size status mask\n')

        w = wcs.WCS(pyfits.open(image_nvss)[0].header)
        with open('spidx_cat.txt', 'a') as cat:
            for s in xrange(1,number_of_blobs+1):
                #if s % 10 == 0:
                #    print '.',
                idx = (blobs == s) # faster than np.where()
                flux_nvss = np.sum(data_nvss[idx])/area
                flux_tgss = np.sum(data_tgss[idx])/area
                rms_nvss = np.nanmean(data_rms_nvss[idx]) * np.sqrt( np.sum(idx)/area )
                rms_tgss = np.nanmean(data_rms_tgss[idx]) * np.sqrt( np.sum(idx)/area )
                snr_nvss = flux_nvss/rms_nvss
                if snr_nvss<0: snr_nvss = 0
                snr_tgss = flux_tgss/rms_tgss
                if snr_tgss<0: snr_tgss = 0

                # get rid of edges or partially covered images
                if np.isnan(snr_nvss) or np.isnan(snr_tgss): continue

                # good
                #elif flux_nvss>3*rms_nvss and flux_tgss>3*rms_tgss: 
                if snr_nvss+snr_tgss > 5 and (snr_nvss > 2 and snr_tgss > 2):
                    status = 'g'
                    spidx, spidx_err = twopoint_spidx_bootstrap([147.,1400.], [flux_tgss,flux_nvss], [rms_tgss, rms_nvss], niter=1000)

                # upper limit
                #if flux_nvss<3*rms_nvss and flux_tgss>3*rms_tgss:
                elif snr_tgss > 3:
                    status = 'u'
                    flux_nvss = 3*rms_nvss
                    spidx = np.log10(flux_tgss/flux_nvss)/np.log10(147./1400.)
                    spidx_err = -1

                # lower limit
                #elif flux_nvss>3*rms_nvss and flux_tgss<3*rms_tgss:
                elif snr_nvss > 3:
                    status = 'l'
                    flux_tgss = 3*rms_tgss
                    spidx = np.log10(flux_tgss/flux_nvss)/np.log10(147./1400.)
                    spidx_err = -1

                else:
                    print "!",
                    continue
    
                # get coords
                _, _, x, y = center_of_mass(idx)
                ra, dec, _, _ = w.all_pix2world([[y,x,0,0]], 0)[0]

                #print 'Flux NVSS:', flux_nvss, '+/-',rms_nvss,' - TGSS:', flux_tgss, '+/-', rms_nvss, '- status='+status
                cat.write('%.4f %.4f %.5f %.5f %.5f %.5f %.3f %.3f %i %s %s\n' \
                    % (ra, dec, flux_nvss, rms_nvss, flux_tgss, rms_tgss, spidx, spidx_err, np.sum(idx), status, os.path.basename(image_mask)))

        print "done."

    #######################
    # spectral index map
    # TODO: rewrite headers
    image_spidx = wdir+'spidx/'+os.path.basename(image_nvss).replace('NVSS_','spidx_')
    image_spidx_ul = wdir+'spidx_ul/'+os.path.basename(image_nvss).replace('NVSS_','spidxUL_')
    image_spidx_err = wdir+'spidx_err/'+os.path.basename(image_nvss).replace('NVSS_','spidx_').replace('.fits','-err.fits')
    if not os.path.exists(image_spidx) or not os.path.exists(image_spidx_err) or not os.path.exists(image_spidx_ul):
        print "Makign spidx map..."
        snr_nvss = pyfits.getdata(image_nvss, 0)/pyfits.getdata(image_rms_nvss, 0)
        snr_nvss[snr_nvss<0] = 0
        snr_tgss = pyfits.getdata(image_tgss, 0)/pyfits.getdata(image_rms_tgss, 0)
        snr_tgss[snr_tgss<0] = 0
        mask = pyfits.getdata(image_mask, 0)
        idx_g = (snr_nvss > 3) & (snr_tgss > 3) & (mask == 1)
        idx_u = (snr_tgss > 3) & (snr_nvss < 3) & (mask == 1)
        idx_l = (snr_nvss > 3) & (snr_tgss < 3) & (mask == 1)
        print "Number of good pixels: ", idx_g.sum()
        print "Number of upper lim: ", idx_u.sum()
        print "Number of lower lim: ", idx_l.sum()

        data_nvss = np.array( pyfits.getdata(image_nvss, 0) )
        data_rms_nvss =  np.array( pyfits.getdata(image_rms_nvss, 0) )
        data_nvss[idx_u] = 3*data_rms_nvss[idx_u] # set val for upper limits

        data_tgss = np.array( pyfits.getdata(image_tgss, 0) )
        data_rms_tgss =  np.array( pyfits.getdata(image_rms_tgss, 0) )
        data_tgss[idx_l] = 3*data_rms_nvss[idx_l] # set val for lower limits

        data_spidx_g = np.zeros(shape=data_nvss[idx_g].shape)
        data_spidx_g_err = np.zeros(shape=data_nvss[idx_g].shape)

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

remove_duplicates()

print "All done!"
