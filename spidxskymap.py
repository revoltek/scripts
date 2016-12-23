#!/usr/bin/python

import os, sys, glob
import numpy as np
import astropy.io.fits as pyfits
from astropy import wcs
from lofar import bdsm
from scipy.ndimage.measurements import label
from scipy.ndimage.measurements import center_of_mass
from linearfit import twopoint_spidx_bootstrap
from astropy.table import Table

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

    print "Matching catalogues..."
    sources = Table.read('spidx_cat.txt', format='ascii')
    idx, _, _ = match_coordinates_sky(SkyCoord(sources['ra']*u.deg, sources['dec']*u.deg),\
                                      SkyCoord(centers['ra'], centers['dec']))

    print "Removing duplicates..."
    idx_duplicates = [] 
    for i, source in enumerate(sources):
        # check if closest
#        print idx[i], source['mask'], centers[int(idx[i])]['name']
        if source['mask'] != centers[int(idx[i])]['name']:
            idx_duplicates.append(i)
#            print "Removing source ", i
    print "Removing a total of", len(idx_duplicates), "sources."
    sources.remove_rows(idx_duplicates)
    sources.remove_columns('mask')
    sources.write('spidx_cat-nodup.txt', format='ascii.commented_header')

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
    #image_gaus_nvss = image_nvss.replace('.fits','-gaus.fits').replace('NVSS','NVSS/gaus',1)
    image_isl_nvss = image_nvss.replace('.fits','-isl.fits').replace('NVSS','NVSS/isl',1)
    cat_srl_nvss = image_tgss.replace('.fits','-srl.fits').replace('NVSS','NVSS/catalog',1)
    if not os.path.exists(image_rms_nvss) or not os.path.exists(image_isl_nvss):
        c = bdsm.process_image(image_nvss, frequency=1400e6, rms_box=(102,34), advanced_opts=True, group_tol=0.5, thresh_isl=3, thresh_pix=4)
        c.export_image(outfile=image_rms_nvss, img_type='rms', clobber=True)
        #c.export_image(outfile=image_gaus_nvss, img_type='gaus_model', clobber=True)
        c.export_image(outfile=image_isl_nvss, img_type='isl_mask', clobber=True)
        # write source catalog
        c.write_catalog(outfile=cat_srl_nvss, catalog_type='srl', format='fits', clobber=True)
        # clip the gaus images and create island masks
        # we use gaus images because they catch better the extended emission than isl_mask
        #clip_gaus(image_gaus_nvss, image_rms_nvss)

    image_rms_tgss = image_tgss.replace('.fits','-rms.fits').replace('TGSS','TGSS/rms',1)
    #image_gaus_tgss = image_tgss.replace('.fits','-gaus.fits').replace('TGSS','TGSS/gaus',1)
    image_isl_tgss = image_tgss.replace('.fits','-isl.fits').replace('TGSS','TGSS/isl',1)
    cat_srl_tgss = image_tgss.replace('.fits','-srl.fits').replace('TGSS','TGSS/catalog',1)
    if not os.path.exists(image_rms_tgss) or not os.path.exists(image_isl_tgss):
        c = bdsm.process_image(image_tgss, frequency=147e6, rms_box=(102,34), advanced_opts=True, group_tol=0.5, thresh_isl=3, thresh_pix=4)
        c.export_image(outfile=image_rms_tgss, img_type='rms', clobber=True)
        #c.export_image(outfile=image_gaus_tgss, img_type='gaus_model', clobber=True)
        c.export_image(outfile=image_isl_tgss, img_type='isl_mask', clobber=True)
        c.write_catalog(outfile=cat_srl_tgss, catalog_type='srl', format='fits', clobber=True)
        #clip_gaus(image_gaus_tgss, image_rms_tgss)

    image_mask = wdir+'masks/'+os.path.basename(image_nvss).replace('NVSS_','mask_')
    if not os.path.exists(image_mask):
        
        # create a mask with 1 where at least one survey has a detection, 0 otherwise
        print "Making mask..."
        with pyfits.open(image_isl_nvss) as fits_nvss:
            with pyfits.open(image_isl_tgss) as fits_tgss:
                fits_nvss[0].data = ((fits_nvss[0].data == 1) | (fits_tgss[0].data == 1)).astype(float)
                fits_nvss.writeto(image_mask, clobber=True)

        # Separate combined island mask into islands, each has a number
        with pyfits.open(image_mask) as fits_mask:
            blobs, number_of_blobs = label(fits_mask[0].data.astype(int))
        print "# of islands found:", number_of_blobs

#        data_nvss =  np.array( pyfits.getdata(image_nvss, 0) )
#        data_rms_nvss =  np.array( pyfits.getdata(image_rms_nvss, 0) )
#        data_tgss =  np.array( pyfits.getdata(image_tgss, 0) )
#        data_rms_tgss =  np.array( pyfits.getdata(image_rms_tgss, 0) )
#        area = 10.1978092553 # beam area in pixels - with beam: 45"x45" and pixels size: 15"x15"

        # opening catalogue as astropy Table
        if not os.path.exists('spidx_cat.fits'):
            t = Table(names=('ra','dec','S_nvss','S_nvss_err','rms_nvss','S_tgss','S_tgss_err','rms_tgss','spidx','status','mask'),\
                    dtype=('f4', 'i4', 'S2'))
        else:
            t = Table.open('spidx_cat.fits')

#            with open('spidx_cat.txt', 'a') as cat:
#                # deg deg Jy Jy Jy Jy - - pixels -
#                cat.write('#ra dec S_nvss S_nvss_err rms_nvss S_tgss S_tgss_err rms_tgss spidx spidx_err f2p_ratio_nvss f2p_ratio_tgss size status mask\n')
#                cat.write('#fluxes are in mJy\n')

        t_nvss = Table.open(cat_srl_nvss)
        t_nvss['s2n'] = t_nvss['Total_flux']/t_nvss['E_Total_flux']
        t_tgss = Table.open(cat_srl_tgss)
        t_tgss['s2n'] = t_tgss['Total_flux']/t_tgss['E_Total_flux']

        # Cross-match NVSS-TGSS source catalogue
        # match_coordinates_sky() gives an idx of the 2nd catalogue for each source of the 1st catalogue
        idx_match_tgss, sep, _ = match_coordinates_sky(SkyCoord(t_nvss['RA'], t_nvss['DEC']),\
                                          SkyCoord(t_tgss['RA'], t_tgss['DEC']))
        idx_match_tgss = idx_match_nvss[sep<15*u.arcsec]
        #idx_match_nvss, sep, _ = match_coordinates_sky(SkyCoord(t_tgss['RA'], t_tgss['DEC']),\
        #                                  SkyCoord(t_nvss['RA'], t_nvss['DEC']))
        #idx_match_nvss = idx_match_tgss[sep<15*u.arcsec]

        # Check if S/N > 5
        t_tgss_matched = t_tgss[idx_match_tgss[sep<=15*u.arcsec]]
        t_nvss_matched = t_nvss[sep<=15*u.arcsec]
        # TODO: write out matched sources

        # find blob for each unmatched source
        w = wcs.WCS(pyfits.open(image_nvss)[0].header)
        t_tgss_unmatched = t_tgss[idx_match_tgss[sep>15*u.arcsec]]
        t_nvss_unmatched = t_nvss[sep>15*u.arcsec]
        x_nvss, y_nvss, _, _ = w.all_world2pix([[t_tgss_unmatched['RA'],t_tgss_unmatched['DEC'],0,0]], 0)[0]
        x_tgss, y_tgss, _, _ = w.all_world2pix([[t_nvss_unmatched['RA'],t_nvss_unmatched['DEC'],0,0]], 0)[0]

        # For sources with no matches check in which island they are
        for s in no_match:
            # if alone -> add as upper/lower limit

            # if in an island with other unmatched sources
                # if same freq -> add as separate upper/lower limits

                # if different freq -> combine

            # if in an island with matched sources -> combine

            # if in an island with matched and umatched sources -> combine


        for s in xrange(1,number_of_blobs+1):
            idx = (blobs == s) # faster than np.where()
            flux_nvss = np.sum(data_nvss[idx])/area
            flux_tgss = np.sum(data_tgss[idx])/area
            flux_peak_nvss = np.max(data_nvss[idx])
            flux_peak_tgss = np.max(data_tgss[idx])
            rms_nvss = np.nanmean(data_rms_nvss[idx])
            rms_tgss = np.nanmean(data_rms_tgss[idx])
            flux_err_nvss = rms_nvss * np.sqrt( np.sum(idx)/area )
            flux_err_tgss = rms_tgss * np.sqrt( np.sum(idx)/area )
            snr_nvss = flux_nvss/flux_err_nvss
            if snr_nvss<0: snr_nvss = 0
            snr_tgss = flux_tgss/flux_err_tgss
            if snr_tgss<0: snr_tgss = 0

            # get rid of edges or partially covered images
            if np.isnan(snr_nvss) or np.isnan(snr_tgss): continue

            # good
            #elif flux_nvss>3*flux_err_nvss and flux_tgss>3*flux_err_tgss: 
            if snr_nvss+snr_tgss > 5 and (snr_nvss > 2 and snr_tgss > 2):
                status = 'g'
                spidx, spidx_err = twopoint_spidx_bootstrap([147.,1400.], [flux_tgss,flux_nvss], [flux_err_tgss, flux_err_nvss], niter=1000)
                f2p_ratio_nvss = flux_nvss/flux_peak_nvss
                f2p_ratio_tgss = flux_tgss/flux_peak_tgss

            # upper limit
            #if flux_nvss<3*flux_err_nvss and flux_tgss>3*flux_err_tgss:
            elif snr_tgss > 3:
                status = 'u'
                flux_nvss = 3*flux_err_nvss
                spidx = np.log10(flux_tgss/flux_nvss)/np.log10(147./1400.)
                spidx_err = -1
                f2p_ratio_nvss = -1
                f2p_ratio_tgss = flux_tgss/flux_peak_tgss

            # lower limit
            #elif flux_nvss>3*flux_err_nvss and flux_tgss<3*flux_err_tgss:
            elif snr_nvss > 3:
                status = 'l'
                flux_tgss = 3*flux_err_tgss
                spidx = np.log10(flux_tgss/flux_nvss)/np.log10(147./1400.)
                spidx_err = -1
                f2p_ratio_nvss = flux_nvss/flux_peak_nvss
                f2p_ratio_tgss = -1

            else:
                print "!",
                continue

            # get coords
            _, _, x, y = center_of_mass(idx)
            ra, dec, _, _ = w.all_pix2world([[y,x,0,0]], 0)[0]

            #print 'Flux NVSS:', flux_nvss, '+/-',flux_err_nvss,' - TGSS:', flux_tgss, '+/-', flux_err_nvss, '- status='+status
            t.write('%.4f %.4f %.3f %.3f %.3f %.3f %.3f %.3f %.4f %.4f %.2f %.2f %i %s %s\n' \
                % (ra, dec, flux_nvss*1e3, flux_err_nvss*1e3, rms_nvss*1e3, flux_tgss*1e3, flux_err_tgss*1e3, rms_tgss*1e3, \
                spidx, spidx_err, f2p_ratio_nvss, f2p_ratio_tgss, np.sum(idx), status, os.path.basename(image_mask)))

        print "Writing catalogue"
        t.write('spidx-cat.fits')

    sys.exit()

    #######################
    # spectral index map
    # TODO: rewrite headers
    image_spidx = wdir+'spidx/'+os.path.basename(image_nvss).replace('NVSS_','spidx_')
    image_spidx_u = wdir+'spidx_u/'+os.path.basename(image_nvss).replace('NVSS_','spidxU_')
    image_spidx_l = wdir+'spidx_l/'+os.path.basename(image_nvss).replace('NVSS_','spidxL_')
    image_spidx_err = wdir+'spidx_err/'+os.path.basename(image_nvss).replace('NVSS_','spidx_').replace('.fits','-err.fits')
    if not os.path.exists(image_spidx) or not os.path.exists(image_spidx_err) or not os.path.exists(image_spidx_u) or not os.path.exists(image_spidx_l):
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
            fits_nvss[0].data[idx_u] = data_spidx_u
            fits_nvss.writeto(image_spidx_u, clobber=True)

            fits_nvss[0].data[:] = np.nan
            fits_nvss[0].data[idx_l] = data_spidx_l
            fits_nvss.writeto(image_spidx_l, clobber=True)

            fits_nvss[0].data[:] = np.nan
            fits_nvss[0].data[idx_g] = data_spidx_g_err
            fits_nvss.writeto(image_spidx_err, clobber=True)

remove_duplicates()

print "All done!"
