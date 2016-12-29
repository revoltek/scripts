#!/usr/bin/python

import os, sys, glob
import numpy as np
from scipy.ndimage.measurements import label
from scipy.ndimage.measurements import center_of_mass

from lofar import bdsm

import astropy.io.fits as pyfits
from astropy import wcs
from astropy.table import Table
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
import astropy.units as u

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
        c.export_image(outfile=image_isl_nvss, img_type='island_mask', clobber=True)
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
        c.export_image(outfile=image_isl_tgss, img_type='island_mask', clobber=True)
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

        # opening catalogue as astropy Table
        if not os.path.exists('spidx_cat.fits'):
            t = Table(names=('Source_id','RA','E_RA','DEC','E_DEC','Total_flux_NVSS','E_Total_flux_NVSS','Peak_flux_NVSS','E_Peak_flux_NVSS','Total_flux_TGSS','E_Total_flux_TGSS','Peak_flux_TGSS','E_Peak_flux_TGSS','spidx','E_spidx','Rms_NVSS','Rms_TGSS','S_Code','s2n','size','Mask'\
                      dtype=('i','f','f','f','f','f','f','f','f','f','S1','S100'))
        else:
            t = Table.read('spidx_cat.fits')

        t_nvss = Table.read(cat_srl_nvss)
        t_nvss['s2n'] = t_nvss['Total_flux']/t_nvss['E_Total_flux']
        t_nvss['s2n'].unit = ''
        t_tgss = Table.read(cat_srl_tgss)
        t_tgss['s2n'] = t_tgss['Total_flux']/t_tgss['E_Total_flux']
        t_tgss['s2n'].unit = ''

        # Cross-match NVSS-TGSS source catalogue
        # match_coordinates_sky() gives an idx of the 2nd catalogue for each source of the 1st catalogue
        idx_match, sep, _ = match_coordinates_sky(SkyCoord(t_nvss['RA'], t_nvss['DEC']),\
                                          SkyCoord(t_tgss['RA'], t_tgss['DEC']))

        match_sep = 15 # arcsec
        idx_matched_nvss = np.arange(0,len(t_nvss))[sep<=match_sep*u.arcsec] # nvss idx with a match
        idx_matched_tgss = idx_match[sep<=match_sep*u.arcsec] # tgss idx with a match
        idx_unmatched_nvss = np.arange(0,len(t_nvss))[sep>match_sep*u.arcsec] # nvss idx with NOT a match
        idx_unmatched_tgss = [ x for x in np.arange(0,len(t_tgss)) if x not in idx_matched_tgss ] # tgss idx with NOT a match
        assert len(set(idx_matched_tgss)) == len(idx_matched_tgss) # check no sources are cross-matched twice
        print "Matched sources - NVSS: %i/%i - TGSS %i/%i" % ( len(idx_matched_nvss), len(t_nvss), len(idx_matched_tgss), len(t_tgss) )

        # find blob number for each source
        w = wcs.WCS(pyfits.open(image_nvss)[0].header)
        x_nvss, y_nvss, _, _ = w.all_world2pix(t_nvss['RA'],t_nvss['DEC'],np.zeros(len(t_nvss)),np.zeros(len(t_nvss)), 0, ra_dec_order=True)
        x_tgss, y_tgss, _, _ = w.all_world2pix(t_tgss['RA'],t_tgss['DEC'],np.zeros(len(t_tgss)),np.zeros(len(t_tgss)), 0, ra_dec_order=True)
        val_blob_nvss = blobs[np.zeros(len(t_nvss)).astype(int),np.zeros(len(t_nvss)).astype(int),y_nvss.round().astype(int),x_nvss.round().astype(int)]
        val_blob_tgss = blobs[np.zeros(len(t_tgss)).astype(int),np.zeros(len(t_tgss)).astype(int),y_tgss.round().astype(int),x_tgss.round().astype(int)]
        assert (val_blob_nvss != 0).all() and (val_blob_tgss != 0).all() # all sources are into a blob

        # For sources with no matches check in which island they are
        for i, s in enumerate(idx_unmatched_tgss):
            blob = val_blob_tgss[s]
            idx_blob_nvss = np.where(val_blob_nvss == blob) # idx of NVSS sources in same blob of this unmatched source
            idx_blob_tgss = np.where(val_blob_tgss == blob) # idx of TGSS sources in same blob of this unmatched source

            idx_blob_matched_nvss = np.intersect1d(idx_blob_nvss, idx_matched_nvss)
            idx_blob_matched_tgss = np.intersect1d(idx_blob_tgss, idx_matched_tgss)
            idx_blob_unmatched_nvss = np.intersect1d(idx_blob_nvss, idx_unmatched_nvss)
            idx_blob_unmatched_tgss = np.intersect1d(idx_blob_tgss, idx_unmatched_tgss)
            #print len(idx_blob_matched_nvss), len(idx_blob_unmatched_nvss), len(idx_blob_unmatched_tgss) # DEBUG - check in aladin result of cross-match NVSS/TGSS catalogs
            classify(idx, idx_blob_matched_tgss, idx_blob_matched_nvss, idx_blob_unmatched_tgss, idx_blob_unmatched_nvss, t_tgss, t_nvss, 'U')

        def classify(idx, idx_matched_blob_this, idx_matched_blob_that, idx_unmatched_blob_this, idx_unmatched_blob_that, t_this, t_that, limit):
            # if in an island with one or more unmatched sources of same freq -> add as upper/lower limit?
            if len(idx_blob_unmatched_tgss) > 0 and len(idx_blob_unmatched_nvss) == 0 and len(idx_blob_matched_tgss) == 0:
                add_line(t, t_nvss, t_tgss, idx_nvss=None, idx_tgss=idx_blob_unmatched_tgss, srl_type=limit)
                [ idx_unmatched_tgss.pop(this_idx) for this_idx in idx_blob_unmatched_tgss if this_idx != idx ] # remove TGSS

            # if in an island with other unmatched sources of different freq -> combine+remove unmatch
            elif len(idx_blob_unmatched_tgss) > 0 and len(idx_blob_unmatched_nvss) > 0 and len(idx_blob_matched_tgss) == 0:
                add_line(t, t_nvss, t_tgss, idx_nvss=idx_blob_unmatched_nvss, idx_tgss=idx_blob_unmatched_tgss, srl_type='C')
                [ idx_unmatched_nvss.pop(this_idx) for this_idx in idx_blob_unmatched_nvss if this_idx != idx ] # remove NVSS
                [ idx_unmatched_tgss.pop(this_idx) for this_idx in idx_blob_unmatched_tgss if this_idx != idx ] # remove TGSS

            # if in an island with matched and/or umatched sources -> combine+remove match
            elif len(idx_blob_unmatched_tgss) > 0 and len(idx_blob_unmatched_nvss) >= 0 and len(idx_blob_matched_tgss) >= 0:
                add_line(t, t_nvss, t_tgss, idx_nvss=idx_blob_unmatched_nvss+idx_blob_matched_nvss, \
                        idx_tgss=idx_blob_unmatched_tgss+idx_blob_matched_tgss, srl_type='C')
                [ idx_unmatched_nvss.pop(this_idx) for this_idx in idx_blob_unmatched_nvss if this_idx != idx ] # remove NVSS
                [ idx_unmatched_tgss.pop(this_idx) for this_idx in idx_blob_unmatched_tgss if this_idx != idx ] # remove TGSS
                [ idx_matched_nvss.pop(this_idx) for this_idx in idx_blob_matched_nvss if this_idx != idx ] # remove NVSS (match)
                [ idx_matched_tgss.pop(this_idx) for this_idx in idx_blob_matched_tgss if this_idx != idx ] # remove TGSS (match)

            else:
                print "Something went wrong..."
                sys.exit()


        def add_line(t, t_nvss=None, t_tgss=None, idx_nvss=None, idx_tgss=None, srl_type='S'):
            """
            Add one or more lines to table t
            """
            idx = len(t)
            with file.open('idx_isl.txt', 'a') as f:
                f.write(str(idx)+' - NVSS'+' '.join(idx_nvss)+' - TGSS'+' '.join(idx_tgss))

            if idx_nvss != None and idx_tgss != None:
                # Check if sum S/N > 5

                # Source_id RA E_RA DEC E_DEC Total_flux_NVSS E_Total_flux_NVSS Peak_flux_NVSS E_Peak_flux_NVSS Total_flux_TGSS E_Total_flux_TGSS Peak_flux_TGSS E_Peak_flux_TGSS spidx E_spidx Rms_NVSS Rms_TGSS S_Code s2n size Mask
                t.add_row((idx,))

            elif idx_nvss != None:
                # Check if S/N > 5
                pass

            elif idx_tgss != None:
                # Check if S/N > 5
                pass

        add_line(t, t_nvss, t_tgss, idx_matched_nvss, idx_matched_tgss)
        
        print "Writing catalogue"
        t.write('spidx-cat.fits')
        
        sys.exit()

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
