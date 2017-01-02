#!/usr/bin/python

import os, sys, glob
import numpy as np
from scipy.ndimage.measurements import label
from scipy.ndimage.measurements import center_of_mass

from lofar import bdsm

import astropy.io.fits as pyfits
from astropy import wcs
from astropy.table import Table, vstack
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

source_id = 0
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

        # Creating new astropy Table
        t = Table(names=('Source_id','RA','E_RA','DEC','E_DEC','Total_flux_NVSS','E_Total_flux_NVSS','Peak_flux_NVSS','E_Peak_flux_NVSS','Total_flux_TGSS','E_Total_flux_TGSS','Peak_flux_TGSS','E_Peak_flux_TGSS','spidx','E_spidx','Rms_NVSS','Rms_TGSS','S_Code','Num_match','Num_unmatch','s2n','Mask'\
                      dtype=('i','f','f','f','f','f','f','f','f','f','S1','S100'))

        class t_surv():
            def __init__(self, cat, name):
                self.t = Table.read(cat_srl_nvss)
                self.len = len(self.t)
                self.name = name
                self.t['s2n'] = self.t['Total_flux']/self.t['E_Total_flux']
                self.t['s2n'].unit = ''
                x, y, _, _ = w.all_world2pix(self.t['RA'],self.t['DEC'],np.zeros(self.len),np.zeros(self.len), 0, ra_dec_order=True)
                self.x = x.round().astype(int)
                self.y = y.round().astype(int)

        nvss = t_surv(cat_srl_nvss, 'NVSS')
        tgss = t_surv(cat_srl_tgss, 'TGSS')

        # Cross-match NVSS-TGSS source catalogue
        # match_coordinates_sky() gives an idx of the 2nd catalogue for each source of the 1st catalogue
        idx_match, sep, _ = match_coordinates_sky(SkyCoord(nvss.t['RA'], nvss.t['DEC']),\
                                          SkyCoord(tgss.t['RA'], tgss.t['DEC']))

        match_sep = 15 # arcsec
        idx_matched_nvss = np.arange(0,nvss.len)[sep<=match_sep*u.arcsec] # nvss idx with a match
        idx_matched_tgss = idx_match[sep<=match_sep*u.arcsec] # tgss idx with a match
        idx_unmatched_nvss = np.arange(0,nvss.len)[sep>match_sep*u.arcsec] # nvss idx with NOT a match
        idx_unmatched_tgss = [ x for x in np.arange(0,tgss.len) if x not in idx_matched_tgss ] # tgss idx with NOT a match
        assert len(set(idx_matched_tgss)) == len(idx_matched_tgss) # check no sources are cross-matched twice
        print "Matched sources - NVSS: %i/%i - TGSS %i/%i" % ( len(idx_matched_nvss), nvss.len, len(idx_matched_tgss), tgss.len )

        w = wcs.WCS(pyfits.open(image_nvss)[0].header)
        x_nvss, y_nvss, _, _ = nvss.get_pix()
        x_tgss, y_tgss, _, _ = tgss.get_pix()
        val_blob_nvss = blobs[ np.zeros(nvss.len).astype(int), np.zeros(nvss.len).astype(int), nvss.y, nvss.x ] # blob val for ech source
        val_blob_tgss = blobs[ np.zeros(tgss.len).astype(int), np.zeros(tgss.len).astype(int), tgss.y, tgss.x ]
        assert (val_blob_nvss != 0).all() and (val_blob_tgss != 0).all() # all sources are into a blob           

        # For each island figure out what to do
        for blob in xrange(number_of_blobs):
            idx_blob_nvss = np.where(val_blob_nvss == blob) # idx of NVSS sources in this blob
            idx_blob_tgss = np.where(val_blob_tgss == blob) # idx of TGSS sources in this blob

            idx_blob_matched_nvss = np.intersect1d(idx_blob_nvss, idx_matched_nvss)
            idx_blob_matched_tgss = np.intersect1d(idx_blob_tgss, idx_matched_tgss)
            idx_blob_unmatched_nvss = np.intersect1d(idx_blob_nvss, idx_unmatched_nvss)
            idx_blob_unmatched_tgss = np.intersect1d(idx_blob_tgss, idx_unmatched_tgss)
            num_match = len(idx_blob_matched_nvss)+len(idx_blob_matched_tgss)
            num_unmatch = len(idx_blob_unmatched_nvss)+len(idx_blob_unmatched_tgss)
            #print len(idx_blob_matched_nvss), len(idx_blob_unmatched_nvss), len(idx_blob_unmatched_tgss) # DEBUG - check in aladin result of cross-match NVSS/TGSS catalogs

            # add all matched sources (multiple entries)
            if len(idx_blob_unmatched_nvss) == 0 and len(idx_blob_unmatched_tgss) == 0 and len(idx_blob_matched_nvss) == 1: code = 'S' # single
            elif len(idx_blob_unmatched_nvss) == 0 and len(idx_blob_unmatched_tgss) == 0 and len(idx_blob_matched_nvss) > 1: code = 'M' # multiple (all matched)
            else: code = 'C' # multiple (complex -> with unmatched sources in the same island)
            for i in xrange(len(idx_blob_matched_nvss)):
                s_nvss = nvss.t[idx_blob_matched_nvss[i]]
                s_tgss = tgss.t[idx_blob_matched_tgss[i]]
                s2n = s_nvss['Total_flux']/s_nvss['E_Total_flux'] + s_tgss['Total_flux']/s_tgss['E_Total_flux']
                if s2n < 5: 
                    print 'skip Good: s2n==%f' % s2n 
                    continue
                ra = np.average([s_nvss['RA'],s_tgss['RA']], weights=[s_nvss['E_RA']**2,s_tgss['E_RA']**2])
                e_ra = np.sqrt(s_nvss['E_RA']**2+s_tgss['E_RA']**2)
                dec = np.average([s_nvss['DEC'],s_tgss['DEC']], weights=[s_nvss['E_DEC']**2,s_tgss['E_DEC']**2])
                e_dec = np.sqrt(s_nvss['E_DEC']**2+s_tgss['E_DEC']**2)
                spidx, e_spidx = twopoint_spidx_bootstrap([147.,1400.], \
                        [s_tgss['Total_flux'],s_nvss['Total_flux']], [s_tgss['E_Total_flux'], s_nvss['E_Total_flux']], niter=1000)

                # Source_id RA E_RA DEC E_DEC Total_flux_NVSS E_Total_flux_NVSS Peak_flux_NVSS E_Peak_flux_NVSS Total_flux_TGSS E_Total_flux_TGSS Peak_flux_TGSS E_Peak_flux_TGSS spidx E_spidx Rms_NVSS Rms_TGSS S_Code Num_match Num_unmatch s2n Mask
                t.add_row((source_id, ra, e_ra, dec, e_dec, \
                        s_nvss['Total_flux'], s_nvss['E_Total_flux'], s_nvss['Peak_flux'], s_nvss['E_Peak_flux'],\
                        s_tgss['Total_flux'], s_tgss['E_Total_flux'], s_tgss['Peak_flux'], s_tgss['E_Peak_flux'],\
                        spidx, e_spidx, s_nvss['rms'], s_tgss['rms'], code, num_match, num_unmatch, s2n, mask))
                source_id += 1

            # if in an island with one or more unmatched NVSS sources -> add each as lower limit (multiple entries)
            if len(idx_blob_unmatched_nvss) > 0 and len(idx_blob_tgss) == 0:
                for i in xrange(len(idx_blob_unmatched_nvss)):
                    s_nvss = nvss.t[idx_blob_unmatched_nvss[i]]
                    s2n = np.sum(s_nvss['Total_flux']/s_nvss['E_Total_flux'])
                    if s2n < 5: 
                        print 'skip L: s2n==%f' % s2n 
                        continue
                    rms_tgss =  np.mean(np.array( pyfits.getdata(image_rms_tgss, 0) ) ) # no source, get from map
                    spidx = np.log10(s_nvss['Total_flux']/(3*rms_tgss))/np.log10(1400./147.)

                    t.add_row((source_id, s_nvss['RA'], s_nvss['E_RA'], s_nvss['DEC'], s_nvss['E_DEC'], \
                            s_nvss['Total_flux'], s_nvss['E_Total_flux'], s_nvss['Peak_flux'], s_nvss['E_Peak_flux'],\
                            0, 0, 0, 0,\
                            spidx, 0, s_nvss['rms'], rms_tgss, 'L', num_match, num_unmatch, s2n, mask))
                    source_id += 1

            # if in an island with one or more unmatched TGSS sources -> add each as upper limit (multiple entries)
            elif len(idx_blob_unmatched_tgss) > 0 and len(idx_blob_nvss) == 0:
                for i in xrange(len(idx_blob_unmatched_nvss)):
                    s_tgss = tgss.t[idx_blob_unmatched_tgss[i]]
                    s2n = np.sum(s_tgss['Total_flux']/s_tgss['E_Total_flux'])
                    if s2n < 5: 
                        print 'skip U: s2n==%f' % s2n 
                        continue
                    rms_nvss =  np.mean(np.array( pyfits.getdata(image_rms_nvss, 0) ) ) # no source, get from map
                    spidx = np.log10(s_tgss['Total_flux']/(3*rms_nvss))/np.log10(147./1400.)

                    t.add_row((source_id, s_tgss['RA'], s_tgss['E_RA'], s_tgss['DEC'], s_tgss['E_DEC'], \
                            0, 0, 0, 0,\
                            s_tgss['Total_flux'], s_tgss['E_Total_flux'], s_tgss['Peak_flux'], s_tgss['E_Peak_flux'],\
                            spidx, 0, s_tgss['rms'], rms_nvss, 'U', num_match, num_unmatch, s2n, mask))
                    source_id += 1

            # any other case with unmatched sources -> combine island (1 entry)
            elif len(idx_blob_unmatched_nvss) > 0 or len(idx_blob_unmatched_tgss) > 0:
                s_nvss = nvss.t[idx_blob_matched_nvss+idx_blob_unmatched_nvss]
                s_tgss = tgss.t[idx_blob_matched_tgss+idx_blob_unmatched_tgss]
                s_stack = vstack([s_nvss, s_tgss])
                s2n = np.sum(s_nvss['Total_flux'])/np.sqrt(np.sum(s_nvss['E_Total_flux']**2)) + \
                      np.sum(s_tgss['Total_flux'])/np.sqrt(np.sum(s_tgss['E_Total_flux']**2))

                if s2n < 5: 
                    print 'skip Good: s2n==%f' % s2n 
                    continue
                ra = np.average(s_stack['RA'], weights=s_stack['E_RA']**2)
                e_ra = np.sqrt(np.sum(s_stack['E_RA']**2))
                dec = np.average(s_stack['DEC'], weights=s_stack['E_DEC']**2)
                e_dec = np.sqrt(np.sum(s_stack['E_DEC']**2))
                spidx, e_spidx = twopoint_spidx_bootstrap([1400.,147.], \
                        [np.sum(s_nvss['Total_flux']), np.sum(s_tgss['Total_flux'])], \
                        [np.sqrt(np.sum(s_nvss['E_Total_flux']**2)), np.sqrt(np.sum(s_tgss['E_Total_flux']**2))], niter=1000)

                t.add_row((source_id, ra, e_ra, dec, e_dec, \
                        np.sum(s_nvss['Total_flux']), np.sqrt(np.sum(s_nvss['E_Total_flux']**2)), np.max(s_nvss['Peak_flux']), s_nvss['E_Peak_flux'][np.argmax(s_nvss['Peak_flux'])],\
                        np.sum(s_tgss['Total_flux']), np.sqrt(np.sum(s_tgss['E_Total_flux']**2)), np.max(s_tgss['Peak_flux']), s_tgss['E_Peak_flux'][np.argmax(s_tgss['Peak_flux'])],\
                        spidx, e_spidx, np.mean(s_nvss['rms']), np.mean(s_tgss['rms']), 'I', num_match, num_unmatch, s2n, mask))
                source_id += 1


        # adding a raw require the entire cat to be rewritten in memory
        # to speed up a catalogue is generated for each mask and appended to the existing file
        print "Writing concat catalogue"
        t_complete = Table.read('spidx-cat.fits')
        vstack([t_complete, t]).write('spidx-cat.fits')
        t_complete.close()
        
        sys.exit()

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

        data_nvss = np.array( pyfits.getdata(image_nvss, 0) )
        data_rms_nvss =  np.array( pyfits.getdata(image_rms_nvss, 0) )
        data_tgss = np.array( pyfits.getdata(image_tgss, 0) )
        data_rms_tgss =  np.array( pyfits.getdata(image_rms_tgss, 0) )

        snr_nvss = data_nvss/data_rms_nvss
        snr_nvss[snr_nvss<0] = 0
        snr_tgss = data_tgss/data_rms_tgss
        snr_tgss[snr_tgss<0] = 0
        mask = pyfits.getdata(image_mask, 0)
        idx_g = (snr_nvss > 3) & (snr_tgss > 3) & (mask == 1)
        idx_u = (snr_tgss > 3) & (snr_nvss < 3) & (mask == 1)
        idx_l = (snr_nvss > 3) & (snr_tgss < 3) & (mask == 1)
        print "Number of good pixels: ", idx_g.sum()
        print "Number of upper lim: ", idx_u.sum()
        print "Number of lower lim: ", idx_l.sum()

        data_nvss[idx_u] = 3*data_rms_nvss[idx_u] # set val for upper limits
        data_tgss[idx_l] = 3*data_rms_nvss[idx_l] # set val for lower limits

        #data_spidx_g = np.zeros(shape=data_nvss[idx_g].shape)
        #data_spidx_g_err = np.zeros(shape=data_nvss[idx_g].shape)

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
