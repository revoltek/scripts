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

from lib_linearfit import twopoint_spidx_bootstrap

import warnings
warnings.filterwarnings("ignore") 

# options
#wdir = '/home/fdg/phd/works/spidxskymap/sw_test/'
wdir = '/data2/spidxskymap/'
#area = 10.1978092553 # beam area in pixels - with beam: 45"x45" and pixels size: 15"x15"
area = 0.0125*0.0125*np.pi/(4*np.log(2.)) # beam area in deg

def remove_duplicates(file_cat='spidx-cat.txt'):
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
    centers = Table([[],[],[]], names=('Mask','RA','DEC'), dtype=['S100',float,float])
    centers['RA'].unit = 'deg'
    centers['DEC'].unit = 'deg'
    for i, mask_file in enumerate(glob.glob('masks/mask*.fits')):
        with pyfits.open(mask_file) as fits:
            head = fits[0].header
            ra = head['CRVAL1']
            dec = head['CRVAL2']
            centers.add_row([os.path.basename(mask_file), ra, dec])

    print "Matching catalogues..."
    sources = Table.read('spidx-cat.fits', format='fits')
    idx, _, _ = match_coordinates_sky(SkyCoord(sources['RA']*u.deg, sources['DEC']*u.deg),\
                                      SkyCoord(centers['RA'], centers['DEC']))

    print "Removing duplicates..."
    idx_duplicates = [] 
    for i, source in enumerate(sources):
        # check if closest
#        print idx[i], source['Mask'], centers[int(idx[i])]['Mask']
        if source['Mask'] != centers[int(idx[i])]['Mask']:
            idx_duplicates.append(i)
#            print "Removing source ", i
    print "Removing a total of", len(idx_duplicates), "sources over", len(sources)
    sources.remove_rows(idx_duplicates)

    print "Add unique Isl_id..."
    # add unique island idx based on combination of Mask name and blob idx
    last_isl_id = 0
    for mask in set(sources['Mask']):
        # first cycle add 0 to blob_id, next cycles add highest isl_id from previous cycle (blob_ids are >0)
        incr = np.max(sources['Isl_id'][ np.where(sources['Mask'] == mask) ])
        sources['Isl_id'][ np.where(sources['Mask'] == mask) ] += last_isl_id
        last_isl_id += incr

    sources.remove_columns('Mask')
    sources['Source_id'] = range(len(sources)) # set id after removal
    sources.write('spidx-cat-nodup.fits', format='fits', overwrite=True)


def clip_gaus(img, rmsimg, nsigma=1):
    """
    Use local rms from rmsimg and clip fitsimage
    to nsigma times this value. Then return the mask.
    """
    with pyfits.open(rmsimg) as fits:
        rms_data = fits[0].data

    with pyfits.open(img) as fits:
        data = fits[0].data
        mask = (data > nsigma*rms_data).astype(int)
        fits[0].data = mask
        fits.writeto(img, clobber=True)
    
    return img


def find_closest_isl(x, y, blobs):
    """
    find closest non-0 pixel in blobs to (x,y)
    can be speed up removing cycles to create dists
    """
    x = float(x)
    y = float(y)
    dists = np.zeros(shape=blobs.shape)
    for j in xrange(blobs.shape[0]):
        for i in xrange(blobs.shape[1]):
            if blobs[j,i] == 0:
                dists[j,i] = np.inf
            else: 
                dists[j,i] = np.sqrt( (i-x)**2 + (j-y)**2 )
    blob_y, blob_x = np.unravel_index(dists.argmin(), dists.shape)
    assert dists[blob_y, blob_x] < 5
    print " - associated with pix:", blob_x, blob_y
    return blobs[blob_y, blob_x]



class t_surv():
    """
    A class to store survey information
    """
    def __init__(self, cat):
        self.w = wcs.WCS(pyfits.open(image_nvss)[0].header)
        self.t = Table.read(cat)

        # arbitrary remove bad detections, it happens a couple of times
        self.t.remove_rows(np.where(self.t['E_Total_flux'] == 0.)[0])

        self.len = len(self.t)
        #self.t['s2n'] = self.t['Total_flux']/self.t['E_Total_flux']
        #self.t['s2n'].unit = ''
        x, y, _, _ = self.w.all_world2pix(self.t['RA'],self.t['DEC'],np.zeros(self.len),np.zeros(self.len), 0, ra_dec_order=True)
        self.x = x.round().astype(int)
        self.y = y.round().astype(int)


# assume NVSS/TGSS files have same name apart from initial part
images_nvss = glob.glob(wdir+'NVSS/*fits')

for image_nvss in images_nvss:
    image_tgss = image_nvss.replace('NVSS','TGSS')
    if not os.path.exists(image_tgss):
        print "TGSS file: ", image_tgss, "does not exist, continue."
        continue

    print "-- Make catlog on:", image_nvss, image_tgss

    # source finder
    image_rms_nvss = image_nvss.replace('.fits','-rms.fits').replace('NVSS','NVSS/rms',1)
    #image_gaus_nvss = image_nvss.replace('.fits','-gaus.fits').replace('NVSS','NVSS/gaus',1)
    image_isl_nvss = image_nvss.replace('.fits','-isl.fits').replace('NVSS','NVSS/isl',1)
    cat_srl_nvss = image_nvss.replace('.fits','-srl.fits').replace('NVSS','NVSS/catalog',1)
    if not os.path.exists(cat_srl_nvss) or not os.path.exists(image_isl_nvss) or not os.path.exists(image_rms_nvss):# or not os.path.exists(image_gaus_nvss):
        c = bdsm.process_image(image_nvss, frequency=1400e6, rms_box=(102,34), advanced_opts=True, group_tol=0.5, thresh_isl=3, thresh_pix=4)
        c.export_image(outfile=image_rms_nvss, img_type='rms', clobber=True)
        #c.export_image(outfile=image_gaus_nvss, img_type='gaus_model', clobber=True)
        c.export_image(outfile=image_isl_nvss, img_type='island_mask', clobber=True)
        c.write_catalog(outfile=cat_srl_nvss, catalog_type='srl', format='fits', clobber=True)
        # clip the gaus images and create island masks
        # we use gaus images because they catch better the extended emission than isl_mask
        #clip_gaus(image_gaus_nvss, image_rms_nvss)

    image_rms_tgss = image_tgss.replace('.fits','-rms.fits').replace('TGSS','TGSS/rms',1)
    #image_gaus_tgss = image_tgss.replace('.fits','-gaus.fits').replace('TGSS','TGSS/gaus',1)
    image_isl_tgss = image_tgss.replace('.fits','-isl.fits').replace('TGSS','TGSS/isl',1)
    cat_srl_tgss = image_tgss.replace('.fits','-srl.fits').replace('TGSS','TGSS/catalog',1)
    if not os.path.exists(cat_srl_tgss) or not os.path.exists(image_isl_tgss) or not os.path.exists(image_rms_tgss):# or not os.path.exists(image_gaus_tgss):
        c = bdsm.process_image(image_tgss, frequency=147e6, rms_box=(102,34), advanced_opts=True, group_tol=0.5, thresh_isl=3, thresh_pix=4)
        c.export_image(outfile=image_rms_tgss, img_type='rms', clobber=True)
        #c.export_image(outfile=image_gaus_tgss, img_type='gaus_model', clobber=True)
        c.export_image(outfile=image_isl_tgss, img_type='island_mask', clobber=True)
        c.write_catalog(outfile=cat_srl_tgss, catalog_type='srl', format='fits', clobber=True)
        #clip_gaus(image_gaus_tgss, image_rms_tgss)

    image_mask = wdir+'masks/'+os.path.basename(image_nvss).replace('NVSS_','mask_')
    if not os.path.exists(image_mask):
        
        # create a mask with 1 where at least one survey has a detection, 0 otherwise
        print "Making mask:", image_mask
        with pyfits.open(image_isl_nvss) as fits_nvss:
            with pyfits.open(image_isl_tgss) as fits_tgss:
                fits_nvss[0].data = ((fits_nvss[0].data == 1) | (fits_tgss[0].data == 1)).astype(float)
                fits_nvss.writeto(image_mask, clobber=True)

        # Separate combined island mask into islands, each has a number
        with pyfits.open(image_mask) as fits_mask:
            blobs, number_of_blobs = label(fits_mask[0].data.astype(int).squeeze(), structure=[[1,1,1],[1,1,1],[1,1,1]])
        print "# of islands found:", number_of_blobs

        # Creating new astropy Table
        t = Table(names=('RA','DEC','Total_flux_NVSS','E_Total_flux_NVSS','Peak_flux_NVSS','E_Peak_flux_NVSS','Rms_NVSS','Total_flux_TGSS','E_Total_flux_TGSS','Peak_flux_TGSS','E_Peak_flux_TGSS','Rms_TGSS','Spidx','E_Spidx','s2n','S_code','Num_match','Num_unmatch_NVSS','Num_unmatch_TGSS','Isl_id','Mask'),\
                  dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','S1','i4','i4','i4','i4','S100'))

        nvss = t_surv(cat_srl_nvss)
        tgss = t_surv(cat_srl_tgss)

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

        val_blob_nvss = blobs[ nvss.y, nvss.x ] # blob val for each source
        val_blob_tgss = blobs[ tgss.y, tgss.x ]
        
        # if source slightly outside island: find closest island
        for s in np.where(val_blob_nvss == 0)[0]:
            print "NVSS ORPHAN SOURCE:", image_mask, nvss.x[s], nvss.y[s],
            val_blob_nvss[s] = find_closest_isl(nvss.x[s], nvss.y[s], blobs)

        for s in np.where(val_blob_tgss == 0)[0]:
            print "TGSS ORPHAN SOURCE:", image_mask, tgss.x[s], tgss.y[s],
            val_blob_tgss[s] = find_closest_isl(tgss.x[s], tgss.y[s], blobs)

        assert (val_blob_nvss != 0).all() and (val_blob_tgss != 0).all() # all sources are into a blob           

        # For each island figure out what to do - number 0 is no-island
        for blob in xrange(1,number_of_blobs):
            idx_blob_nvss = np.where(val_blob_nvss == blob) # idx of NVSS sources in this blob
            idx_blob_tgss = np.where(val_blob_tgss == blob) # idx of TGSS sources in this blob

            idx_blob_matched_nvss = np.intersect1d(idx_blob_nvss, idx_matched_nvss)
            idx_blob_matched_tgss = idx_match[idx_blob_matched_nvss] # this preserves order within a blob (important only for matched sources in the M case)
            idx_blob_unmatched_nvss = np.intersect1d(idx_blob_nvss, idx_unmatched_nvss)
            idx_blob_unmatched_tgss = np.intersect1d(idx_blob_tgss, idx_unmatched_tgss)
            num_match = len(idx_blob_matched_nvss) # same of len(idx_blob_matched_tgss)
            num_unmatch_NVSS = len(idx_blob_unmatched_nvss)
            num_unmatch_TGSS = len(idx_blob_unmatched_tgss)
            
            # add all matched sources (multiple entries)
            if num_unmatch_NVSS == 0 and num_unmatch_TGSS == 0 and num_match == 1: code = 'S' # single
            elif num_unmatch_NVSS == 0 and num_unmatch_TGSS == 0 and num_match > 1: code = 'M' # multiple (all matched)
            else: code = 'C' # multiple (complex -> with unmatched sources in the same island)
            for i in xrange(num_match):
                s_nvss = nvss.t[idx_blob_matched_nvss[i]]
                s_tgss = tgss.t[idx_blob_matched_tgss[i]]
                s2n = np.sqrt( (s_nvss['Peak_flux']/s_nvss['Isl_rms'])**2 + (s_tgss['Peak_flux']/s_tgss['Isl_rms'])**2 )
                if s2n < 5: 
                    print 'skip Good: s2n==%f' % s2n 
                    continue

                ra = np.average([s_nvss['RA'],s_tgss['RA']], weights=[1/s_nvss['E_RA']**2,1/s_tgss['E_RA']**2])
                dec = np.average([s_nvss['DEC'],s_tgss['DEC']], weights=[1/s_nvss['E_DEC']**2,1/s_tgss['E_DEC']**2])
                spidx, e_spidx = twopoint_spidx_bootstrap([147.,1400.], \
                        [s_tgss['Total_flux'],s_nvss['Total_flux']], [s_tgss['E_Total_flux'], s_nvss['E_Total_flux']], niter=1000)

                t.add_row((ra, dec, \
                        s_nvss['Total_flux'], s_nvss['E_Total_flux'], s_nvss['Peak_flux'], s_nvss['E_Peak_flux'], s_nvss['Isl_rms'], \
                        s_tgss['Total_flux'], s_tgss['E_Total_flux'], s_tgss['Peak_flux'], s_tgss['E_Peak_flux'], s_tgss['Isl_rms'], \
                        spidx, e_spidx, s2n, code, num_match, num_unmatch_NVSS, num_unmatch_TGSS, blob, os.path.basename(image_mask)))

            # if in an island with one or more unmatched NVSS sources -> add each as lower limit (multiple entries)
            if num_unmatch_NVSS > 0 and len(idx_blob_tgss[0]) == 0:
                code = 'L'
                for i in xrange(num_unmatch_NVSS):
                    s_nvss = nvss.t[idx_blob_unmatched_nvss[i]]
                    s2n = np.sum(s_nvss['Peak_flux']/s_nvss['Isl_rms'])
                    if s2n < 5: 
                    #    print 'skip L: s2n==%f' % s2n 
                        continue
                    x, y, _, _ = nvss.w.all_world2pix(s_nvss['RA'],s_nvss['DEC'],0,0, 0, ra_dec_order=True)
                    rms_tgss = pyfits.getdata(image_rms_tgss, 0)[0,0,int(np.around(y)),int(np.around(x))] # no source, get from map

                    # correction for extended (estimated using island size) sources
                    ext_corr = np.sqrt( (s_nvss['Maj']*s_nvss['Min']*np.pi/(4*np.log(2.)))/area )
                    if ext_corr > 1: rms_tgss *= ext_corr
                    
                    # for few sources on the edges
                    if np.isnan(rms_tgss):
                        continue
                    spidx = np.log10(s_nvss['Total_flux']/(3*rms_tgss))/np.log10(1400./147.)

                    t.add_row((s_nvss['RA'], s_nvss['DEC'], \
                            s_nvss['Total_flux'], s_nvss['E_Total_flux'], s_nvss['Peak_flux'], s_nvss['E_Peak_flux'], s_nvss['Isl_rms'], \
                            0, 0, 0, 0, rms_tgss, \
                            spidx, 0, s2n, code, num_match, num_unmatch_NVSS, num_unmatch_TGSS, blob, os.path.basename(image_mask)))

            # if in an island with one or more unmatched TGSS sources -> add each as upper limit (multiple entries)
            elif num_unmatch_TGSS > 0 and len(idx_blob_nvss[0]) == 0:
                code = 'U'
                for i in xrange(num_unmatch_TGSS):
                    s_tgss = tgss.t[idx_blob_unmatched_tgss[i]]
                    s2n = np.sum(s_tgss['Peak_flux']/s_tgss['Isl_rms'])
                    if s2n < 5: 
                    #    print 'skip U: s2n==%f' % s2n 
                        continue
                    x, y, _, _ = tgss.w.all_world2pix(s_tgss['RA'],s_tgss['DEC'],0,0, 0, ra_dec_order=True)
                    rms_nvss = pyfits.getdata(image_rms_nvss, 0)[0,0,int(np.around(y)),int(np.around(x))] # no source, get from map

                    # correction for extended (estimated using island size) sources
                    ext_corr = np.sqrt( (s_tgss['Maj']*s_tgss['Min']*np.pi/(4*np.log(2.)))/area )
                    if ext_corr > 1: rms_nvss *= ext_corr

                    # for few sources on the edges
                    if np.isnan(rms_nvss):
                        continue
                    spidx = np.log10(s_tgss['Total_flux']/(3*rms_nvss))/np.log10(147./1400.)

                    t.add_row((s_tgss['RA'], s_tgss['DEC'], \
                            0, 0, 0, 0, rms_nvss, \
                            s_tgss['Total_flux'], s_tgss['E_Total_flux'], s_tgss['Peak_flux'], s_tgss['E_Peak_flux'], s_tgss['Isl_rms'], \
                            spidx, 0, s2n, code, num_match, num_unmatch_NVSS, num_unmatch_TGSS, blob, os.path.basename(image_mask)))

            # any other case with unmatched sources -> combine island (1 entry)
            elif num_unmatch_NVSS > 0 or num_unmatch_TGSS > 0:
                code = 'I'
                s_nvss = vstack([nvss.t[idx_blob_matched_nvss], nvss.t[idx_blob_unmatched_nvss]])
                s_tgss = vstack([tgss.t[idx_blob_matched_tgss], tgss.t[idx_blob_unmatched_tgss]])
                s_stack = vstack([s_nvss, s_tgss])
                s2n = np.sum(s_nvss['Peak_flux'])/np.sqrt(np.sum(s_nvss['Isl_rms']**2)) + \
                      np.sum(s_tgss['Peak_flux'])/np.sqrt(np.sum(s_tgss['Isl_rms']**2))

                if s2n < 5: 
                #    print 'skip I: s2n==%f' % s2n 
                    continue
                ra = np.average(s_stack['RA'], weights=1/s_stack['E_RA']**2)
                dec = np.average(s_stack['DEC'], weights=1/s_stack['E_DEC']**2)
                spidx, e_spidx = twopoint_spidx_bootstrap([1400.,147.], \
                        [np.sum(s_nvss['Total_flux']), np.sum(s_tgss['Total_flux'])], \
                        [np.sqrt(np.sum(s_nvss['E_Total_flux']**2)), np.sqrt(np.sum(s_tgss['E_Total_flux']**2))], niter=1000)

                t.add_row((ra, dec, \
                        np.sum(s_nvss['Total_flux']), np.sqrt(np.sum(s_nvss['E_Total_flux']**2)), np.max(s_nvss['Peak_flux']), s_nvss['E_Peak_flux'][np.argmax(s_nvss['Peak_flux'])], np.mean(s_nvss['Isl_rms']), \
                        np.sum(s_tgss['Total_flux']), np.sqrt(np.sum(s_tgss['E_Total_flux']**2)), np.max(s_tgss['Peak_flux']), s_tgss['E_Peak_flux'][np.argmax(s_tgss['Peak_flux'])], np.mean(s_tgss['Isl_rms']), \
                        spidx, e_spidx, s2n, code, num_match, num_unmatch_NVSS, num_unmatch_TGSS, blob, os.path.basename(image_mask)))

            # DEBUG - check in aladin result of cross-match NVSS/TGSS catalogs
            #print "Match NVSS:", len(idx_blob_matched_nvss), \
            #      "Match TGSS:", len(idx_blob_matched_tgss), \
            #      "Un-Match NVSS:", len(idx_blob_unmatched_nvss), \
            #      "Un-Match TGSS:", len(idx_blob_unmatched_tgss), \
            #      "- code:", code

        # adding a raw require the entire cat to be rewritten in memory
        # to speed up a catalogue is generated for each mask and appended to the existing file
        print "Writing catalogue"
        if os.path.exists('spidx-cat.fits'):
            t_complete = Table.read('spidx-cat.fits')
            # remove existing rows of this mask (important for re-runs)
            idx_remove = np.where(t_complete['Mask'] == os.path.basename(image_mask))[0]
            if len(idx_remove) > 0:
                print "Deleting %i rows for mask %s." % (len(idx_remove), os.path.basename(image_mask))
                t_complete.remove_rows(idx_remove)
            vstack([t_complete, t]).write('spidx-cat.fits', overwrite=True)
            del t_complete
        else:
            t.write('spidx-cat.fits')

remove_duplicates()

for image_nvss in images_nvss:
    image_tgss = image_nvss.replace('NVSS','TGSS')
    if not os.path.exists(image_tgss):
        print "TGSS file: ", image_tgss, "does not exist, continue."
        continue

    print "-- Make spidx map of:", image_nvss, image_tgss

    image_mask = wdir+'masks/'+os.path.basename(image_nvss).replace('NVSS_','mask_')
    image_rms_nvss = image_nvss.replace('.fits','-rms.fits').replace('NVSS','NVSS/rms',1)
    image_rms_tgss = image_tgss.replace('.fits','-rms.fits').replace('TGSS','TGSS/rms',1)
    # new names
    image_spidx = wdir+'spidx/'+os.path.basename(image_nvss).replace('NVSS_','spidx_')
    image_spidx_u = wdir+'spidx_u/'+os.path.basename(image_nvss).replace('NVSS_','spidxU_')
    image_spidx_l = wdir+'spidx_l/'+os.path.basename(image_nvss).replace('NVSS_','spidxL_')
    image_spidx_err = wdir+'spidx_err/'+os.path.basename(image_nvss).replace('NVSS_','spidx_').replace('.fits','-err.fits')
    if not os.path.exists(image_spidx) or not os.path.exists(image_spidx_err) or not os.path.exists(image_spidx_u) or not os.path.exists(image_spidx_l):
        print "Makign spidx map..."

        data_nvss = np.array( pyfits.getdata(image_nvss, 0) ).squeeze()
        data_rms_nvss =  np.array( pyfits.getdata(image_rms_nvss, 0) ).squeeze()
        data_tgss = np.array( pyfits.getdata(image_tgss, 0) ).squeeze()
        data_rms_tgss =  np.array( pyfits.getdata(image_rms_tgss, 0) ).squeeze()

        snr_nvss = data_nvss/data_rms_nvss
        snr_nvss[snr_nvss<0] = 0
        snr_tgss = data_tgss/data_rms_tgss
        snr_tgss[snr_tgss<0] = 0
        mask = pyfits.getdata(image_mask, 0).squeeze()
        idx_g = (snr_nvss > 3) & (snr_tgss > 3) & (mask == 1)
        idx_u = (snr_tgss > 3) & (snr_nvss < 3) & (mask == 1)
        idx_l = (snr_nvss > 3) & (snr_tgss < 3) & (mask == 1)
        print "Number of good pixels: ", idx_g.sum()
        print "Number of upper lim: ", idx_u.sum()
        print "Number of lower lim: ", idx_l.sum()

        data_spidx_g, data_spidx_g_err = twopoint_spidx_bootstrap([147.,1400.], [data_tgss[idx_g],data_nvss[idx_g]], \
                    [data_rms_tgss[idx_g],data_rms_nvss[idx_g]], niter=1000)
        # upper limits
        data_spidx_u = np.log10(data_tgss[idx_u]/(3*data_rms_nvss[idx_u]))/np.log10(147./1400.)
        # lower limits
        data_spidx_l = np.log10((3*data_rms_tgss[idx_l])/data_nvss[idx_l])/np.log10(147./1400.)

        # write spidx, spidx+upper/lower limits and spidx error map
        with pyfits.open(image_nvss) as fits_nvss:
            data = fits_nvss[0].data.squeeze() # drops the size-1 axes
            header = fits_nvss[0].header
            new_wcs = wcs.WCS(fits_nvss[0].header).celestial
            new_header = new_wcs.to_header()

            data[:] = np.nan
            data[idx_g] = data_spidx_g
            new_fh = pyfits.PrimaryHDU(data=data, header=new_header)
            new_fh.writeto(image_spidx, clobber=True)

            data[:] = np.nan
            data[idx_u] = data_spidx_u
            new_fh = pyfits.PrimaryHDU(data=data, header=new_header)
            new_fh.writeto(image_spidx_u, clobber=True)


            data[:] = np.nan
            data[idx_l] = data_spidx_l
            new_fh = pyfits.PrimaryHDU(data=data, header=new_header)
            new_fh.writeto(image_spidx_l, clobber=True)

            data[:] = np.nan
            data[idx_g] = data_spidx_g_err
            new_fh = pyfits.PrimaryHDU(data=data, header=new_header)
            new_fh.writeto(image_spidx_err, clobber=True)

print "All done!"
