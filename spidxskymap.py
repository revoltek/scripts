#!/usr/bin/python

import os, sys, glob
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
import lofar.bdsm as bdsm
from linearfit import linear_fit_bootstrap


def add_to_coordfile(coord_file, fits_file):
    header = pyfits.getheader(fits_file, 0)
    ra = header['CRVAL1']
    dec = header['CRVAL2']
    print ra, dec
    needle = fits_file+' '+str(ra)+' '+str(dec)+'\n'
    with open(coord_file, 'a+') as file:
        for line in file:
            if needle in line:
                break
        else: # not found, we are at the eof
            file.write(needle) # append missing data

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

# create file catalogues
print "Read NVSS coords."
if not os.path.exists('NVSS_coords.txt'):
    # read all NVSS and write positions
    nvss_files = glob.glob('NVSS/*fits')
    for nvss_file in nvss_files:
        add_to_coordfile('NVSS_coords.txt', nvss_file)
print "Read TGSS coords."
if not os.path.exists('TGSS_coords.txt'):
    # read all TGSS and write positions
    tgss_files = glob.glob('TGSS/*fits')
    for tgss_file in tgss_files:
        add_to_coordfile('TGSS_coords.txt', tgss_file)

# create pairs
print "Create Pairs."
types = np.dtype({'names':['name', 'ra', 'dec'], 'formats':['S100',np.float,np.float]})
coord_nvss = np.loadtxt('NVSS_coords.txt', comments='#', delimiter=' ', unpack=False, converters={}, dtype=types)
coord_tgss = np.loadtxt('TGSS_coords.txt', comments='#', delimiter=' ', unpack=False, converters={}, dtype=types)

cat_nvss = SkyCoord(coord_nvss['ra'], coord_nvss['dec'], frame="fk5", unit="deg")
cat_tgss = SkyCoord(coord_tgss['ra'], coord_tgss['dec'], frame="fk5", unit="deg")
# index is a list of index for the closest cat_nvss position of each cat_tgss entry
index, _, _ = cat_nvss.match_to_catalog_sky(cat_tgss)

# for each pair do bdsm
for nvss, tgss in zip(coord_nvss, coord_tgss[index]):
    print "Work on", nvss, tgss

    # source finder
    nvss_image_rms = nvss['name'].replace('.fits','-rms.fits')
    nvss_image_gaus = nvss['name'].replace('.fits','-gaus.fits')
    if not os.path.exists(nvss_image_rms) or not os.path.exists(nvss_image_gaus):
        c = bdsm.process_image(nvss['name'], frequency=1400e6, group_by_isl=True, mean_map='zero')
        c.export_image(outfile=nvss_image_rms, img_type='rms', clobber=True)
        c.export_image(outfile=nvss_image_gaus, img_type='gaus_model', clobber=True)

    tgss_image_rms = tgss['name'].replace('.fits','-rms.fits')
    tgss_image_gaus = tgss['name'].replace('.fits','-gaus.fits')
    if not os.path.exists(tgss_image_rms) or not os.path.exists(tgss_image_gaus):
        c = bdsm.process_image(tgss['name'], frequency=147e6, group_by_isl=True, mean_map='zero')
        c.export_image(outfile=tgss_image_rms, img_type='rms', clobber=True)
        c.export_image(outfile=tgss_image_gaus, img_type='gaus_model', clobber=True)

    # clip the gaus images and create island masks
    # we use gaus images because they catch better the extended emission than isl_mask
    clip_gaus(nvss_image_gaus, nvss_image_rms)
    clip_gaus(tgss_image_gaus, tgss_image_rms)

    # do an OR between the two island masks
    image_mask = 'masks/'+os.path.basename(nvss_image_gaus).replace('NVSS_','mask_')
    if not os.path.exists(image_mask):
        nvss_fits = pyfits.open(nvss_image_gaus)
        tgss_fits = pyfits.open(tgss_image_gaus)
        data = ((nvss_fits[0].data == 1) | (tgss_fits[0].data == 1)).astype(float)
        nvss_fits[0].data = data
        nvss_fits.writeto(image_mask, clobber=True)

    # spectral index map
    print "Makign spidx map..."
    image_spidx = 'spidx/'+os.path.basename(nvss['name']).replace('NVSS_','spidx_')
    image_spidx_err = 'spidx/'+os.path.basename(nvss['name']).replace('NVSS_','spidx_').replace('.fits','-err.fits')
    if not os.path.exists(image_spidx) or not os.path.exists(image_spidx_err):
        idx = (pyfits.getdata(nvss['name'], 0) > 3 * pyfits.getdata(nvss_image_rms, 0)) | (pyfits.getdata(tgss['name'], 0) > 3 * pyfits.getdata(tgss_image_rms, 0))
        idx_g = (pyfits.getdata(nvss['name'], 0) > 3 * pyfits.getdata(nvss_image_rms, 0)) & (pyfits.getdata(tgss['name'], 0) > 3 * pyfits.getdata(tgss_image_rms, 0))
        idx_u = (pyfits.getdata(nvss['name'], 0) <= 3 * pyfits.getdata(nvss_image_rms, 0)) & (pyfits.getdata(tgss['name'], 0) > 3 * pyfits.getdata(tgss_image_rms, 0))
        idx_l = (pyfits.getdata(nvss['name'], 0) > 3 * pyfits.getdata(nvss_image_rms, 0)) & (pyfits.getdata(tgss['name'], 0) <= 3 * pyfits.getdata(tgss_image_rms, 0))

        data_nvss = np.array( pyfits.getdata(nvss['name'], 0) )
        data_nvss_rms =  np.array( pyfits.getdata(nvss_image_rms, 0) )
        data_nvss[idx_u] = 3*data_nvss_rms[idx_u] # set val for upper limits

        data_tgss = np.array( pyfits.getdata(tgss['name'], 0) )
        data_tgss_rms =  np.array( pyfits.getdata(tgss_image_rms, 0) )
        data_tgss[idx_l] = 3*data_nvss_rms[idx_l] # set val for lower limits

        # to be quick update only those pixels that will eventually be used 
        data_nvss_rms[idx] = 0.434*data_nvss_rms[idx]/data_nvss[idx]
        data_tgss_rms[idx] = 0.434*data_tgss_rms[idx]/data_tgss[idx]
        data_nvss[idx] = np.log10( data_nvss[idx] )
        data_tgss[idx] = np.log10( data_tgss[idx] )

        data_spidx_g = np.zeros(shape=data_nvss[idx_g].shape)
        data_spidx_g_err = np.zeros(shape=data_nvss[idx_g].shape)
        for d in xrange(len(data_nvss[idx_g])):
            if i % 10 == 0: print '.',
            data_spidx_g[d], _, data_spidx_g_err[d], _ = linear_fit_bootstrap(np.log10([147.,1400.]), \
                    [data_tgss[idx_g][d],data_nvss[idx_g][d]], [data_tgss_rms[idx_g][d],data_nvss_rms[idx_g][d]], niter=100)
        # upper limits
        data_spidx_u = (data_tgss[idx_u]-data_nvss[idx_u])/np.log10(147./1400.)
        # lower limits
        data_spidx_l = (data_tgss[idx_l]-data_nvss[idx_l])/np.log10(147./1400.)

        # write spidx and error map
        nvss_fits = pyfits.open(nvss['name'])
        nvss_fits[0].data[:] = np.nan
        nvss_fits[0].data[idx_g] = data_spidx_g
        nvss_fits[0].data[idx_u] = data_spidx_u
        nvss_fits[0].data[idx_l] = data_spidx_l
        nvss_fits.writeto(image_spidx, clobber=True)
        nvss_fits[0].data[:] = np.nan
        nvss_fits[0].data[idx_g] = data_spidx_g_err
        nvss_fits[0].data[idx_u][:] = -1
        nvss_fits[0].data[idx_l][:] = -2
        nvss_fits.writeto(image_spidx_err, clobber=True)
        nvss_fits.close()

    # re-run the source finder
    print "Making catalogue..."
    c = bdsm.process_image(nvss['name'], frequency=1400e6, group_by_isl=True, detection_image=image_mask, beam=(1.25000e-02, 1.25000e-02, 0.0), mean_map='zero', rms_map=False, rms_value=0.1)
    cat = nvss['name'].replace('.fits','.cat')
    if c.sources == []:
        open(cat, 'a').close()
    else:
        c.write_catalog(outfile=cat, format='ascii', catalog_type='srl', clobber=True)

    c = bdsm.process_image(tgss['name'], frequency=147e6, group_by_isl=True, detection_image=image_mask, beam=(1.25000e-02, 1.25000e-02, 0.0), mean_map='zero', rms_map=False, rms_value=0.1)
    cat = tgss['name'].replace('.fits','.cat')
    if c.sources == []:
        open(cat, 'a').close()
    else:
        c.write_catalog(outfile=cat, format='ascii', catalog_type='srl', clobber=True)


