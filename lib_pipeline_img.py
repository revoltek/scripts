#!/usr/bin/python

import os, sys
import logging
import numpy as np


#def size_from_facet(img, c_coord, pixsize):
#    """
#    Given an image, a new centre find the smallest image size which cover the whole image.
#    img = CASA-image name
#    c_coord = [ra,dec] in degrees, the wanted image centre
#    pixsize = in arcsec, the final image will have this pixels size, so a rescaling might be needed
#    """
#    import pyrap.images
#    img = pyrap.images.image(img)
#    c = img.coordinates()
#    # assumes images in a standard casa shape
#    assert c.get_axes()[2] == ['Declination', 'Right Ascension']
#    # assume same increment in image axes
#    assert abs(c.get_increment()[2][0]) == abs(c.get_increment()[2][1])
#    cen_y, cen_x = img.topixel([1,1,c_coord[1]*np.pi/180., c_coord[0]*np.pi/180.])[2:]
#    max_y, max_x = img.shape()[2:]
#    max_dist = max(max(cen_x, max_x - cen_x), max(cen_y, max_y - cen_y))
#    max_dist = max_dist * abs(c.get_increment()[2][0])*180/np.pi*3600 / pixsize
#    if max_dist > 6400: return 6400
#    # multiply distance *2 (so to have the image size) and add 100% to be conservative
#    max_dist = (max_dist*2)*2
#    goodvalues = np.array([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512, 256, 128])
#    shape = min(goodvalues[np.where(goodvalues>=max_dist)])
#    del img
#    return shape

def flatten(f, channel=0, freqaxis=0):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """
    from astropy import wcs

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        return f[0].header,f[0].data

    w = wcs.WCS(f[0].header)
    wn=wcs.WCS(naxis=2)

    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]

    header = wn.to_header()
    header["NAXIS"]=2
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    slice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            slice.append(np.s_[:],)
        elif i==freqaxis:
            slice.append(channel)
        else:
            slice.append(0)

    # slice=(0,)*(naxis-2)+(np.s_[:],)*2
    return header,f[0].data[slice]


def size_from_reg(filename, region, coord, pixscale, pad=1.2):
    """
    find the minimum image size to cover a certain region
    pad: multiplicative factor on the final size
    """
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    import pyregion

    # open fits
    fits = pyfits.open(filename)
    header, data = flatten(fits)

    # extract mask
    r = pyregion.open(region)
    mask = r.get_mask(header=header, shape=data.shape)
    y,x = mask.nonzero()
    min_ra_pix = np.min(x)
    min_dec_pix = np.min(y)
    max_ra_pix = np.max(x)
    max_dec_pix = np.max(y)
#    print 'min ra - min dec (pix)', min_ra_pix, min_dec_pix
#    print 'max ra - max dec (pix)', max_ra_pix, max_dec_pix

    # to degrees
    w = pywcs.WCS(fits[0].header)
    min_ra, min_dec = w.all_pix2world(min_ra_pix, min_dec_pix, 0, 0, 0, ra_dec_order=True)
    max_ra, max_dec = w.all_pix2world(max_ra_pix, max_dec_pix, 0, 0, 0, ra_dec_order=True)
#    print 'min ra - min dec (sky)', min_ra, min_dec
#    print 'max ra - max dec (sky)', max_ra, max_dec

    max_ra_size = 2*max([abs(coord[0] - min_ra), abs(coord[0] - max_ra)])*3600 # arcsec
    max_dec_size = 2*max([abs(coord[1] - min_dec), abs(coord[1] - max_dec)])*3600 # arcsec

    return int(pad*max([max_ra_size, max_dec_size])/pixscale)

 
def scale_from_ms(ms):
    """
    Get the pixel scale in arcsec for a full-res image.
    It is 1/3 of the max resolution assuming zenit observation.
    """
    from pyrap.tables import *
    import numpy as np

    c = 299792458

    t = table(ms, ack=False)
    col = t.getcol('UVW')
    maxdist = 0
    #mindist=np.inf

    t = table(ms+'/SPECTRAL_WINDOW', ack=False)
    wavelenght = c/t.getcol('REF_FREQUENCY')[0]
    #print 'Wavelenght:', wavelenght,'m (Freq: '+str(t.getcol('REF_FREQUENCY')[0]/1.e6)+' MHz)'

    for u,v,w in col:
        dist = np.sqrt(u*u+v*v)
        if dist > maxdist: maxdist = dist
    #    if dist < mindist and dist != 0.0: mindist = dist

    #print 'Min scale: ~',wavelenght/maxdist*(180/np.pi)*3600, 'arcsec'
    return int(round(wavelenght/maxdist*(180/np.pi)*3600/3.)) # arcsec

def blank_image_fits(filename, maskname, outfile=None, inverse=False, blankval=0.):
    """
    Set to "blankval" all the pixels inside the given region
    if inverse=True, set to "blankval" pixels outside region.

    filename: fits file
    region: ds9 region
    outfile: output name
    inverse: reverse region mask
    blankval: pixel value to set
    """
    import astropy.io.fits as pyfits
    import pyregion

    if outfile == None: outfile = filename

    with pyfits.open(maskname) as fits:
        mask = fits[0].data
    
    if inverse: mask = ~(mask.astype(bool))

    with pyfits.open(filename) as fits:
        data = fits[0].data

        assert mask.shape == data.shape # mask and data should be same shape

        sum_before = np.sum(data)
        data[mask] = blankval
        print "Sum of values: %f -> %f" % (sum_before, np.sum(data))
        fits.writeto(outfile, clobber=True)

 
def blank_image_reg(filename, region, outfile=None, inverse=False, blankval=0.):
    """
    Set to "blankval" all the pixels inside the given region
    if inverse=True, set to "blankval" pixels outside region.

    filename: fits file
    region: ds9 region
    outfile: output name
    inverse: reverse region mask
    blankval: pixel value to set
    """
    import astropy.io.fits as pyfits
    import pyregion

    if outfile == None: outfile = filename

    # open fits
    with pyfits.open(filename) as fits:
        origshape = fits[0].data.shape
        header, data = flatten(fits)
        # extract mask
        r = pyregion.open(region)
        mask = r.get_mask(header=header, shape=data.shape)
        if inverse: mask = ~mask
        data[mask] = blankval
        # save fits
        fits[0].data = data.reshape(origshape)
        fits.writeto(outfile, clobber=True)


def get_coord_centroid(filename, region):
    """
    Get centroid coordinates from an image and a region
    filename: fits file
    region: ds9 region
    """
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    import pyregion
    from scipy.ndimage.measurements import center_of_mass

    fits = pyfits.open(filename)
    header, data = flatten(fits)

    # extract mask and find center of mass
    r = pyregion.open(region)
    mask = r.get_mask(header=header, shape=data.shape)
    dec_pix, ra_pix = center_of_mass(mask)
    
    # convert to ra/dec in angle
    w = pywcs.WCS(fits[0].header)
    #w = w.celestial # needs newer astropy
    ra, dec = w.all_pix2world(ra_pix, dec_pix, 0, 0, 0, ra_dec_order=True)

    fits.close()
    return float(ra), float(dec)




