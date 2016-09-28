#!/usr/bin/python

import os, sys
import logging
import numpy as np


def size_from_facet(img, c_coord, pixsize):
    """
    Given an image, a new centre find the smallest image size which cover the whole image.
    img = CASA-image name
    c_coord = [ra,dec] in degrees, the wanted image centre
    pixsize = in arcsec, the final image will have this pixels size, so a rescaling might be needed
    """
    import pyrap.images
    img = pyrap.images.image(img)
    c = img.coordinates()
    # assumes images in a standard casa shape
    assert c.get_axes()[2] == ['Declination', 'Right Ascension']
    # assume same increment in image axes
    assert abs(c.get_increment()[2][0]) == abs(c.get_increment()[2][1])
    cen_y, cen_x = img.topixel([1,1,c_coord[1]*np.pi/180., c_coord[0]*np.pi/180.])[2:]
    max_y, max_x = img.shape()[2:]
    max_dist = max(max(cen_x, max_x - cen_x), max(cen_y, max_y - cen_y))
    max_dist = max_dist * abs(c.get_increment()[2][0])*180/np.pi*3600 / pixsize
    if max_dist > 6400: return 6400
    # multiply distance *2 (so to have the image size) and add 100% to be conservative
    max_dist = (max_dist*2)*2
    goodvalues = np.array([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512, 256, 128])
    shape = min(goodvalues[np.where(goodvalues>=max_dist)])
    del img
    return shape

def blank_image(filename, region, inverse=False, blankval=0.):
    """
    Set to "blankval" all the pixels inside the given region
    if inverse=True, set to "blankval" pixels outside region.
    """
    import pyfits
    import pyregion

    # open fits
    fits = pyfits.open(filename)
    hdu = fits[0]
    data = fits[0].data
    # extract mask
    mask = region.get_mask(hdu, shape=np.shape(data))
    # set pixel to blankval
    data[mask] = blankval
    # save fits
    fits[0].data = data
    fits.writeto(filename)
