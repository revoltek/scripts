#!/usr/bin/env python3

import os, sys, argparse, logging
import numpy as np
from astropy.io import fits as pyfits
from lib_linearfit import linear_fit, linear_fit_bootstrap
from lib_fits import AllImages
logging.root.setLevel(logging.DEBUG)

spidx_low = 'spidx1.fits'
spidx_high = 'spidx2.fits'

def get_args():
    parser = argparse.ArgumentParser(description='Make curvature maps, e.g. curvaturemap.py  *fits')
    parser.add_argument('images', nargs='+', help='List of images to use for curvature')
    parser.add_argument('--region', dest='region', type=str, help='Ds9 region to restrict analysis')
    parser.add_argument('--output', dest='output', default='curvature.fits', type=str, help='Name of output mosaic (default: curvature.fits)')

    args = parser.parse_args()

    # check input
    if len(args.images) != 2:
        logging.error('Requires exactly 2 images.')
        sys.exit(1)

    return args


if __name__ == '__main__':
    args = get_args()
    
    all_images = AllImages(args.images)
    regrid_hdr = all_images.regrid_common(region=args.region, action='regrid_header')
    xsize = regrid_hdr['NAXIS1']
    ysize = regrid_hdr['NAXIS2']

    data_low = all_images[0].img_data
    data_high = all_images[1].img_data

    bamj, bmin, bpa = all_images[0].get_beam()
    regrid_hdr['BMAJ'] = bamj; regrid_hdr['BMIN'] = bmin; regrid_hdr['BPA'] = bpa
    regrid_hdr.add_history(' '.join(sys.argv))

    pyfits.writeto(args.output, data_low-data_high, regrid_hdr, overwrite=True)

