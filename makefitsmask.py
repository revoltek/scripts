#!/usr/bin/python
# Take a fits template and a ds9 region and create a mask
# Usage: makefitsmask.py template.fits ds9.reg
# output: output.fits

import os, sys
import radioflux2 as radioflux
import pyregion
import numpy as np
from astropy.io import fits

rm = radioflux.radiomap(sys.argv[1])
region = pyregion.open(sys.argv[2]).as_imagecoord(rm.headers[0])
mask = region.get_mask(hdu=rm.f,shape=np.shape(rm.d[0]))
print("Masking %.2f %% of the map." % (100.*np.sum(mask)/np.size(mask)))
hdulist = fits.open(sys.argv[1])
hdulist[0].data = mask.astype(np.uint8)
hdulist.writeto('output.fits')
