#!/usr/bin/python

import bdsf
import sys

file = sys.argv[1]

img = bdsf.process_image(file, advanced_opts=True, blank_zeros=True, interactive=True, thresh_pix=6., thresh_isl=4.5)
img.write_catalog(format='fits',catalog_type='srl', clobber=True)
img.write_catalog(format='ds9',catalog_type='gaul', clobber=True)
img.export_image(outfile=file+'_gaus_resid.fits', img_type='gaus_resid',clobber=True)
img.export_image(outfile=file+'_gaus_model.fits', img_type='gaus_model',clobber=True)
img.export_image(outfile=file+'_rms.fits', img_type='rms',clobber=True)
img.export_image(outfile=file+'_mean.fits', img_type='mean',clobber=True)
