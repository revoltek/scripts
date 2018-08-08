#!/usr/bin/python

import bdsf
import sys

f = sys.argv[1]

img = bdsf.process_image(f, advanced_opts=True, detection_image=f.replace('-pb',''), interactive=True, thresh_pix=5., thresh_isl=3., \
        adaptive_rms_box=True, rms_box_bright=(100,30), adaptive_thresh=10.)
img.write_catalog(format='fits',catalog_type='srl', clobber=True)
img.write_catalog(format='ds9',catalog_type='srl', clobber=True)
img.write_catalog(format='ds9',catalog_type='gaul', clobber=True)
img.export_image(outfile=f+'_gaus_resid.fits', img_type='gaus_resid',clobber=True)
img.export_image(outfile=f+'_gaus_model.fits', img_type='gaus_model',clobber=True)
img.export_image(outfile=f+'_rms.fits', img_type='rms',clobber=True)
img.export_image(outfile=f+'_mean.fits', img_type='mean',clobber=True)
