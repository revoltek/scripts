#!/usr/bin/python

import lofar.bdsm as bdsm
import sys

file = sys.argv[1]

img = bdsm.process_image(file, advanced_opts=True, blank_zeros=True, interactive=True, thresh_pix=6., thresh_isl=4.5)
img.write_catalog(format='fits',catalog_type='srl', clobber=True) # for bdsm2cat.py
img.write_catalog(format='bbs', bbs_patches=None,catalog_type='gaul', clobber=True) # for sagecal to trnsform
img.write_catalog(format='ds9',catalog_type='gaul', clobber=True) # for ds9
img.export_image(outfile=file+'_gaus_resid.fits', img_type='gaus_resid',clobber=True)
img.export_image(outfile=file+'_gaus_model.fits', img_type='gaus_model',clobber=True)
img.export_image(outfile=file+'_rms.fits', img_type='rms',clobber=True)
img.export_image(outfile=file+'_mean.fits', img_type='mean',clobber=True)

#~/scripts/bdsm2cat.py --snr 1 --sr 30 --phase_center "187.705833, 12.391111" awimager-145.pybdsm.srl.fits
