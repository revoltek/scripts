#!/usr/bin/env python3

import os, sys, glob

hdr_final = 'mosaic-rms.hdr'
tbl_proj  = 'mosaic-rms.tbl'
dir_raw = 'mosaic-rms'
dir_proj = 'mosaic-rms/projected'
fits_output = 'mosaic-rms.fits'

if not os.path.exists(dir_raw): 
    print('Missing directory with raw files: %s' % dir_raw)
    sys.exit()

if not os.path.exists(dir_proj): os.system('mkdir %s' % dir_proj)

# project all files on the same pixel space
fits_raws = glob.glob(dir_raw+'/*fits')
for fits_raw in fits_raws:
    fits_proj = dir_proj+'/'+os.path.basename(fits_raw).replace('.fits','-proj.fits')
    print (fits_proj)
    if not os.path.exists(fits_proj):
        os.system('mProjectPP -d 1 %s %s %s' % (fits_raw, fits_proj, hdr_final))

#os.system('rm %s/*area*' % dir_proj)
# make info table for mAdd
os.system('mImgtbl %s %s' % (dir_proj, tbl_proj))
# add images
os.system('mAdd -d 1 -p %s %s %s %s' % (dir_proj, tbl_proj, hdr_final, fits_output))
