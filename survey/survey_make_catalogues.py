#!/usr/bin/python3
# make_catalogues.py *mosaic.fits [run into mosaics dir]

import bdsf
import os, sys, glob

infiles = sys.argv[1:]
restfrq = 54e6

if not os.path.exists('catalogues'): os.system('mkdir catalogues')
print("All files:", infiles)

for infile in infiles:
    print('Working on %s' % infile)
    img = bdsf.process_image(infile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=50, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=False, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None, frequency=restfrq)

    catprefix = 'catalogues/'+infile.replace('.fits','')
    img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True',clobber=True)
    img.write_catalog(outfile=catprefix +'.gaus.fits',catalog_type='gaul',format='fits',correct_proj='True',clobber=True)
    img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)
    img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
    img.export_image(outfile=catprefix +'.pybdsmmask.fits',img_type='island_mask',img_format='fits',clobber=True)
    img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True',clobber=True)
