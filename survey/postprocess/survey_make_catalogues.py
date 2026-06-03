#!/usr/bin/python3
# make_catalogues.py *mosaic.fits [run into mosaics dir]

import bdsf
import os, glob, logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

restfrq = 54e6
dir_mosaics = "/homes/fdg/storage/surveytgts/mosaics/healpix/"
dir_catalogues = "/homes/fdg/storage/surveytgts/catalogues/"
infiles = sorted(glob.glob(dir_mosaics+'*.fits'))

if not os.path.exists(dir_catalogues): os.makedirs(dir_catalogues)
logger.info('All files: %s', infiles)

for infile in infiles:
    catprefix = dir_catalogues+os.path.basename(infile).replace('.fits','')

    if os.path.exists(catprefix+'.cat.fits'):
        logger.info('Skipping %s (catalogue already exists)', infile)
        continue

    logger.info('Working on %s', infile)
    img = bdsf.process_image(infile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=50, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, atrous_do=False, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None, frequency=restfrq)
    img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True',clobber=True)
    img.write_catalog(outfile=catprefix +'.gaus.fits',catalog_type='gaul',format='fits',correct_proj='True',clobber=True)
    img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)
    img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
    img.export_image(outfile=catprefix +'.pybdsmmask.fits',img_type='island_mask',img_format='fits',clobber=True)
    img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True',clobber=True)
