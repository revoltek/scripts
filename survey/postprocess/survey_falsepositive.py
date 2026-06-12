#!/usr/bin/env python3
# Run on the whole survey, e.g. in the mosaic-res dir:
#   survey_falsepositive.py *mosaic.fits

import bdsf
import argparse
import os, sys, re
from astropy.io import fits

parser = argparse.ArgumentParser(description='Completness testing')
parser.add_argument('-m','--mosaicdir',dest='mosaicdir',action='store',
                    default=None, required=True,
                    help='Directory containing the mosaic FITS files')
parser.add_argument('-c','--catdir',dest='catdir',action='store',
                    default=None, required=True,
                    help='Directory containing BDSF RMS and mean maps, named as <mosaic>.pybdsm.rms.fits / <mosaic>.pybdsm.mean.fits')
args = parser.parse_args()

restfrq = 54e6
os.makedirs(os.path.join(args.catdir, 'catalogues-inv'), exist_ok=True)
os.makedirs(args.mosaicdir.rstrip('/') + '-inv', exist_ok=True)

for infile in sorted(os.listdir(args.mosaicdir)):
    if not infile.endswith('.fits'):
        continue
    infile = os.path.join(args.mosaicdir, infile)
    print('Working on %s' % infile)
    basename = os.path.basename(infile).replace('.fits', '')

    inv_dir = args.mosaicdir.rstrip('/') + '-inv'
    rmsmap  = os.path.relpath(os.path.join(args.catdir, basename + '.rms.fits'), inv_dir)
    meanmap = os.path.relpath(os.path.join(args.catdir, basename + '.mean.fits'), inv_dir)

    if not os.path.exists(os.path.join(inv_dir, rmsmap)):
        print('ERROR: RMS map not found: %s, skipping.' % rmsmap)
        continue
    if not os.path.exists(os.path.join(inv_dir, meanmap)):
        print('ERROR: mean map not found: %s, skipping.' % meanmap)
        continue

    fakefile = os.path.join(args.mosaicdir.rstrip('/') + '-inv', os.path.basename(infile).replace('.fits', '.inverted.fits'))
    fp = fits.open(infile)
    f = fp[0].data
    f *= -1
    fp.writeto(fakefile, overwrite=True)
    fp.close()

    img = bdsf.process_image(fakefile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=50, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, atrous_do=False, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq, rmsmean_map_filename=[meanmap, rmsmap])

    catprefix = os.path.join(args.catdir, 'catalogues-inv', os.path.basename(infile).replace('.fits', ''))
    img.write_catalog(outfile=catprefix + '.cat.fits', catalog_type='srl', format='fits', correct_proj='True', clobber=True)
    img.write_catalog(outfile=catprefix + '.gaus.fits', catalog_type='gaul', format='fits', correct_proj='True', clobber=True)
