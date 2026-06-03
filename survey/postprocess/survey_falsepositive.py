#!/usr/bin/env python3

# to run on the whoel survey, in the mosaic-res dir use
"""
import os, sys, glob

mosaics = sorted(glob.glob('*mosaic.fits'))

for mosaic in mosaics:
    print('Working on:', mosaic)
    basename = mosaic.replace('.fits', '')
    rmsmap = glob.glob('../mosaic-i/%s_pybdsm/*/background/%s.pybdsm.rmsd_I.fits' % (basename, basename))[0]
    meanmap = glob.glob('../mosaic-i/%s_pybdsm/*/background/%s.pybdsm.mean_I.fits' % (basename, basename))[0]
    os.system('~/scripts/survey/survey_falsepositive.py %s --rmsmap %s --meanmap %s' % (mosaic, rmsmap, meanmap))
"""

import bdsf
import argparse
import os, sys, glob, re
from astropy.io import fits

parser = argparse.ArgumentParser(description='Completness testing')
parser.add_argument('file', metavar='FILE', nargs='+',
                   help='Root name of file to process (a residual image)')
parser.add_argument('-r','--rmsmap',dest='rmsmap',action='store',
                    default=None,
                    help='BDSF RMS map (important to re-use the originally created one as the "empty" image has different rms properties)')
parser.add_argument('-e','--meanmap',dest='meanmap',action='store',
                    default=None,
                    help='BDSF mean map')
args = parser.parse_args()

infile = args.file[0]
restfrq = 54e6

if not os.path.exists('catalogues-inv'): os.system('mkdir catalogues-inv')

print('Working on %s' % infile)
fakefile = re.sub('.fits','.fake.fits', infile)
fp = fits.open(infile)
f = fp[0].data
f *= -1
fp.writeto(fakefile, overwrite=True)
fp.close()

img=bdsf.process_image(fakefile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=50, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0,output_opts=True, output_all=True, atrous_do=False,atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None,frequency=restfrq, rmsmean_map_filename=[args.meanmap,args.rmsmap])

catprefix = 'catalogues-inv/'+infile.replace('.fits','')
img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True',clobber=True)
img.write_catalog(outfile=catprefix +'.gaus.fits',catalog_type='gaul',format='fits',correct_proj='True',clobber=True)
