#!/usr/bin/env python3

import glob
import subprocess
import argparse
from astropy.table import Table
from tqdm import tqdm
import os
from astropy_healpix import HEALPix
import astropy.units as u
from astropy.coordinates import SkyCoord

def run(s):
    print(s)
    subprocess.run(s, shell=True, check=True)

parser = argparse.ArgumentParser(description='Concatenate survey catalogues')
parser.add_argument('--catdir', default=os.path.expanduser('~/storage/surveytgts/catalogues/'),
                    help='Directory containing the per-HEALPix catalogue FITS files (default: ~/storage/surveytgts/catalogues/)')
parser.add_argument('--force', action='store_true', default=False,
                    help='Re-process files even if filtered output already exists')
parser.add_argument('--output', default='LoLSS_DR2',
                    help='Base name for the output FITS files (default: LoLSS_DR2)')
args = parser.parse_args()

CATALOGUE_DIR = args.catdir
force = args.force
    
hp = HEALPix(nside=16)
print(hp.npix,'total healpix pixels on sky')
area=hp.pixel_area.to(u.deg**2).value
print('area of one healpix is',area,'sq. deg')

os.chdir(CATALOGUE_DIR)

for cattype in ['cat','gaus']:
    with open(f'{cattype}-catlist.txt','w') as flist:
        for infile in tqdm(glob.glob(f'HP*-mosaicI.{cattype}.fits')):
            pix = infile.split('-')[0][2:]
            outname=infile.replace('.fits','-filtered.fits')
            if not os.path.isfile(outname) or force:
                t = Table.read(infile)
                # select only sources in this healpix pixel (mosaics cover multiple pixels)
                t['HEALPIX'] = hp.lonlat_to_healpix(t['RA'].data*u.deg,t['DEC'].data*u.deg)
                t = t[t['HEALPIX']==int(pix)]
                if len(t)==0:
                    print(f'*** empty table {infile}, skipping ***')
                    continue
                for k in t.colnames:
                    if 'img_plane' in k or k.startswith('Resid') or k.endswith('_max') or k in ['Flag_beam','Isl_mean','Source_id','Isl_id']:
                        del t[k]
                for k in ['E_RA','E_DEC','Maj','E_Maj','Min','E_Min','DC_Maj','E_DC_Maj','DC_Min','E_DC_Min']:
                    t[k].convert_unit_to(u.arcsec)
                for k in ['Peak_flux','E_Peak_flux','Isl_rms']:
                    t[k].convert_unit_to(u.mJy/u.beam)
                for k in ['Total_flux','E_Total_flux','Isl_Total_flux','E_Isl_Total_flux']:
                    t[k].convert_unit_to(u.mJy)
                sc = SkyCoord(t['RA'],t['DEC'],frame='icrs')
                strings = sc.to_string(style='hmsdms',sep='',precision=2)
                ilt = [('LoLJ'+s).replace(' ','')[:-1] for s in strings]
                t.add_column(ilt, name='Source_Name', index=0)
                t.write(outname, overwrite=True)
            flist.write(outname+'\n')
    print('Now running STILTS')
    run(f'stilts tcat in=@{cattype}-catlist.txt out={args.output}-{cattype}.fits lazy=true ocmd="tablename {args.output}_{cattype}; sort RA"')