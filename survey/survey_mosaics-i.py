#!/usr/bin/env python3

# for a list of pointings, get the closest and mosaic them

import os, sys, glob
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS as pywcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from LiLF.surveys_db import SurveysDB

dir_mosaics = "/homes/fdg/storage/surveytgts/mosaics/"
dir_done = "/homes/fdg/storage/surveytgts/done/"

file_MFS_i = 'wideDD-c00-MFS-image.fits' # this is used to isolate file names
nchans = None
#file_suffix = '-wide-cube.fits' # this is used to isolate file names
#nchans = 6

file_beam = 'primarybeam.fits' # this is used to isolate file names
beam_size = 4 # approx beam RADIUS at the blank point, usually 30% of the beam power [deg] - measured at low dec = 24deg
grid_file = '/homes/fdg/storage/allsky-grid.fits'
file_header_template = '/homes/fdg/storage/scripts/survey/LoLSS-bkp/headers_template.hdr'

min_good = 2 # minimum number of good pointings to mosaic

class Pointing():

    def __init__(self, pointing_name, pointings_name_nearby):
        self.pointing_name = pointing_name
        self.pointing_files = [dir_done+p+'/'+file_MFS_i for p in pointings_name_nearby]
        self.beam_files = [dir_done+p+'/'+file_beam for p in pointings_name_nearby]
        self.output_file = dir_mosaics+pointing_name+'-mosaicI.fits'

    def mosaic(self):
        print('mosaic.py --images %s --beams %s --beamcorr --beamcut 0.3 --find_noise --header %s --output %s > %s.log 2>&1' % \
                (' '.join(self.pointing_files), ' '.join(self.beam_files), self.output_file, self.output_file, self.output_file))
        os.system('mosaic.py --images %s --beams %s --beamcorr --beamcut 0.3 --find_noise --header %s --output %s > %s.log 2>&1' % \
                (' '.join(self.pointing_files), ' '.join(self.beam_files), self.output_file, self.output_file, self.output_file))

    def make_empty_mosaic(self, file_header_template, ra, dec):
        print('Creating empty mosaic for %s' % self.output_file)
        with open(file_header_template) as f:
            hdr_str = f.read()
        self.head = fits.Header.fromstring(hdr_str, sep='\n')
        self.head['CRVAL1'] = ra 
        self.head['CRVAL2'] = dec
        self.head['OBJECT'] = self.pointing_name
        data = np.zeros((3500,3500))
        fits.writeto(self.output_file, header=self.head, data=data, overwrite=True)

# calculate matrix of distances
grid = Table.read(grid_file)
grid = grid[grid['dec'] > 24]
coords = SkyCoord(ra=grid['ra'], dec=grid['dec'], frame='icrs')
# Do a “search around” with itself, using a 4° threshold
idx1, idx2, sep2d, _ = coords.search_around_sky(coords, 4 * u.deg)
unique = idx1 <= idx2
idx1 = idx1[unique]
idx2 = idx2[unique]
# Build a mapping: object_id → list of neighbor_ids within 4°
neighbors = {obj_id: [] for obj_id in grid['name']}
for ii, jj in zip(idx1, idx2):
    id_i = grid['name'][ii]
    id_j = grid['name'][jj]
    neighbors[id_i].append(id_j)
    neighbors[id_j].append(id_i)  # since the relation is symmetric

with SurveysDB(survey='lba',readonly=True) as sdb:
    sdb.execute('SELECT id,status FROM fields')
r = sdb.cur.fetchall()
pointings_status = {r[i]['id']: r[i]['status'] for i in range(len(r))}

#pointing_files = sorted(glob.glob(dir_mosaics+file_MFS_i))
# list of pointings to mosaic
for pointing, close_pointings in neighbors.items():
    close_pointings = np.unique(close_pointings) # remove double entry of self-pointing
    #print('Pointing:', pointing, close_pointings)
    # skip not observed pointings
    if pointings_status[pointing+'o'] == 'Not observed' and \
       pointings_status[pointing+'s'] == 'Not observed': continue

    # count how many good pointings we have
    good = 0; bad = 0; good_pointings = []
    for close_pointing in close_pointings:
        if pointings_status[close_pointing+'o'] == 'Done': 
            good += 1
            good_pointings.append(close_pointing+'o')
        elif pointings_status[close_pointing+'o'] == 'Not observed': pass
        else: bad += 1
        if pointings_status[close_pointing+'s'] == 'Done': 
            good += 1
            good_pointings.append(close_pointing+'s')
        elif pointings_status[close_pointing+'s'] == 'Not observed': pass
        else: bad += 1
    
    if good >= min_good:
        print(f'{pointing}: Enough good pointing done ({good}/{good+bad} - {good_pointings}), proceeding.')
    elif bad == 0:
        print(f'{pointing}: All pointings are done ({good}/{good+bad}), proceeding')
    else:
        #print(f'{pointing}: Not all pointings are done ({good}/{good+bad} - {good_pointings}), skipping.')
        continue    

    p = Pointing(pointing, good_pointings)
    if not os.path.exists(p.output_file):
        ra = grid['ra'][grid['name'] == pointing][0]
        dec = grid['dec'][grid['name'] == pointing][0]
        p.make_empty_mosaic(file_header_template, ra, dec)
        p.mosaic()