#!/usr/bin/env python3

# for a list of pointings, get the closest and mosaic them

import os, sys, glob
from astropy.table import Table
from astropy.io import fits as pyfits
from astropy.wcs import WCS as pywcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

file_suffix = '-wide-ddc1.fits' # this is used to isolate file names
# each one needs to have a corresponding "-avgbeam.fits"
cutout_size = 3.3 # size of the final cutout [deg]
beam_size = 3.7 # approx beam RADIUS at the black point, usually 30% of the beam power [deg]
grid_file = 'allsky-grid.fits'

class Pointing():

    def __init__(self, pointing_file):
        self.pointing_file = pointing_file
        self.beam_file = pointing_file.replace(file_suffix,'-avgbeam.fits')
        assert os.path.exists(self.beam_file)
        self.output_file = pointing_file.replace(file_suffix,'-mosaic.fits')
        with pyfits.open(pointing_file) as f:
            head = f[0].header
        ra = head['CRVAL1']
        dec = head['CRVAL2']
        self.wcs = pywcs(head)
        self.coord = SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')

    def add_closest(self, dist):
        print("Max ditance: %f deg" % dist)
        grid['dist'] = self.coord.separation(SkyCoord(grid['ra'],grid['dec'],frame='fk5'))
        self.pointing_closest_files = [f.lower()+file_suffix for f in grid[grid['dist'] < dist*u.deg]['name']]
        print('distances:', grid[grid['dist'] < dist*u.deg])
        
        # track also beam files
        self.beam_closest_files = [pointing_closest_file.replace(file_suffix,'-avgbeam.fits') for pointing_closest_file in self.pointing_closest_files]

    def mosaic(self):
        in_pointings = [self.pointing_file]+self.pointing_closest_files
        in_beams = [self.beam_file]+self.beam_closest_files
        print('%s - mosaiching:' % self.output_file, in_pointings)
        os.system('mosaic.py --images %s --beams %s --beamcorr --beamcut 0.3 --find_noise --header %s --output %s' % \
                (' '.join(in_pointings), ' '.join(in_beams), self.output_file, self.output_file))

    def make_empty_mosaic(self):
        rwcs = pywcs(naxis=2)
        rwcs.wcs.ctype = [self.wcs.wcs.ctype[0],self.wcs.wcs.ctype[1]]
        rwcs.wcs.cdelt = [self.wcs.wcs.cdelt[0],self.wcs.wcs.cdelt[1]]
        rwcs.wcs.crval = [self.coord.ra.deg,self.coord.dec.deg]
        rwcs.wcs.crpix = [2000,2000]
        regrid_hdr = rwcs.to_header()
        regrid_hdr['NAXIS'] = 2
        regrid_hdr['NAXIS1'] = 4000
        regrid_hdr['NAXIS2'] = 4000
        data = np.zeros((4000,4000))
        pyfits.writeto(self.output_file, header=regrid_hdr, data=data, overwrite=True)

pointing_files = sorted(glob.glob('*'+file_suffix))

# open grid file and restrict to available pointings
pointing_names = [p.replace(file_suffix,'') for p in pointing_files]
print('Use pointings:', pointing_names)
grid = Table.read(grid_file)
selection = [p['name'].lower() in pointing_names for p in grid]
grid = grid[selection]
print('Restricting to %i pointings.' % len(grid))

# list of pointings to mosaic
for pointing_file in pointing_files:
    pointing = Pointing(pointing_file)
    # get the closest (assuming worst case: perfectly diagonally aligned pointings)
    dist = np.sqrt(2.)/2*cutout_size+beam_size
    pointing.add_closest(dist)
    pointing.make_empty_mosaic()
    pointing.mosaic()

