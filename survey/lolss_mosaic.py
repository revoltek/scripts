#!/usr/bin/env python3

# for a list of pointings, get the closest and mosaic them

import os, sys, glob
from astropy.table import Table

file_suffix = '-wide-ddc1.fits' # this is used to isolate file names
# each one needs to have a corresponding "-avgbeam.fits"
cutout_size = 3.3 # size of the final cutout [deg]
beam_size = # approx beam width at the black point, usually 30% of the beam power [deg]
grid_file = 'allsky-grid.fits'

class Pointing():
    def __init__(self, pointing_file):
        self.pointing_file = pointing_file
        self.beam_file = pointing_file.replace(file_suffix,'-avgbeam.fits')
        assert os.system.exists(self.beam_file)
        self.output_file = pointing_file.replace(file_suffix,'-mosaic.fits')
        with fits.open(pointing_file) as f:
            head = f[0].header
            ra = head['CRVAL1']
            dec = head['CRVAL2']
        self.coord = SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')

    def add_closest(self, dist):
        coord = SkyCoord(grid['ra']*u.deg, grid['dec']*u.deg, frame='fk5')
        grid['dist'] = coord.separation(SkyCoord(ras,decs,frame='fk5'))
        self.pointings_closest_files = [f+file_suffix for f in grid[grid['dist'] < dist*u.deg]['name']]
        # track also beam files
        self.beam_closest_files = [pointings_closest_file.replace(file_suffix,'-avgbeam.fits') for pointings_closest_file in pointings_closest_files]

    def mosaic(self):
        in_pointings = [self.pointing_file]+self.pointing_closest_files
        in_beams = [self.beam_file]+self.beam_closest_files
        print('%s - mosaiching:' % self.output_file, in_pointings)
        os.system('mosaic.py --images %s --beams %s --beamcorr --find_noise --output %s' %(' '.join(in_pointings), ' '.join(in_beams), self.output_file))

    def cutout(self):
        pass

pointing_files = glob.glob('*'+file_suffix)

# open grid file and restrict to available pointings
pointing_names = [p.replce(file_suffix,'') for p in pointing_files]
print('Use pointings:', pointing_names)
grid = Table.read(grid_file)
selection [p['name'] in pointing_names for p in grid]
grid = grid[selection]
print('Restricting to %i pointings.' % len(grid))

# list of pointings to mosaic
for pointings_file in pointing_files:
    pointing = Pointing(pointing_file)
    # get the closest (assuming worst case: perfectly diagonally aligned pointings)
    dist = np.sqrt(2.)/2*cutout_size*beam_size
    pointing.add_closest(dist)
    pointing.mosaic()
    pointing.cutout()
