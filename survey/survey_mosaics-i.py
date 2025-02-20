#!/usr/bin/env python3

# for a list of pointings, get the closest and mosaic them

import os, sys, glob
from astropy.table import Table
from astropy.io import fits as pyfits
from astropy.wcs import WCS as pywcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

file_suffix = '-wide-v.fits' # this is used to isolate file names
nchans = None
#file_suffix = '-wide-cube.fits' # this is used to isolate file names
#nchans = 6
cutout_size = 3.3 # size of the final cutout [deg]
beam_size = 3.7 # approx beam RADIUS at the blank point, usually 30% of the beam power [deg]
grid_file = '../allsky-grid.fits'
beamdir = '../beams/'
stokes = 1

class Pointing():

    def __init__(self, pointing_file, channel=None, stokes=0):
        self.pointing_file = pointing_file
        if channel is not None:
            self.beam_file = beamdir+pointing_file.replace(file_suffix,'-beam%04i.fits' % channel)
            self.output_file = pointing_file.replace(file_suffix,'-mosaic-chan%02i.fits' % channel)
        else:
            self.beam_file = beamdir+pointing_file.replace(file_suffix,'-avgbeam.fits')
            self.output_file = pointing_file.replace(file_suffix,'-mosaic.fits')
        if not os.path.exists(self.beam_file):
            sys.exit('Missing %s' % self.beam_file)

        with pyfits.open(pointing_file) as f:
            self.head = f[0].header
        if channel:
            assert self.head['CTYPE4']  == 'FREQ'
            self.head['RESTFREQ'] = self.head['CRVAL4']+self.head['CDELT4']*channel
        ra = self.head['CRVAL1']
        dec = self.head['CRVAL2']
        self.wcs = pywcs(self.head)
        self.coord = SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')
        self.channel = channel
        self.stokes = stokes

    def add_closest(self, dist):
        print("Max ditance: %f deg" % dist)
        grid['dist'] = self.coord.separation(SkyCoord(grid['ra'],grid['dec'],frame='fk5'))
        self.pointing_closest_files = [f.lower()+file_suffix for f in grid[grid['dist'] < dist*u.deg]['name']]
        print('distances:', grid[grid['dist'] < dist*u.deg])
        
        # track also beam files
        if self.channel is not None:
            self.beam_closest_files = [beamdir+pointing_closest_file.replace(file_suffix,'-beam%04i.fits' % self.channel) for pointing_closest_file in self.pointing_closest_files]
        else:
            self.beam_closest_files = [beamdir+pointing_closest_file.replace(file_suffix,'-avgbeam.fits') for pointing_closest_file in self.pointing_closest_files]

    def mosaic(self):
        in_pointings = [self.pointing_file]+self.pointing_closest_files
        in_beams = [self.beam_file]+self.beam_closest_files
        if self.channel is not None:
            print('mosaic.py --shift --images %s --beams %s --beamcorr --beamcut 0.3 --find_noise --header %s --use_channel %i --use_stokes %i --output %s > %s.log 2>&1' % \
                (' '.join(in_pointings), ' '.join(in_beams), self.output_file, self.channel, self.stokes, self.output_file, self.output_file))
            os.system('mosaic.py --shift --images %s --beams %s --beamcorr --beamcut 0.3 --find_noise --header %s --use_channel %i --use_stokes %i --output %s > %s.log 2>&1' % \
                (' '.join(in_pointings), ' '.join(in_beams), self.output_file, self.channel, self.stokes, self.output_file, self.output_file))
        else:
            print('mosaic.py --shift --images %s --beams %s --beamcorr --beamcut 0.3 --find_noise --header %s --use_stokes %i --output %s > %s.log 2>&1' % \
                (' '.join(in_pointings), ' '.join(in_beams), self.output_file, self.stokes, self.output_file, self.output_file))
            os.system('mosaic.py --shift --images %s --beams %s --beamcorr --beamcut 0.3 --find_noise --header %s --use_stokes %i --output %s > %s.log 2>&1' % \
                (' '.join(in_pointings), ' '.join(in_beams), self.output_file, self.stokes, self.output_file, self.output_file))

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
        # inherit some headers
        for h in ['RESTFRQ','BTYPE','BUNIT','BMAJ','BMIN','BPA']:
            regrid_hdr[h] = self.head[h]
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
if nchans:
    for chan in range(nchans):
        os.system('rm *__beam.fits') # be sure we do not re-use the old beam
        for pointing_file in pointing_files:
            pointing = Pointing(pointing_file, channel=chan, stokes=stokes)
            if not os.path.exists(pointing.output_file):
                # get the closest (assuming worst case: perfectly diagonally aligned pointings)
                dist = np.sqrt(2.)/2*cutout_size+beam_size
                pointing.add_closest(dist)
                pointing.make_empty_mosaic()
                pointing.mosaic()
else:
    for pointing_file in pointing_files:
        pointing = Pointing(pointing_file, channel=None, stokes=stokes)
        if not os.path.exists(pointing.output_file):
            # get the closest (assuming worst case: perfectly diagonally aligned pointings)
            dist = np.sqrt(2.)/2*cutout_size+beam_size
            pointing.add_closest(dist)
            pointing.make_empty_mosaic()
            pointing.mosaic()
