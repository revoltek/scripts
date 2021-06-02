#!/usr/bin/env python3

# for a list of pointings, get the closest and mosaic them

import os, sys, glob

filesuffix = '-wide-ddc1.fits' # this is used to isolate file names
# each one needs to have a corresponding "-avgbeam.fits"
dist = 4 # dist in deg to search for closest pointings

class Pointing():
    def __init__(self, pointing_file):
        self.pointing_file = pointing_file
        self.beam_file = pointing_file.replace(filesuffix,'-avgbeam.fits')
        assert os.system.exists(self.beam_file)
        self.output_file = pointing_file.replace(filesuffix,'-mosaic.fits')

    def add_closest(self, dist):
        skycov = 
        # track also beam files
        self.beam_closest_files = [pointings_closest_file.replace(filesuffix,'-avgbeam.fits') for pointings_closest_file in pointings_closest_files]

    def mosaic(self):
        in_pointings = [self.pointing_file]+self.pointing_closest_files
        in_beams = [self.beam_file]+self.beam_closest_files
        os.system('mosaic.py --images %s --beams %s --beamcorr --find_noise --output %s' %(' '.join(in_pointings), ' '.join(in_beams), self.output_file))

    def cutout(self):
        pass


# list of pointings to mosaic
for pointings_file in glob.glob('*'+filesuffix):
    pointing = Pointing(pointing_file)
    # get the closest
    pointing.add_closest(dist)
    pointing.mosaic()
    pointing.cutout()
