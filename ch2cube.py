#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from astropy.io import fits
import sys, os
import argparse
import numpy as np

class Fits(object):
    def __init__(self, filename):
        self.filename = filename
        self.data, self.header = fits.getdata(filename, 0, header=True)
        # first two axes must be ra/dec
        assert (('RA' in self.header['CTYPE1']) and ('DEC' in self.header['CTYPE2'])) or \
               (('RA' in self.header['CTYPE2']) and ('DEC' in self.header['CTYPE1']))

        self.nx = self.header['NAXIS1']
        self.ny = self.header['NAXIS2']
        
        if self.header['CTYPE3'] == 'FREQ':
            self.freq = self.header['CRVAL3']
            self.deltafreq = self.header['CDELT3']
        elif self.header['CTYPE4'] == 'FREQ':
            self.freq = self.header['CRVAL4']
            self.deltafreq = self.header['CDELT4']
        else:
            raise "Incompatible axes format in %s." % filename

def main(outfile, fitsfiles, clobber=False):
    
    # create fits objects
    fits_list = []
    for fitsfile in fitsfiles:
        fits_list.append(Fits(fitsfile))

        assert fits_list[0].nx == fits_list[-1].nx # assume same x size
        assert fits_list[0].ny == fits_list[-1].ny # assume same y size
        # not true in many wsclean cases!
        #assert fits_list[0].deltafreq == fits_list[-1].deltafreq # assume same freq size

    # Order by incresing freq
    freqs = [fits.freq for fits in fits_list]
    fits_list = [fits for (freq, fits) in sorted(zip(freqs, fits_list))]
    freqs = sorted(freqs)
    print("Frequencies:", freqs)

    # Get cube dimensions from first image
    nx = fits_list[0].nx
    ny = fits_list[0].ny
    nv = len(fits_list)

    # Create the cube, only stokes I
    cube = np.empty([1, nv, ny, nx], dtype=float)
    
    for i, fits in enumerate(fits_list):
        cube[0,i,:,:] = fits.data

    hdulist = fits.PrimaryHDU(cube)
    hdulist.header = fits_list[0].header.copy()
    hdulist.writeto(outfile, clobber=clobber)
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('filelist', help="List of fits files to combine into a cube.", 
                        nargs=argparse.REMAINDER)
    parser.add_argument('-o', '--outfile', 
                        help="Output cube name.", required=True)
    parser.add_argument('-c', '--clobber', 
                        help="Ovserwrite output file?",
                        action='store_true', default=False)

    args = parser.parse_args()
    
    main(args.outfile, args.filelist, args.clobber)
