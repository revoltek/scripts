#!/usr/bin/env python

from astropy.io import fits
import sys, os
import argparse
import numpy as np

class Ms(object):
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

def main(outfile, ms_files, clobber=False):
    
    # create ms objects
    mss_list = []
    for ms_file in ms_files:
        mss_list.append(Ms(ms_file))

        assert mss_list[0].nx == mss_list[-1].nx # assume same x size
        assert mss_list[0].ny == mss_list[-1].ny # assume same y size
        # not true in many wsclean cases!
        #assert mss_list[0].deltafreq == mss_list[-1].deltafreq # assume same freq size

    # Order by incresing freq
    freqs = [ms.freq for ms in mss_list]
    mss_list = [ms for (freq, ms) in sorted(zip(freqs, mss_list))]
    freqs = sorted(freqs)
    print "Frequencies:", freqs

    # Get cube dimensions from first image
    nx = mss_list[0].nx
    ny = mss_list[0].ny
    nv = len(mss_list)

    # Create the cube, only stokes I
    cube = np.empty([1, nv, ny, nx], dtype=float)
    
    for i, ms in enumerate(mss_list):
        cube[0,i,:,:] = ms.data

    hdulist = fits.PrimaryHDU(cube)
    hdulist.header = mss_list[0].header.copy()
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
