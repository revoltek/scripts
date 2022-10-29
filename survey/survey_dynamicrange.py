#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2022 - Francesco de Gasperin
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

import os, sys, glob
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.io import fits
import pyregion

def isolated(tab, dist):
    # reduce to isolated sources - nothing closer than $dist arcsec
    idx_match, sep, _ = match_coordinates_sky(SkyCoord(tab['RA'], tab['DEC']),\
                                              SkyCoord(tab['RA'], tab['DEC']), nthneighbor=2)
    idx_match = np.arange(0,len(tab))[sep>dist*u.arcsec]
    return tab[idx_match]

def do_dynrng(mosaics_coord, catfile, outfile):

    with open(outfile, 'w') as f:
        f.write("#id flux rms... \n")

    cat = Table.read(catfile)
    print("Inital len:", len(cat))
    cat = isolated(cat, 560)
    print("After isolation:", len(cat))
    cat = cat[cat['Maj']<30]
    print("Final len:", len(cat))
    for s, sou in enumerate(cat):

        ra = sou['RA']
        dec = sou['DEC']
        # create region
        region_strings = ['fk5;panda(%f,%f,0,359.9,1,%f",%f",1)' % (ra,dec,dist,dist+30) for dist in np.arange(30,630,60)]
        print(region_strings)

        # debug
        #for region_string in region_strings:
        #    print(region_string[4:])

        # !!! open right mosaic
        img = args.resdir+'/'+sou['Mosaic_id']+'-mosaic.fits'
        img_data, img_hdr = fits.getdata(img, 0, header=True)

        print('(%04i/%04i) Working on %s (%s)' % (s, len(cat), sou['Source_name'], img))
        outstr = "%s %s " % (sou['Source_name'], sou['Total_flux'])
        for region_string in region_strings:
            sl = pyregion.parse(region_string)
            mask = sl.get_mask(header=img_hdr, shape=img_data.shape)
            #print(np.sum(mask))
    
            # extract rms
            rms = np.sqrt(np.nanmean(np.square(img_data[mask])))
    
            # add values
            outstr += "%s " % rms

        outstr += '\n'
        
        # save file
        with open(outfile, 'a') as f:
            f.write(outstr)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Estimate the dynamic range as a function of flux and distance from sources.')
    parser.add_argument('--resdir', type=str, help='Residuals dir.')
    parser.add_argument('--catfile', type=str, help='Source catalogue.')
    parser.add_argument('--outfile', type=str, help='Output file.')
    args = parser.parse_args()

    mosaics = glob.glob(args.resdir+'/*mosaic.fits')
    mosaics_coord = {'file':[], 'ra':[], 'dec':[]}
    print('Preparing mosaic images')
    for mosaic in mosaics:
        hdr = fits.getheader(mosaic)
        ra = hdr['CRVAL1']
        dec = hdr['CRVAL2']
        mosaics_coord['file'].append(mosaic)
        mosaics_coord['ra'].append(ra)
        mosaics_coord['dec'].append(dec)

    do_dynrng(mosaics_coord, args.catfile, args.outfile)
