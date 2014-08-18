#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013 - Francesco de Gasperin
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

#./sdss_countmap.py sdss_file zmin zmax

# to get the sdss_file: http://skyserver.sdss3.org/dr10/en/tools/search/sql.aspx
# select top 120000 p.ra, p.dec, s.z, s.zErr
# from PhotoObj as p, Photoz as s
# where
# p.dec between 30.30 and 31.00
# and p.ra between 229.69 and 230.66
# and s.z between 0.03 and 0.2
# and s.zErr < 0.07
# and s.objID=p.objid

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pyfits, pywcs
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=73, Om0=0.27)

# parameters
sdss_file = sys.argv[1]
z = 0.1259
cluster_coords = [157.93151, 35.036615]
cluster_size = 1.03 # virial radius assumed cluster size in Mpc

def LdFromZ(z): return cosmo.luminosity_distance(z).value

types = np.dtype({'names':['ra','dec','z', 'z_err'], 'formats':[np.float,np.float,np.float,np.float]})
data = np.loadtxt(sdss_file, comments='#', unpack=False, converters={}, dtype=types, usecols=(0,1,2,3))

# create fits file (1 pixel is 100 kpc)
cell_size = cosmo.arcsec_per_kpc_proper(z).value * 100.
print "Cell size in arcesc:", cell_size
npix = int(cosmo.arcsec_per_kpc_proper(z).value * 2*cluster_size*1e3 / cell_size)
print "Number of pixel of final image:", npix, "x", npix
fits_data = np.zeros([npix,npix])
# set WCS
wcs = pywcs.WCS(naxis=2)
wcs.wcs.crval = cluster_coords
wcs.wcs.cdelt = np.array([-cell_size/3600., cell_size/3600.])
wcs.wcs.crpix = [npix/2, npix/2]
wcs.wcs.ctype = ["RA---SIN", "DEC--SIN"]
wcs.wcs.equinox = 2000.0
#wcs.wcs.set_pv([(2, 1, 45.0)])

# fill fits file
good_z = []
for s in data:
    # exclude object with bad z_err
    if s['z_err'] < 0: continue
    # exclude objects too far
    if LdFromZ(s['z']-s['z_err']) - LdFromZ(z) > cluster_size: continue
    # exclude objects too close
    if LdFromZ(z) - LdFromZ(s['z']+s['z_err']) > cluster_size: continue


    ra_pix, dec_pix = wcs.wcs_sky2pix(s['ra'],s['dec'], 1)
    ra_pix = int(np.round(ra_pix))
    dec_pix = int(np.round(dec_pix))
    if ra_pix < 0 or ra_pix >= npix or dec_pix < 0 or dec_pix >= npix: continue # galaxy too far
    fits_data[ra_pix][dec_pix] += 1
    good_z.append(s['z'])

print "Total number of good galaxies: ",len(good_z)

# put fits_data and wcs
header = wcs.to_header()
hdu = pyfits.PrimaryHDU(fits_data, header=header)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('sdss_distrib.fits', clobber=True)

# make z-distrib
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(wspace=0)
ax = fig.add_subplot(110)
ax.set_xlabel(r'z')
ax.set_ylabel(r'Number of galaxies')
ax.label_outer()
ax.hist(good_z, bins=15., color='white')
ax.axvline(z, c='r', ls='--', lw=2)
fig.savefig('z-distrib.pdf', bbox_inches='tight')
fig.clf()
