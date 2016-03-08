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

#./xray-annumi-sub.py name.fits

import os, sys
import numpy as np
import pyfits, pywcs
import matplotlib.pyplot as plt
import scipy.interpolate as sp
import scipy.ndimage.measurements as sm

filename = sys.argv[1] # fits file of X-ray observation
centre = 'centroid' # in deg or 'max' or 'centroid'
ann_separation = 0.003 # in deg

# read file
hdulist = pyfits.open(filename, mode='readonly')
wcs = pywcs.WCS(hdulist[0].header)
data = hdulist[0].data
print "Image size is: ", data.shape

# create sampling point for function

# get the max
if centre == 'max': centre_p = np.squeeze(np.where(data == data.max()))
# get the centroid
elif centre == 'centroid': centre_p = sm.center_of_mass(data)
elif type(centre) is list: centre_p = wcs.wcs_sky2pix([centre], 1)[0] # centre in pixels
else: sys.exit('centre must be [coord,inates],\'max\' or \'centroid\'')

print "Central pixel is", centre_p
x, y = np.indices((wcs.naxis2, wcs.naxis1))
r = np.sqrt( (x-centre_p[0])**2 + (y-centre_p[1])**2 ) # radial distance of each pixel

# create ann_radii from the pixels size
assert abs(hdulist[0].header['CDELT1']) == abs(hdulist[0].header['CDELT2'])
print "Annului are", ann_separation/abs(hdulist[0].header['CDELT1']), 'pixels wide'
ann_radii = np.arange(0, np.max(r), ann_separation/abs(hdulist[0].header['CDELT1'])) # CDELT in in deg

# cycle on annuli
ann_means_val = []; ann_std_val = []
ann_means_radii = []; ann_std_radii = []
# at radius 0 we choose the central pixel
ann_means_radii.append(0.)
ann_means_val.append(data[centre_p[0]][centre_p[1]])
ann_std_radii.append(0.)
ann_std_val.append(0.)
for i in xrange(1,len(ann_radii)):
    s = np.where( np.logical_and(ann_radii[i-1] < r, r < ann_radii[i]) )
    if s == []: break # safety, if the max radius was too large
    ann_means_val.append(np.mean(data[s]))
    ann_std_val.append(np.std(data[s]))
    ann_means_radii.append(np.mean(r[s]))
    ann_std_radii.append(np.std(r[s]))
ann_means_radii.append(np.max(r))
ann_means_val.append(0.)
ann_std_radii.append(0.)
ann_std_val.append(0.)

# create function
print min(ann_means_radii)
f = sp.interp1d(ann_means_radii, ann_means_val, kind='linear') # 'cubic'?

# plot function
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(wspace=0)
ax = fig.add_subplot(110)
ax.set_xlabel(r'r [pixel]')
ax.set_ylabel(r'Mean')
ax.label_outer()
ax.errorbar(ann_means_radii, ann_means_val, xerr=ann_std_radii, yerr=ann_std_val)
ax.plot(ann_means_radii, f(ann_means_radii), '--')
fig.savefig('azavg_funct.pdf', bbox_inches='tight')
fig.clf()

# create 2D aximut-symmetric function
def azavg_creator(x, y):
    r = np.sqrt( (x-centre_p[0])**2 + (y-centre_p[1])**2 )
    return f(r)
azavg = np.fromfunction(azavg_creator, [wcs.naxis2, wcs.naxis1])
print azavg.shape

# subtract and write a new file
hdulist[0].data = data / azavg
hdulist.writeto(filename.replace('fits','annsub.fits'), clobber = True)
hdulist[0].data = azavg
hdulist.writeto(filename.replace('fits','annmean.fits'), clobber = True)
hdulist.close()
