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

#./sobel.py name.img

import sys, os
import pyrap.images
import pyrap.tables
from scipy import ndimage
import numpy as np

img_f = sys.argv[1].replace('/','')

# img: a 2-D array
def sobel(img):
	dx = ndimage.sobel(img, 0)  # horizontal derivative
	dy = ndimage.sobel(img, 1)  # vertical derivative
	mag = np.hypot(dx, dy)  # magnitude
	mag *= 100. / np.max(mag)  # normalize
    	return mag


img = pyrap.images.image(img_f)
pixels = np.squeeze(img.getdata())#[0][0]
print "Image shape:", pixels.shape

pixels = ndimage.gaussian_filter(pixels, sigma=3)
pixels_s = sobel(pixels)

# Write sobel filtered pixel data
img.saveas(img_f + ".sobel")
img = pyrap.images.image(img_f + ".sobel")
img.putdata(pixels_s)
del img
