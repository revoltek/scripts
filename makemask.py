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

# create a mask out of an object, masking everything > gs_cut

import os
import numpy
import pyfits
import pyrap.images
import sys

image_name = 'L99284-smooth_8000-4_flag_spidx.restored.corr'
mask_name  = image_name + '.mask'
os.system('cp -r ' + image_name + ' ' + mask_name)

gaus_name  = 'pybdsm_out.img'

gs_cut = 1e-6

img = pyrap.images.image(mask_name)
pixels = numpy.copy(img.getdata())
pixels_mask = 0.*numpy.copy(pixels)

#print numpy.shape(pixels_gs)
#sys.exit()
for i in xrange(len(pixels[0][0])):
  if numpy.mod(i,100) == 0:
    print i
  for j in xrange(len(pixels[0][0][i])):  
    if (pixels_gs[0,0,i,j] > gs_cut):
        pixels_mask[0,0,i,j] = 1.0

#for i in xrange(len(pixels[0][0])):
#  if numpy.mod(i,100) == 0:
#    print i
#  for j in xrange(len(pixels[0][0][i])):  
#    if numpy.isnan(pixels_mask[0,0,i,j]):
#       pixels_mask[0,0,i,j] = 0.0
    
  	
img.putdata(pixels_mask)
