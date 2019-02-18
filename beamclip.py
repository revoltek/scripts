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

# put to 0 the small numbers in the corners of an image created with
# awimager - Bugs/comments: mrbell@mpa-garching.mpg.de
# USAGE: ./beamclip.py image beamimage clip_level

import sys, os
import numpy as np
import pyrap.images
import optparse

usage = "usage: %prog [options] image"
opt = optparse.OptionParser(usage)
opt.add_option('-l', '--clip_level', help="Beam level below which to clip the image. [0.001]", default=0.001, type='float')
opt.add_option('-v', '--clip_value', help="Value to insert in the clipped region. [0.]", default=0., type='float')
opt.add_option('-b', '--beam_image', help="Name of the beam image to use. Uses the \'spheroid_cut_im\' with the same root as the input image by default.", default=None, type='string')
opts, args = opt.parse_args()

cutoff = opts.clip_level
val = opts.clip_value

# List of images to average.
image = args[0]

if opts.beam_image is None:
	ndx=image.find('img')
	beam = image[:ndx]+'img0.spheroid_cut_im'
else:
	beam = opts.beam_image

# Read pixel data of first image.
img = pyrap.images.image(image)
pixels = img.getdata()[0,0]
bm = pyrap.images.image(beam)
pb = bm.getdata()[0,0]

imX = pixels.shape[0]
pbX =  pb.shape[0]
delta = pbX - imX

pixels[pb[delta/2:len(pb[0])-delta/2,delta/2:len(pb[1])-delta/2]<cutoff] = val

# Write averaged pixel data
img.putdata(pixels)
del img
