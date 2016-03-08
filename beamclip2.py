#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys
import numpy
import pyfits
import pyrap.images


image_name = sys.argv[1]
outname    = image_name + '.pbblk'
os.system('cp -r ' + image_name + ' ' + outname)

# image of the averaged primary beam => .avgpb
pbim_name  = sys.argv[2]
# image of the spherodical cut => .spheroid_cut_im
#spher_name = 'L99284-smooth_8000-4_flag_spidx0.spheroid_cut_im'

pb_cut = 0.03
#sp_cut = 0.03

img = pyrap.images.image(outname)
pixels = numpy.copy(img.getdata())

pb  = pyrap.images.image(pbim_name)
pixels_pb = pb.getdata()

#sp  = pyrap.images.image(spher_name)
#pixels_sp = sp.getdata()


idx1 = numpy.where(pixels_pb < pb_cut)
#idx2 = numpy.where(pixels_sp < sp_cut)

pixels[idx1] = numpy.nan
#pixels[idx2] = numpy.nan

#for i in xrange(len(pixels[0][0])):
#  if numpy.mod(i,25) == 0:
#    print i
#  for j in xrange(len(pixels[0][0][i])):
#    if (pixels_pb[0,0,i,j] < pb_cut) or (pixels_sp[0,0,i,j] < sp_cut):
#       pixels[0,0,i,j] = numpy.nan


	
img.putdata(pixels)
#img.saveas(outname)




