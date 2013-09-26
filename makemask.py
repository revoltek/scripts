#!/usr/bin/python
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
