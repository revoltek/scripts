#!/usr/bin/python
# put to 0 the small numbers in the corners of an image created with awimager

import sys, os
import numpy as np
import pyrap.images
import pyrap.tables

image = sys.argv[1]

img = pyrap.images.image(image)
pixels = img.getdata()

blanked = 0
tot = 0
for i in xrange(len(pixels[0][0])):
	sys.stdout.write('\r')
	sys.stdout.write("[%-20s] %d%%" % ('='*(i*20/len(pixels[0][0])), i*100/len(pixels[0][0])))
	sys.stdout.flush()
	for j in xrange(len(pixels[0][0][i])):
		tot += 1
		if abs(pixels[0,0,i,j]) < 1e-7:
			pixels[0,0,i,j] = 0.0
			blanked += 1

print "\nBalnked "+str(blanked*100/tot)+"% of pixels"

# Write pixel data
try:
	img.putdata(pixels)
	del img
except:
	print "WARNING: Using output.fits as output file"
	os.remove(image)
	img.saveas('output.img')
	img = pyrap.images.image('output.img')
	img.putdata(pixels)
	img.tofits(image)
	del img
	import shutil
	shutil.rmtree('output.img')

print "done."
sys.stdout.flush()
