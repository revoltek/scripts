#!/usr/bin/python
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
pixels = img.getdata()[0][0]

pixels = ndimage.gaussian_filter(pixels, sigma=3)
pixels_s = sobel(pixels)

# Write sobel filtered pixel data
img.saveas(img_f + ".sobel")
img = pyrap.images.image(img_f + ".sobel")
img.putdata(pixels_s)
del img
