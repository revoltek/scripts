#!/usr/bin/python
# Extract rms from a series of images, used for image_average

import sys
import optparse
import pyrap.images
import pyrap.images.coordinates
from scipy import odr
import numpy as np
import math

opt = optparse.OptionParser(usage="%prog images", version="%prog 0.1")
opt.add_option('-o', '--outfile', help='Output rms file [default = rms.txt]', default='rms.txt')
opt.add_option('-m', '--mask', help='Mask tp compute rms [default = use inner quarter]', default=None)
(options, imglist) = opt.parse_args()
outfile = options.outfile
mask = options.mask
print "Output file = "+outfile
sys.stdout.flush()

# Read images one by one.
print "Reading files:"
sys.stdout.flush()
values = []
rmsvalues = []
for name in imglist:
    try:
      image = pyrap.images.image(name)
      values.append(np.array(image.getdata()[0][0]))
      print ".",
      sys.stdout.flush()
    except:
      print "ERROR: error accessing iamges data, probably wrong name or data format"
      exit(1)

print ""
values = np.array(values)
imgsizeX = values.shape[1]
imgsizeY = values.shape[2]

# Read mask image
if (mask != None):
  print "Reading mask = "+mask
  image = pyrap.images.image(mask)
  maskval = np.array(image.getdata()[0][0])
else:
  print "Assuming mask = inner quarter"
  maskval = np.zeros(values[0].shape)
  # set the inner quarter to be the mask
  maskval[len(maskval)/4:3*len(maskval)/4,len(maskval[0])/4:3*len(maskval[0])/4] = 1

sys.stdout.flush()
# Make regression pixel by pixel
print "Calculating RMSs"
for t in xrange(len(imglist)):
  rms = np.std(values[t][np.where(maskval == 1)])
  rmsvalues.append(rms)
  print ".",
  sys.stdout.flush()

rmsvalues = np.array(rmsvalues)
np.savetxt(outfile, np.array([imglist,rmsvalues]).transpose(), fmt='%s %s')
  
print "Done."
