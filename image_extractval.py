#!/usr/bin/python
# Extract freq flux rms from a series of images

import sys, os, math
import optparse
import numpy as np
import pyrap.images
import pyrap.tables

opt = optparse.OptionParser(usage="%prog images", version="%prog 0.1")
opt.add_option('-o', '--outfile', help='Output rms file [default = rms.txt]', default='fluxes.dat')
opt.add_option('-r', '--rmsmask', help='Mask tp compute rms [default = use all]', default=None)
opt.add_option('-f', '--fluxmask', help='Mask tp compute flux [default = use inner half]', default=None)
(options, imglist) = opt.parse_args()
outfile = options.outfile
rmsmask = options.rmsmask
fluxmask = options.fluxmask
if len(imglist) == 0: sys.exit("Missing images.")
print "Output file = "+outfile

# Read images one by one.
print "Reading files:"
values = []
rmsvalues = []
fluxvalues = []
freqs = []
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

# Read RMS mask image
if (rmsmask != None):
  print "Reading RMS mask = "+rmsmask
  image = pyrap.images.image(rmsmask)
  rmsmaskval = np.array(image.getdata()[0][0])
else:
  print "Assuming RMS mask = all"
  rmsmaskval = np.ones(values[0].shape)

# Read flux mask image
if (fluxmask != None):
  print "Reading flux mask = "+fluxmask
  image = pyrap.images.image(fluxmask)
  fluxmaskval = np.array(image.getdata()[0][0])
else:
  print "Assuming flux mask = inner quarter"
  fluxmaskval = np.zeros(values[0].shape)
  # set the inner half to be the mask
  fluxmaskval[len(fluxmaskval)/4:3*len(fluxmaskval)/4,\
        len(fluxmaskval[0])/4:3*len(fluxmaskval[0])/4] = 1

sys.stdout.flush()
print "Calculating fluxes and RMSs"
for i, img in enumerate(imglist):
  rms = np.std(values[i][np.where(rmsmaskval == 1)])
  rmsvalues.append(rms)
  t = pyrap.tables.table(img)
  bmaj = t.getkeyword('imageinfo.restoringbeam.major.value')
  bmin = t.getkeyword('imageinfo.restoringbeam.minor.value')
  assert t.getkeyword('imageinfo.restoringbeam.major.unit') == 'arcsec'
  assert t.getkeyword('imageinfo.restoringbeam.minor.unit') == 'arcsec'
  npix = len( np.array(np.where(fluxmaskval == 1)[0]) )
  dpix = abs(t.getkeyword('coords.direction0.cdelt')*(180/np.pi)*3600) # in arcsec
  assert dpix[0] == dpix[1] # increment is equal on both sides of the pixels
  flux = np.sum(values[i][np.where(fluxmaskval == 1)])/((1.1331*bmaj*bmin)/(dpix[0]**2))
  fluxvalues.append(flux)
  freq = t.getkeyword('coords.worldreplace2')
  freqs.append(freq)
  print ".",

rmsvalues = np.array(rmsvalues)
fluxvalues = np.array(fluxvalues)
freqs = np.array(freqs)
np.savetxt(outfile, np.array([freqs,fluxvalues,rmsvalues]).transpose(), fmt='%f %f %f')
  
print "Done."
