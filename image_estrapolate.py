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

import sys, os
import optparse
import pyrap.images
import pyrap.images.coordinates
from scipy import odr
import numpy as np
import math

# solve memory leak
import matplotlib
matplotlib.use('Agg')
from matplotlib.cbook import report_memory
import pylab

def f(B, x):
      return B[0]*x + B[1]

def linear_fit(x, y, xerr=None, yerr=None):
    from scipy import odr
    linear = odr.Model(f)
    if xerr == None: xerr = np.ones(len(x))
    if yerr == None: yerr = np.ones(len(y))
    for i,e in enumerate(yerr):
	    if e == 0: yerr[i] = 1
    mydata = odr.Data(x, y, wd=1/xerr, we=1/yerr)
    myodr = odr.ODR(mydata, linear, beta0=[-1., 0.])
    myoutput = myodr.run()
    return(myoutput.beta[0],myoutput.beta[1],myoutput.sd_beta[0],myoutput.sd_beta[1])
    #return(myoutput.beta[0],myoutput.beta[1],myoutput.sd_beta[0],myoutput.sd_beta[1])

# if values is > 5 sigma for each map then OK, otherwise do not perform regression
def badrms(values, rmsvalues):
  for t in range(len(values)):
    if (values[t] < 5*rmsvalues[t]): 
      return False # TURN TO TRUE TO ACTIVATE
  return False

opt = optparse.OptionParser(usage="%prog images", version="%prog 0.3")
opt.add_option('-o', '--outimg', help='Output spctral index image [default = spidx.img]', default='estrap.img')
opt.add_option('-m', '--maskimg', help='Mask image tells image_spidx.py where perform linear regression')
opt.add_option('-r', '--rmsfile', help='RMSs file with one entry per line entry like "filename RMS_VALUE"')
opt.add_option('-f', '--freq', help='Frequency to estrapolate [Hz]' )
(options, imglist) = opt.parse_args()
outimg = options.outimg
rmsfile = options.rmsfile
maskimg = options.maskimg
try: freq = float(options.freq)
except: exit("ERROR: wrong final frequency")
print("Output file = "+outimg)
sys.stdout.flush()

# Read RMS values from file
if (rmsfile != None):
	print("Reading RMS file: "+rmsfile)
	try:
	    rmsdata = np.loadtxt(rmsfile, comments='#', dtype=np.dtype({'names':['file','rms'], 'formats':['S50',float]}))
	except IOError:
	    print("ERROR: error opening RMSs file, probably a wring name/format")
	    exit(1)

# Read images one by one.
values = []
frequencies = []
rmsvalues = []
for name in imglist:
    print("--- Reading file: "+name)
    try:
      image = pyrap.images.image(name)
      # workaround for getting correct axes
      index = np.where(np.array(image.coordinates().get_names())=='direction')[0][0]
      imgdata = np.array(image.getdata())
      for i in range(index):
        imgdata = imgdata[0]
      values.append(imgdata)
      frequencies.append(image.coordinates().get_referencevalue()[0])
    except:
      print("ERROR: error accessing iamges data, probably wrong name or data format")
      exit(1)

    print("Freq: ", image.coordinates().get_referencevalue()[0])
    if (rmsfile != None):
      try:
	rmsvalues.append([rms for (file, rms) in rmsdata  if file == name][0])
        print("RMS: ", [rms for (file, rms) in rmsdata  if file == name][0])
      except IndexError:
	print("ERROR: error accessing RMSs data, probably wrong names in the file")
	exit(1)
    else:
	rmsvalues.append(0)
    sys.stdout.flush()

values = np.array(values)
frequencies = np.array(frequencies)
rmsvalues = np.array(rmsvalues)
imgsizeX = values.shape[1]
imgsizeY = values.shape[2]

# var used if dozone=T
zoneval4reg=np.zeros(len(imglist))
zonepixelnum=0

# Read mask image
if (maskimg):
  print("Reading mask-image = "+maskimg)
  sys.stdout.flush()
  image = pyrap.images.image(maskimg)
  maskimgval = np.array(image.getdata()[0][0])
else:
  maskimgval = np.ones(imgsizeX*imgsizeY)
  maskimgval = maskimgval.reshape(imgsizeX, imgsizeY)

print("")
print("Performing regression", end=' ')
sys.stdout.flush()

# Make regression pixel by pixel
estrap = np.empty(imgsizeX*imgsizeY)
err = np.empty(imgsizeX*imgsizeY)
estrap = estrap.reshape(imgsizeX, imgsizeY)
err = err.reshape(imgsizeX, imgsizeY)
val4reg = np.empty(len(imglist))
for i in range(0, imgsizeX):
  for j in range(0, imgsizeY):
    val4reg=values.transpose()[j][i]
    # if the mask is set to 0 or rms check is true then skip this pixel
    if (maskimgval[i][j] == 0 or (rmsfile != None and badrms(val4reg, rmsvalues))):
      a = np.nan # TODO: set a mask, not to nan
      b = np.nan # TODO: set a mask, not to nan
      sa = np.nan # TODO: set a mask, not to nan
    else:
      yerr=0.434*rmsvalues/val4reg
      (a, b, sa, sb) = linear_fit(x=np.log10(frequencies), y=np.log10(val4reg), yerr=yerr)

    estrap[i][j] = np.power(10,(a*np.log10(freq)+b))
    err[i][j] = -100*sa/a # error in %
  print(".", end=' ')
  sys.stdout.flush()

# Write data
print("Writing spidx data.")
estrapimg = pyrap.images.image(imglist[0])
estrapimg.saveas(outimg)
estrapimg = pyrap.images.image(outimg)
estrapimg.putdata(estrap)
#estrapimg.tofits(outimg + ".fits")
del estrapimg
os.popen('patchCasaFreq '+outimg+' '+str(freq))
rmsimg = pyrap.images.image(imglist[0])
rmsimg.saveas(outimg + "-rms")
rmsimg = pyrap.images.image(outimg + "-rms")
rmsimg.putdata(err)
#rmsimg.tofits(outimg + "-rms.fits")
del rmsimg

print("Done.")
sys.stdout.flush()
