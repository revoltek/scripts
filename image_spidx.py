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

# in case of pyrap RuntimeError: Failed Always Assert "coords[nc-1] != 0"
# try to convert the map in fits, and then back to MS using pyrap "saveas"

import sys, os
import optparse
import numpy as np
import pyrap.images
import pyrap.tables
from progressbar import ProgressBar
from linearfit import *

opt = optparse.OptionParser(usage="%prog images", version="%prog 0.3")
opt.add_option('-o', '--outimg', help='Output spctral index image [default = spidx.img]', default='spidx.img')
opt.add_option('-m', '--fluxmask', help='Mask image tells image_spidx.py where perform linear regression', default=None)
opt.add_option('-e', '--fluxerr', help='A percentage error to be added to the flux [default=0]', default=0)
opt.add_option('-r', '--rmsmask', help='RMSs mask', default=None)
opt.add_option('-s', '--Nsigma', help='Remove a pixel if it is below this sigma in any image [default = 1]', default=1, type='float')
opt.add_option('-k', '--rescale', help='File with rescaling values to multiply [imgname rescaling] [default = None]', default=None)
(options, imglist) = opt.parse_args()
if len(imglist) < 2: sys.exit('Not enough images provided')
outimg = options.outimg
fluxmask = options.fluxmask
fluxerr = options.fluxerr
rmsmask = options.rmsmask
n_sigma = options.Nsigma
rescalefile = options.rescale
print "Output file: "+outimg

# Read images one by one.
values = []
frequencies = []
rmsvalues = []
for name in imglist:
    print "--- Reading file: "+name
    image = pyrap.images.image(name)
    # workaround for getting correct axes
    index = np.where(np.array(image.coordinates().get_names())=='direction')[0][0]
    val = np.array(image.getdata())
    # images mask have 0 for good data and 1 for bad! we use the opposite!
    maskval = ~np.array(image.getmask())

    # rescaling if required
    if rescalefile != None:
        rescaleData = np.loadtxt(rescalefile, comments='#', dtype=np.dtype({'names':['img','rescale'], 'formats':['S100',float]}))
        if name in rescaleData['img']:
            r = rescaleData['rescale'][np.where(rescaleData['img'] == name)][0]
            val *= r
            print "Rescaling", name, "flux by", r

    # remove nested useless axes
    for i in xrange(index):
        val = val[0]
        maskval = maskval[0]
    values.append(val)

    # get mask from image file and fluxmask file
    if fluxmask != None:
        image = pyrap.images.image(fluxmask)
        fluxmaskval = np.array(image.getdata()[0][0])
        del image
        # Merge given mask with images masks
        maskval = np.logical_and(fluxmaskval, maskval)

    # get frequencies
    t = pyrap.tables.table(name, ack=False)
    try:
        frequencies.append(t.getkeyword('coords.worldreplace2')[0])
    except:
        #frequencies.append(image.coordinates().get_referencevalue()[0])
        frequencies.append(t.getkeyword('coords.worldreplace1')[0])
    t.close()
    print "Freq: ", frequencies[-1]

    # mask negative values
    maskval[np.where(val < 0.)] = 0

    # Read RMS values from rmsmask file
    if rmsmask != None:
        image = pyrap.images.image(rmsmask)
        rmsmaskval = np.array(image.getdata()[0][0])
        del image
        linearized_val = val[np.where(rmsmaskval == 1)]
        rms = np.sqrt(sum(n*n for n in linearized_val)/len(linearized_val))
        print "Rms:", rms
        maskval[np.where(val < n_sigma*rms)] = 0
        rmsvalues.append(rms)
    else:
        rmsvalues.append(0.)

    # initialize global mask
    if name == imglist[0]: globalmaskval = maskval
    else: globalmaskval = np.logical_and(globalmaskval, maskval)

values = np.array(values)
frequencies = np.array(frequencies)
rmsvalues = np.array(rmsvalues)
n, imgsizeX, imgsizeY = values.shape

# Make regression pixel by pixel
spidx = np.empty(imgsizeX*imgsizeY)
err = np.empty(imgsizeX*imgsizeY)
spidx = spidx.reshape(imgsizeX, imgsizeY)
err = err.reshape(imgsizeX, imgsizeY)
val4reg = np.empty(len(imglist))
pbar = ProgressBar(maxval=imgsizeX)
for i in range(0, imgsizeX):
  for j in range(0, imgsizeY):
    val4reg=values.transpose()[j][i]
    # if the mask is set to 0 or rms check is true then skip this pixel
    if globalmaskval[i][j] == False:
      a = np.nan
      sa = np.nan
    else:
      # err of log of quadrature sum of rms and fluxscale-error
      yerr=0.434*np.sqrt(rmsvalues**2+(fluxerr*val4reg/100)**2)/val4reg
      (a, b, sa, sb) = linear_fit(x=np.log10(frequencies), y=np.log10(val4reg), yerr=yerr)
#      (a, b, sa, sb) = linear_fit_odr(x=np.log10(frequencies), y=np.log10(val4reg), yerr=yerr)
#      if i == 100 and j == 100:
#          plotlogax({'flux':val4reg, 'freq':frequencies, 'rms':rmsvalues}, 'testodr.png')

    spidx[i][j] = a
    #err[i][j] = -100*sa/a # error in %
    err[i][j] = sa # error in std dev
  pbar.update(i)
pbar.finish()

print "\n\n"
  
# Write data (go back to mask convenction: 1=bad, 0=good)
print "Writing spidx data."
spidximg = pyrap.images.image(imglist[-1])
spidximg.saveas(outimg)
spidximg = pyrap.images.image(outimg)
spidximg.putdata(spidx)
spidximg.putmask(~globalmaskval)
#spidximg.tofits(outimg + ".fits")
del spidximg
rmsimg = pyrap.images.image(imglist[-1])
rmsimg.saveas(outimg + "-rms")
rmsimg = pyrap.images.image(outimg + "-rms")
rmsimg.putdata(err)
rmsimg.putmask(~globalmaskval)
#rmsimg.tofits(outimg + "-rms.fits")
del rmsimg

print "Done."
