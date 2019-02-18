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
import pyrap.images
import pyrap.tables
import numpy as np
import optparse

opt = optparse.OptionParser(usage="%prog images", version="%prog 1.0")
opt.add_option('-o', '--outimg', help='Output averaged image [default = averaged_weight]', default='averaged_weight')
opt.add_option('-r', '--rmsfile', help='rms file [must be in the form "image_name rms"]', default=None)
opt.add_option('-k', '--skip_keyword', help='Don\'t set the beam/freq keywords [default = True]', action="store_false")
(options, imglist) = opt.parse_args()
outimg = options.outimg
rmsfile = options.rmsfile
do_keyword = not options.skip_keyword
do_rms = False
print("Output file = "+outimg)

# Read RMS values from file
if (rmsfile != None):
    do_rms = True
    print("Reading RMSs:")
    try:
	    rmsdata = np.loadtxt(rmsfile, comments='#', dtype=np.dtype({'names':['file','rms'], 'formats':['S50',float]}))
    except IOError:
	    print("ERROR: error opening RMSs file, probably a wring name/format")
	    exit(1)

    rmsvalues=[]
    for name in imglist:
        rmsvalues.append([rms for (file, rms) in rmsdata if file == name][0])

# Read pixel data of first image.
img = pyrap.images.image(imglist[0])
if do_rms: pixels = img.getdata()/rmsvalues[0]
else: pixels = img.getdata()
# freq & beam infor to averaged them
freqs, ref1, ref2 = img.coordinates().get_referencevalue()
if do_keyword:
    bmajs = img.imageinfo()['restoringbeam']['major']['value']
    bmins = img.imageinfo()['restoringbeam']['minor']['value']
    pas = img.imageinfo()['restoringbeam']['positionangle']['value']
print("Averaging freq: "+str(freqs)+"(max:"+str(np.max(pixels))+",min:"+str(np.min(pixels))+",rms:"+str(np.std(pixels))+")")
sys.stdout.flush()

# Add pixel data of other images one by one.
for i, name in enumerate(imglist[1:]):
    tmp = pyrap.images.image(name)
    if do_rms: pixels_i = tmp.getdata()/rmsvalues[i]
    else: pixels_i = tmp.getdata()
    print("Averaging freq:",str(tmp.coordinates().get_referencevalue()[0])+"(max:"+str(np.max(pixels_i))+",min:"+str(np.min(pixels_i))+",rms:"+str(np.std(pixels_i))+")")
    pixels += pixels_i
    freqs += tmp.coordinates().get_referencevalue()[0]
    if do_keyword:
        bmajs += tmp.imageinfo()['restoringbeam']['major']['value']
        bmins += tmp.imageinfo()['restoringbeam']['minor']['value']
        pas += tmp.imageinfo()['restoringbeam']['positionangle']['value']
    del tmp
    sys.stdout.flush()

freq_mean = freqs / float(len(imglist))
if do_keyword:
    bmaj_mean = bmajs / float(len(imglist))
    bmin_mean = bmins / float(len(imglist))
    pa_mean = pas / float(len(imglist))
    print("Mean freq:",str(freq_mean))
    print("Mean beam (maj, min, pa):",str(bmaj_mean),str(bmin_mean),str(pa_mean))

# Write averaged pixel data
img.saveas(outimg + ".img")
img = pyrap.images.image(outimg + ".img")
if do_rms: 
    rmsvalues = [1/i for i in rmsvalues]
    img.putdata(pixels / sum(rmsvalues))
else:
    img.putdata(pixels / float(len(imglist)))
del img

# upgrading frequency values
os.popen('patchCasaFreq '+outimg+".img "+str(freq_mean))
if do_keyword:
    t=pyrap.tables.table(outimg+'.img', readonly=False)
    t.putkeyword('imageinfo.restoringbeam.major.value',bmaj_mean)
    t.putkeyword('imageinfo.restoringbeam.minor.value',bmin_mean)
    t.putkeyword('imageinfo.restoringbeam.positionangle.value',pa_mean)

print("done.")
sys.stdout.flush()
