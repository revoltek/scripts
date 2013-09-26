#!/usr/bin/python

import sys
import os
import numpy
out = sys.stdout
args = sys.argv[5:]
img = args[0]

# REGRID
# regrid an image based on this image increment
ia.open('~/virgo/scale/gridmodel.image')
cs1 = ia.coordsys()
incr1 = cs1.increment()
refpix1 = cs1.referencepixel()['numeric']
shape1=ia.shape()
ref1=cs1.referencevalue()
ia.close()

ia.open(img)
cs2 = ia.coordsys()
cs2.setincrement(incr1)
cs2.setreferencepixel(value=refpix1)
ref2=cs2.referencevalue()
ref2['numeric'][0] = ref1['numeric'][0]
ref2['numeric'][1] = ref1['numeric'][1]
cs2.setreferencevalue(ref2)
imrr = ia.regrid(outfile=img+'-regridded', csys=cs2.torecord(), shape=shape1, method='cubic')
ia.close()
