#!/usr/bin/python

import sys
import os
import numpy
out = sys.stdout
args = sys.argv[5:]
img = args[0]
rainc = args[1]
decinc = args[2]

# SHIFT
# shift an image

ia.open(img)
cs = ia.coordsys()
ref=cs.referencevalue()
ref['numeric'][0] += numpy.float(rainc)
ref['numeric'][1] += numpy.float(decinc)
cs.setreferencevalue(ref)
imrr = ia.regrid(outfile=img+'-shifted', csys=cs.torecord(), shape=ia.shape(), method='cubic')
ia.close()
ia.open(img+'-shifted')
ref['numeric'][0] -= numpy.float(rainc)
ref['numeric'][1] -= numpy.float(decinc)
cs.setreferencevalue(ref)
ia.setcoordsys(cs.torecord())
ia.close()
