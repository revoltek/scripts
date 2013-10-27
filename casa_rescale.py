#!/usr/bin/python

import sys
import os
import numpy
out = sys.stdout
args = sys.argv[5:]
img = args[0]

# RESCALE
ia.open(img)
cs = ia.coordsys()
freq = cs.restfrequency()['value'][0]
# core
#f = lambda p, x: p[1]*x+p[0]
#pfinal=[7.02864678, -0.55034487] # computed from VLA 325mhz and 1.4ghz at core (50")
#pfinal=[6.88863054,-0.55142578] # computed from VLA 325mhz and 1.4ghz at core (75")
#flux=imstat(imagename=img,region='M87-1024-1-core.rgn')['flux'][0]
# whole source
f = lambda p, x: numpy.log10(p[0])+p[1]*x
pfinal=[1228.36475782, -0.789306499434] # computed with rescaled fit on the whole source (ref freq: 150 MHz)
flux=imstat(imagename=img,region='M87-1024-1.rgn')['flux'][0]
print "Flux: "+str(flux)
print "Expected flux: "+str(10**f(pfinal, numpy.log10(freq/150e6)))
rescale=10**(f(pfinal, numpy.log10(freq/150e6)))/flux #check the 150mhz
print "Correction: "+str(rescale)
immath(imagename=img,mode='evalexpr',outfile=img+'-rescaled',expr=str(rescale)+"*IM0")
o = open('rescale-values.txt','a')
o.write(img+' '+str(freq)+' '+str(rescale)+'\n')
o.close()
