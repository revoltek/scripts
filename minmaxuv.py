#!/usr/bin/python

from pyrap.tables import *
import sys
import numpy as np
import math

c = 299792458

t = table(sys.argv[1])
col = t.getcol('UVW')
maxdist = 0; mindist=np.inf

t = table(sys.argv[1]+'/SPECTRAL_WINDOW')

wavelenght = c/t.getcol('REF_FREQUENCY')[0]
print 'Wavelenght:', wavelenght,'m (Freq: '+str(t.getcol('REF_FREQUENCY')[0]/1e6)+' MHz)'

for u,v,w in col:
        dist = math.sqrt(u*u+v*v)
        if dist > maxdist: maxdist = dist
        if dist < mindist and dist != 0.0: mindist = dist

print 'MaxUVdist (wavelenght): ', maxdist/wavelenght
print 'MaxUVdist (meters): ', maxdist
print 'MinUVdist (wavelenght): ', mindist/wavelenght
print 'MinUVdist (meters): ', mindist
print 'Max Scale: ~',wavelenght/mindist*(180/np.pi), 'deg'
print 'Min scale: ~',wavelenght/maxdist*(180/np.pi)*3600, 'arcsec'

