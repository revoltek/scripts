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
print 'Wavelenght:', wavelenght,'m'

count = 0
for u,v,w in col:
	dist = math.sqrt(u*u+v*v)
	if dist > maxdist: maxdist = dist
	if dist < mindist and dist != 0.0:
		mindist = dist
		#print dist, u,v,w

print 'MaxUVdist (in unit of wavelenght): ', maxdist/wavelenght
print 'MaxUVdist (in meters): ', maxdist
print 'MinUVdist (in unit of wavelenght): ', mindist/wavelenght
print 'MinUVdist (in meters): ', mindist
print 'Resolution: ~',wavelenght/maxdist*(180/np.pi)*3600, 'arcsec'
