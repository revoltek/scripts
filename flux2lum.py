#!/usr/bin/python
# flux2lum.py flux (in Jy)

import sys, os
import numpy as np

flux = float(sys.argv[1]) # in Jy
dist = 45.6 # in Mpc

# 1 pc = 3.08567758e16 m
print flux * 1e-26 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2),
print "W/Hz"
print flux * 1e-26 * 1e7 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2),
print "erg/s/Hz"
print "(For a distance of "+str(dist)+" Mpc)"
