#!/usr/bin/python

# Compute the Mach number given the radio injection spectral index
# Usage: mach_from_inj.py alpha alpha_err

import sys, os
import numpy as np

a = abs(float(sys.argv[1]))
a_e = abs(float(sys.argv[2]))

M = np.sqrt( (2*a+3) / (2*a-1) )
M_e = a_e * (2*(2*a-1)-2*(2*a+3))/(2*(2*a-1)**2*np.sqrt((2*a+3)/(2*a-1)))

print M, "+/-", M_e
