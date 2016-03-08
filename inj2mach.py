#!/usr/bin/python

# Compute the Mach number given the radio injection spectral index
# Usage: mach_from_inj.py alpha +alpha_err -alpha_err

import sys, os
import numpy as np

a = abs(float(sys.argv[1]))
a_ep = abs(float(sys.argv[2]))
a_em = abs(float(sys.argv[3]))

M = np.sqrt( (2*a+3) / (2*a-1) )
M_ep = abs( a_ep * (2*(2*a-1)-2*(2*a+3))/(2*(2*a-1)**2*np.sqrt((2*a+3)/(2*a-1))) )
M_em = abs( a_em * (2*(2*a-1)-2*(2*a+3))/(2*(2*a-1)**2*np.sqrt((2*a+3)/(2*a-1))) )

print M, "+", M_ep, "-", M_em
