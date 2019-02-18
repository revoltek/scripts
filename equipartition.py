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

# calcolate the equipartition value using beck & krause 2005

import math
import numpy as py
import scipy.constants as const

# (alpha) synchrotron spectral index
a = 0.87
# Synchrotron intensity [erg/s/cm^2/Hz/sterad]
flux = 2e-03 # Jy
bmaj = 30 # arcsec
bmin = 27 # arcsec
beamarea2arcsec = 1/(bmaj*bmin*1.1331) # from beam area to arcsec**2
arcsec2sterad = 4.25e10 # from arcsec**2 to sterad
Inu_jy = flux*beamarea2arcsec*arcsec2sterad
Inu = Inu_jy*1e-23 # Jy -> cgs
# Frequency [Hz]
nu = 323.e6
# pathlenght [cm]
l = 2.7e23
# proton to electron number density ratio
K0 = 100.

# const
e = 4.80320425e-10 # [statC] elementary charge
pi = const.pi # Pi
me = const.m_e*1e3 # [g] electrion mass
c = const.c*1e2 # [cm] speed of light
Ep = 1.5033e-3 # [erg] proton spectral brake

c1 = (3*e)/(4*pi*me**3*c**5)
def c2(a):
    gamma = 2.*a+1
    return (1/4.) * c3 * (gamma + 7/3.) / ( gamma + 1) * math.gamma((3*gamma-1)/12.) * math.gamma((3*gamma+7)/12.)
c3 = math.sqrt(3.) * e**3 / ( 4*pi*me*c**2 )
c4 = 2/3.

# Magnetic field [G]
B = ( ( 4*pi * (2.*a+1) * (K0+1) * Inu * Ep**(1.-2*a) * (nu/(2*c1))**a ) / ( (2.*a-1) * c2(a) * l * c4 ) )**(1./(a+3))
print("B_eq", B*1e6, "uG")
print("B_min", B*((a+1)/2.)**(1./(a+3))*1e6, "uG")
# 2 times the magnetic field energy density times a volume of a sphere with diameter l
print("E_mag", B**2/(8*pi)*4/3.*pi*(l/2.)**3, "erg")
print("E_tot", 2*B**2/(8*pi)*4/3.*pi*(l/2.)**3, "erg")
print("p", B**2/(8*pi)/3., "dyne/cm^2")
