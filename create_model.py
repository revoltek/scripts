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
#
# Create a fake skymodel, useful for tests

import numpy as np
from lib_coordinates_mode import *

# number of point sourcer
Nsources = 10
# rmax of the spiral in deg
Rmax = 5.
# number of angles
Nang = 8.

# Center of the spiral
#virgo
#12:30:49.420000
rah=12.
ram=30.
ras=49.42
#+012.23.28.000000
decd=12.
decm=23.
decs=28.

# print necessary lines
print "# (Name, Type, Ra, Dec, I, Q, U, V, ReferenceFrequency='74e6', SpectralIndex='[]') = format"
print "C0, POINT, 12:30:49.42, +12.23.28.00, 10.0, 0.0, 0.0, 0.0, 136.0e+06, [0]"
print "center, POINT, 12:30:49.42, +12.23.28.00, 0.0, 0.0, 0.0, 0.0, 136.0e+06, [0]"

radeg = hmstora(rah,ram,ras)
decdeg = dmstodec(decd,decm,decs)

Rinc = Rmax / float(Nsources)
for i in xrange(Nsources):
	# compute radii values
	R = (i+1)*Rinc
	# compute angles value
	ang = i*(2*np.pi)/Nang
	# convert from polar to cartesian
	X=R*np.cos(ang)
	Y=R*np.sin(ang)
	# convert to RA,DEC
	(rah,ram,ras) = ratohms(radeg+X)
	(decd,decm,decs) = dectodms(decdeg+Y)
    	ra = str(rah)+":"+str(ram)+":"+str(ras)
    	dec = "+"+str(decd)+"."+str(decm)+"."+str(decs)
    	print "C"+str(i+1)+", POINT, "+str(ra)+", "+str(dec)+", 1.0, 0.0, 0.0, 0.0, 136.0e+06, [0]"
