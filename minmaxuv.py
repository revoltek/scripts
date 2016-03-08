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
print 'Wavelenght:', wavelenght,'m (Freq: '+str(t.getcol('REF_FREQUENCY')[0]/1.e6)+' MHz)'

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

