#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2017 - Francesco de Gasperin
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

# Usage: flag_weight_to_zero.py vis.MS
# Set to 0 the weight of flagged data - useful for dysco

import os, sys
import numpy as np
import pyrap.tables as pt

msfile = sys.argv[1]

# open input/output MS
ms = pt.table(msfile, readonly=False, ack=False)
 
all_weights = ms.getcol('WEIGHT_SPECTRUM')
print('Sum weights before: %f' % np.sum(all_weights))
all_flags = ms.getcol('FLAG')
print('Resetting %f values' % np.sum(all_flags))
all_weights[all_flags] = 0
print('Sum weights after: %f' % np.sum(all_weights))
ms.putcol('WEIGHT_SPECTRUM', all_weights)

# needs recent casacore
#pt.taql("update $ms set WEIGHT_SPECTRUM[FLAG]=0")

ms.close()
