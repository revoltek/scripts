#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
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

print(M, "+", M_ep, "-", M_em)
