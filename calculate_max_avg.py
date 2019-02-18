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

import numpy as np

TEC_low = 0.1 # low TEC val
TEC_high = 0.7 # high TEC val

# calculate max avg in freq for TEC problems (this is important in calibrating)
freq1 = 35e6
freq2 = 35.2e6
print("Freqs: ", freq1, freq2, "MHz")
print("from", ((8.44797245e9*TEC_low/freq1) - (8.44797245e9*TEC_low/freq2))*180/np.pi, "deg")
print("to", ((8.44797245e9*TEC_high/freq1) - (8.44797245e9*TEC_high/freq2))*180/np.pi, "deg")
print("This numbers must be << 100 deg")
