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

import os, sys
import casacore.tables as pt
from astropy.time import Time

def get_timestep(ms):
    with pt.table(ms, ack = False) as t:
        times = sorted(set(t.getcol('TIME')))
    print("%s: Time step %i seconds (total timesteps: %i)." % (ms, times[1]-times[0], len(times)))

def get_freq(ms):
    """
    Get chan frequencies in Hz
    """
    with pt.table(ms + "/SPECTRAL_WINDOW", ack = False) as t:
        freqs = t.getcol("CHAN_FREQ")[0] * 1e-6 # MHz
        nchan = t.getcol("NUM_CHAN")[0]
    min_freq = min(freqs)
    max_freq = max(freqs)
    interval_freq = (max_freq-min_freq)/nchan
    print("%s: Freq range: %.3f MHz - %.3f MHz, interval: %.3f MHz (total channels: %i)" % (ms,min_freq,max_freq,interval_freq,nchan) )


for ms in sys.argv[1:]:
    get_timestep(ms)
    get_freq(ms)
