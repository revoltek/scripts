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
import numpy as np

def get_obs(ms):
    with pt.table(ms + '/OBSERVATION', ack=False) as t:
        obsid =  int(t.getcell("LOFAR_OBSERVATION_ID", 0))
        pathFieldTable = ms + "/FIELD"
        nameField = (pt.taql("select NAME from $pathFieldTable")).getcol("NAME")[0]
        print("%s: Observation field: %s" % (ms, nameField))
        print("%s: Observation ID %s" % (ms, obsid))

def get_timestep(ms):
    with pt.table(ms, ack = False) as t:
        times = sorted(set(t.getcol('TIME')))
    print("%s: Time step %i seconds (total timesteps: %i)." % (ms, times[1]-times[0], len(times)))
    time = Time( times[0]/86400, format='mjd')
    print("%s: Starting time %s" % (ms, str(time.iso)))

def get_freq(ms):
    """
    Get chan frequencies in Hz
    """
    with pt.table(ms + "/SPECTRAL_WINDOW", ack = False) as t:
        freqs = t.getcol("CHAN_FREQ")[0] * 1e-6 # MHz
        nchan = t.getcol("NUM_CHAN")[0]
        chan_bandwidth = t.getcol("CHAN_WIDTH")[0][0] * 1e-6 # MHz
    min_freq = min(freqs-chan_bandwidth/2.)
    max_freq = max(freqs+chan_bandwidth/2.)
    mean_freq = np.mean(freqs)
    bandwidth = max_freq-min_freq
    print("%s: Freq range: %f MHz - %f MHz (bandwidth: %f MHz, mean freq: %f MHz)" % (ms, min_freq, max_freq, bandwidth, mean_freq))
    print("%s: Channels: %i ch (bandwidth: %f MHz)" % (ms, nchan, chan_bandwidth) )

def get_dir(ms):
    """
    Get phase centre
    """
    field_no = 0
    ant_no   = 0
    with pt.table(ms + "/FIELD", ack = False) as field_table:
        direction = field_table.getcol("PHASE_DIR")
    RA        = direction[ant_no, field_no, 0]
    Dec       = direction[ant_no, field_no, 1]

    if (RA < 0):
        RA += 2 * np.pi

    print("%s: Phase centre: %f, %f (deg)" % (ms, np.degrees(RA), np.degrees(Dec)))

    with pt.table(ms+'/FIELD', readonly=True, ack=False) as t:
        code = t.getcell('CODE',0)
    if code == '':
        with pt.table(ms+'/OBSERVATION', readonly=True, ack=False) as t:
            code = t.getcell('LOFAR_TARGET',0)[0]
    code = code.lower().replace(' ','_')
    print("%s: Field: %s" % (ms, code))

def get_cols(ms):
    """
    get non-default colummns
    """
    with pt.table(ms, ack = False) as table:
        colnames = table.colnames()

    default_colnames = ['UVW','FLAG_CATEGORY','WEIGHT','SIGMA','ANTENNA1','ANTENNA2','ARRAY_ID','DATA_DESC_ID','EXPOSURE',
                        'FEED1','FEED2','FIELD_ID','FLAG_ROW','INTERVAL','OBSERVATION_ID','PROCESSOR_ID','SCAN_NUMBER',
                        'STATE_ID','TIME','TIME_CENTROID','DATA','FLAG','WEIGHT_SPECTRUM']

    non_default_colnames =[]
    for colname in colnames:
        if colname not in default_colnames:
            non_default_colnames.append(colname)
    if len(non_default_colnames) != 0:
        print(f'{ms}: Non-default data columns: {", ".join(non_default_colnames)}.')
    else:
        print(f'{ms}: Table does only contain default columns.')

    for colname in colnames:
        if 'DATA' in colname and colname != 'DATA_DESC_ID':
            with pt.table(ms, ack=False) as t:
                kws = t.getcolkeywords(colname)
                if 'LOFAR_APPLIED_BEAM_MODE' in kws:
                    print(f'{ms}: [{colname}] Beam mode: {kws["LOFAR_APPLIED_BEAM_MODE"]} ('
                          f'{np.degrees(kws["LOFAR_APPLIED_BEAM_DIR"]["m0"]["value"]):.4f}, '
                          f'{np.degrees(kws["LOFAR_APPLIED_BEAM_DIR"]["m1"]["value"]):.4f}).')
                else:
                    print(f'{ms}: [{colname}] Beam mode: No beam applied.')


for ms in sys.argv[1:]:
    if not os.path.exists(ms):
        print("ERROR: missing ms %s" % ms)
        sys.exit()
    get_obs(ms)
    get_timestep(ms)
    get_freq(ms)
    get_dir(ms)
    get_cols(ms)
