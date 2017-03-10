#!/usr/bin/python
# usage: plot_beamcorr.py <MS>

import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import casacore.tables as pt

msname = sys.argv[1]
print "Plotting beam of "+msname

# get time
t = pt.table(msname)
times = t.getcol('TIME')[::600] # every 10 min
t.close()

# get direction
t = pt.table(msname+'/POINTING')
direction = t.getcol('TARGET')[0][0]
print "Direction:", direction
t.close()

# get stations
t = pt.table(msname+'/ANTENNA')
stations = t.getcol('NAME')
t.close()

import lofar.stationresponse as st
s = st.stationresponse( msname=msname, inverse=True, useElementResponse=True, useArrayFactor=False, useChanFreq=False )
s.setDirection(direction[0],direction[1])

f_a, axa_a = plt.subplots(6, 7, sharey=True, sharex=True, figsize=(15,12))
f_p, axa_p = plt.subplots(6, 7, sharey=True, sharex=True, figsize=(15,12))
axa_a = axa_a.flatten()
axa_p = axa_p.flatten()
for idx_station in xrange(len(stations)):
    print "Working on station: %s", stations[idx_station]
    beam_matrix = []
    for time in times:
        beam_matrix_ch0 = s.evaluateStation(time=time,station=idx_station)[0,:,:] # only first chan
        beam_matrix.append(beam_matrix_ch0)

    # beam_matrix is [time:pol1:pol2]
    re = np.real(beam_matrix)
    im = np.imag(beam_matrix)

    amp = np.sqrt((re**2)+(im**2))
    axa_a[idx_station].plot(times,amp[:,0,0],'og',markersize=2, markeredgecolor='none')
    axa_a[idx_station].plot(times,amp[:,1,1],'ob',markersize=2, markeredgecolor='none')

    ph = np.arctan2(re,im)
    axa_p[idx_station].plot(times,ph[:,0,0],'og',markersize=2, markeredgecolor='none')
    axa_p[idx_station].plot(times,ph[:,1,1],'ob',markersize=2, markeredgecolor='none')

    # beauty
    axa_a[idx_station].set_xlim(xmin=np.min(times), xmax=np.max(times))
    axa_p[idx_station].set_xlim(xmin=np.min(times), xmax=np.max(times))
    axa_a[idx_station].set_xlabel('Time')
    axa_p[idx_station].set_xlabel('Time')
    axa_a[idx_station].set_ylabel('Beam amp')
    axa_p[idx_station].set_ylabel('Beam phase')
    axa_a[idx_station].set_title(stations[idx_station])
    axa_p[idx_station].set_title(stations[idx_station])

f_a.subplots_adjust(hspace=0, wspace=0)
f_p.subplots_adjust(hspace=0, wspace=0)
f_a.savefig(msname.replace('.MS','_beam_a.png'))
f_p.savefig(msname.replace('.MS','_beam_p.png'))
