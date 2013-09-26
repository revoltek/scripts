#!/usr/bin/python
# read a MS and then rescale all the channels to subtract
# a -0.7 spidx from all of them, the central freq is taken
# as a reference

# Usage: ./retive_spidx_sub.py dataset.MS

import pyrap.tables as pt
import sys, cmath
import numpy as np

# open MS and get channels
t = pt.table(sys.argv[1]+'/SPECTRAL_WINDOW',ack=False,readonly=True)
freq_ref = t.getcell('REF_FREQUENCY', 0)
freq_chans = t.getcell('CHAN_FREQ', 0)
print "Central freq:", freq_ref

t = pt.table(sys.argv[1],ack=False,readonly=False)
data = t.getcol('CORRECTED_DATA')
for i, freq_chan in enumerate(freq_chans):
	print "Working on chan:", freq_chan,
	# find factor
	fact = 10.**(0.7 * np.log(freq_chan/freq_ref))
	print "(factor:", fact, ")", " - ",
	amp = np.absolute(data[:,i,:])
	ph = np.angle(data[:,i,:])
	print np.mean(amp), " -> ",
	amp *= fact
	print np.mean(amp)
	data[:,i,:] = amp * np.exp(1j*ph)

# write MS
t.putcol('CORRECTED_DATA', data)
t.close()
