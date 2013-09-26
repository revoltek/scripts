#!/usr/bin/python

import numpy as np

TEC_low = 0.1 # low TEC val
TEC_high = 0.8 # high TEC val

# calculate max avg in freq for TEC problems (this is important in calibrating)
freq1 = 35e6
freq2 = 35.2e6
print "Freqs: ", freq1, freq2, "MHz"
print "from", ((8.44797245e9*TEC_low/freq1) - (8.44797245e9*TEC_low/freq2))*180/np.pi, "deg"
print "to", ((8.44797245e9*TEC_high/freq1) - (8.44797245e9*TEC_high/freq2))*180/np.pi, "deg"
print "This numbers must be << 100 deg"
