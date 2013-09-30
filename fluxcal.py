#!/usr/bin/python
# fluxcal.py flux (in MHz)

import sys, os
import numpy as np

print "Flux following Scaife and Heald 2012."
freq = float(sys.argv[1]) # in MHz

par = {'3C48':[ 64.768,-0.387,-0.420,0.181],
       '3C147':[66.738,-0.022,-1.012,0.549],
       '3C186':[83.084,-0.699,-0.110],
       '3C286':[27.477,-0.158,0.032,-0.180],
       '3C295':[97.763,-0.582,-0.298,0.583,-0.363 ],
       '3C380':[77.352,-0.767]}

def prod(sou, freq, i=1):
    """
    freq in MHz
    """
    if i == len(par[sou]): return 1
    return 10**(par[sou][i] * (np.log10(freq/150.))**i) * prod(sou, freq, i+1)

def flux(sou, freq):
    return par[sou][0]*prod(sou, freq)

print "3C48: ",flux('3C48', freq)
print "3C147: ",flux('3C147', freq)
print "3C186: ",flux('3C186', freq)
print "3C286: ",flux('3C286', freq)
print "3C295: ",flux('3C295', freq)
print "3C380: ",flux('3C380', freq)
