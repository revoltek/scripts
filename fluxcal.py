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

# fluxcal.py flux (in MHz)

import sys, os
import numpy as np

print "Flux following Scaife and Heald 2012 (for A-team see source code)."
freq = float(sys.argv[1]) # in MHz

par = {'3C48':[ 64.768,-0.387,-0.420,0.181],
       '3C147':[66.738,-0.022,-1.012,0.549],
       '3C196':[83.084,-0.699,-0.110],
       '3C286':[27.477,-0.158,0.032,-0.180],
       '3C295':[97.763,-0.582,-0.298,0.583,-0.363 ],
       '3C380':[77.352,-0.767],
       'VirA':[1226.,-0.79], # deGasperin 2012
       'CygA':[10690, -0.67, -0.240, 0.021], # McKean 2016
       'CasA':[11733, -0.770], # Baars 77
       'TauA':[1838, -0.299], # Baars 77
       'Sun':[150e3, 3]} # http://www-old.ias.ac.in/jarch/jaa/9/225-229.pdf

def prod(sou, freq, i=1):
    """
    freq in MHz
    """
    if i == len(par[sou]): return 1
    return 10**(par[sou][i] * (np.log10(freq/150.))**i) * prod(sou, freq, i+1)

def flux(sou, freq):
    return par[sou][0]*prod(sou, freq)

print "Flux at", str(freq), 'MHz'
print "3C48: ",flux('3C48', freq)
print "3C147: ",flux('3C147', freq)
print "3C196: ",flux('3C196', freq)
print "3C286: ",flux('3C286', freq)
print "3C295: ",flux('3C295', freq)
print "3C380: ",flux('3C380', freq)
print "Virgo A: ",flux('VirA', freq)
print "Cygnus A: ",flux('CygA', freq)
print "Cassiopeia A: ~",flux('CasA', freq)
print "Taurus A: ~",flux('TauA', freq)
print "Sun (avg): ~",flux('Sun', freq)
if freq < 39.5: print "At this freq Jupiter might be >1e6 Jy."
