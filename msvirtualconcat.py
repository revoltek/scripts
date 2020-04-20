#!/usr/bin/env python
# usage: msvirtualconcat.py ms1 ms2...
# cocnat.MS is the virtual concatenation of mss passed as arguments

import os, sys
import casacore.tables

MSlist = sys.argv[1:]
MSconcat = 'concat.MS'

casacore.tables.msutil.msconcat(MSlist, MSconcat, concatTime = False)
