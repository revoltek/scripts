#!/usr/bin/python

import os, sys
import casacore.tables as pt

t = pt.table(sys.argv[1])
times = sorted(set(t.getcol('TIME')))
t.close()
os.system('msoverview in=%s' % sys.argv[1])
print("Time step %i seconds." % (times[1]-times[0]))
