#!/usr/bin/python

import pyrap.tables as pt
import sys

t=pt.table(sys.argv[1])
pos = t.getcol('POSITION')
print '# observatory='+sys.argv[1]
print '# coordsys=XYZ'
for ant in pos:
	        print ant[0], ant[1], ant[2], 45.

