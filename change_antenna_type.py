#!/usr/bin/python
import sys
import pyrap.tables

msname=sys.argv[1]

t = pyrap.tables.table(msname)
t_ant = pyrap.tables.table(t.getkeyword('ANTENNA'),readonly=False)
col = t_ant.getcol( 'MOUNT')
print col
for i in range(len(col)):
	col[i] = 'ALT-AZ'
print col
t_ant.putcol('MOUNT',col)
