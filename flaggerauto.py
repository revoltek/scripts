import numpy
import pyrap.tables as pt
import optparse

#This script flags all autocorrelations in a measurement set.

op = optparse.OptionParser()
op.add_option('--msid','-m',action='store',help='MS name',type='string')
o, arg = op.parse_args()

t = pt.table(o.msid,readonly=False)
t1 = t.query('ANTENNA1 = ANTENNA2', columns='ANTENNA1, ANTENNA2, FLAG')

for t2 in t1.iter(["ANTENNA1","ANTENNA2"]):
 ant1 = t2.getcell("ANTENNA1",0)
 ant2 = t2.getcell("ANTENNA2",0)
 print 'Baseline = '+str(ant1)+' - '+str(ant2)
 flag = t2.getcol("FLAG")
 ones = numpy.ones(t2.getcol("FLAG")[0].shape,dtype='bool')
 for x in range(len(flag)):
  t2.putcell(columnname="FLAG",rownr=x,value=ones)

print 'Flagging of Autocorrelations Complete'
