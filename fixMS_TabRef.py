#!/usr/bin/python
#
# Script to add missing tables in the MS to conform with CASA standards 
# and to make a fix in the MeasInfo in the spectral window table.
#
# Authors: M. Iacobelli
#
# Date: Mar 04, 2016
#

import os,sys
import numpy as np, pyrap.tables as pt
from numpy import uint32

if (len(sys.argv) <= 1):
  print " \n This script will add a line in the STATE and PROCESSOR tables, fixing an issue" 
  print "in the SPECTRAL_WINDOW table that causes backward incompatibility with msoverview."
  print "To use it, please provide an MS name. \n"
  exit(0)

msin=sys.argv[1]

msin_spw=msin+"/SPECTRAL_WINDOW"
try:
   t=pt.table(msin_spw, readonly=False,ack=False)
except:
   print "Cannot find or open",msin_spw
   exit(1)

for colnm in ['CHAN_FREQ','REF_FREQUENCY']:
  tc=t.col(colnm)
  meas=tc.getkeyword('MEASINFO')
  tc.putkeyword('MEASINFO-sav', meas) #save the original keywords
  TabRefTypes_indx=[] ; TabRefCodes_indx=[]
  for indx in range(len(meas)):
    if meas.items()[indx][0] == 'TabRefTypes': TabRefTypes_indx.append(indx)
    if meas.items()[indx][0] == 'TabRefCodes': TabRefCodes_indx.append(indx)
  if len(TabRefTypes_indx) == 0:
    meas['TabRefTypes']=['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP', 'CMB'] ; print 'Column %s: Added TabRefTypes' % colnm
  if len(TabRefCodes_indx) == 0: 
    meas['TabRefCodes']=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8], dtype=uint32) ; print 'Column %s: Added TabRefCodes' % colnm
  tc.putkeyword('MEASINFO', meas)

try:
  t.flush() 
  t.close()
except:
  print "Failed to update", msin
  exit(1)

if len(TabRefTypes_indx) == 0 or len(TabRefCodes_indx) == 0:
  print "Done fixing",msin
else:
  print "None fixing",msin
