#!/usr/bin/env python
#TobiaC 2011
import optparse
import numpy
import sys
import pyrap.tables as pt
from pyrap.quanta import quantity

def checkfile(inms):
  if inms == '':
     print 'Error: give an input MS'
     sys.exit()
  #Check feed polarization
  t = pt.table(inms,ack=False)
  feed = pt.table(t.getkeyword('FEED'),ack=False)
  poltyp=feed.getcell('POLARIZATION_TYPE',0) 
  feed.close()
  t.close() 
  if options.reverse == True:
    if poltyp[0] != 'R' and poltyp[0] != 'L':
       print "WARNING: Data is not from circularly polarized feed but I'm converting a column from circular"
  else: 
    if poltyp[0] != 'X' and poltyp[0] != 'Y':
       print "WARNING: Data is not from linearly polarized feed but I'm converting a column from linear."

def setupiofiles(inms, outms, incolumn, outcolumn):
  if outms == None:
     outms = inms
  if inms != outms :
     t = pt.table(inms)
     t.copy(outms, True, True)
     t.close()
     print "Finished copying"
  # create output column if doesn't exist
  to = pt.table(outms, readonly=False)
  if not outcolumn in to.colnames():
      ti = pt.table(inms)
      desc = ti.getcoldesc('DATA')
      desc['name'] = outcolumn
      dminfo = ti.getdminfo(incolumn)
      dminfo['NAME'] = outcolumn
      to.addcols(desc,dminfo)
      ti.close()
  to.close()
  return outms

def mslin2circ(incol, outcol, outms):
  tc = pt.table(outms,readonly=False)
  dataXY = tc.getcol(incol)
  I=numpy.complex(0.0,1.0)
  dataRL = 0.5* numpy.transpose(numpy.array([
           +dataXY[:,:,0]-I*dataXY[:,:,1]+I*dataXY[:,:,2]+dataXY[:,:,3],
           +dataXY[:,:,0]+I*dataXY[:,:,1]+I*dataXY[:,:,2]-dataXY[:,:,3],
           +dataXY[:,:,0]-I*dataXY[:,:,1]-I*dataXY[:,:,2]-dataXY[:,:,3],
           +dataXY[:,:,0]+I*dataXY[:,:,1]-I*dataXY[:,:,2]+dataXY[:,:,3]]),
           (1,2,0))
  tc.putcol(outcol,dataRL)

  # Merge flags
  flagXY = tc.getcol('FLAG')
  print "Initial flags:", numpy.count_nonzero(flagXY)
  for time in xrange(flagXY.shape[0]):
      for chan in xrange(flagXY.shape[1]):
          flagXY[time][chan] = numpy.count_nonzero(flagXY[time][chan]) > 0
  print "Final flags:", numpy.count_nonzero(flagXY)
  tc.putcol('FLAG',flagXY)

  #Change metadata information to be circular feeds
  feed = pt.table(tc.getkeyword('FEED'),readonly=False)
  for tpart in feed.iter('ANTENNA_ID'):
      tpart.putcell('POLARIZATION_TYPE',0,['R','L'])

  polariz = pt.table(tc.getkeyword('POLARIZATION'),readonly=False)
  polariz.putcell('CORR_TYPE',0,[5,6,7,8])
  tc.close()

def mscirc2lin(incol, outcol, outms):
  tc = pt.table(outms,readonly=False)
  dataRL = tc.getcol(incol)
  I=numpy.complex(0.0,1.0)
  dataXY = 0.5* numpy.transpose(numpy.array([
              +dataRL[:,:,0]+dataRL[:,:,1]+dataRL[:,:,2]+dataRL[:,:,3],
           I*(+dataRL[:,:,0]-dataRL[:,:,1]+dataRL[:,:,2]-dataRL[:,:,3]),
           I*(-dataRL[:,:,0]-dataRL[:,:,1]+dataRL[:,:,2]+dataRL[:,:,3]),
              +dataRL[:,:,0]-dataRL[:,:,1]-dataRL[:,:,2]+dataRL[:,:,3] ]),
           (1,2,0))
  tc.putcol(outcol,dataXY)

  # Merge flags
  flagXY = tc.getcol('FLAG')
  print "Initial flags:", numpy.count_nonzero(flagXY)
  for time in xrange(flagXY.shape[0]):
      for chan in xrange(flagXY.shape[1]):
          flagXY[time][chan] = numpy.count_nonzero(flagXY[time][chan]) > 0
  print "Final flags:", numpy.count_nonzero(flagXY)
  tc.putcol('FLAG',flagXY)

  #Change metadata information to be circular feeds
  feed = pt.table(tc.getkeyword('FEED'),readonly=False)
  for tpart in feed.iter('ANTENNA_ID'):
      tpart.putcell('POLARIZATION_TYPE',0,['X','Y'])

  polariz = pt.table(tc.getkeyword('POLARIZATION'),readonly=False)
  polariz.putcell('CORR_TYPE',0,[9,10,11,12])
  tc.close()


def updatehistory(outms):
  #Update history to show that this script has modified original data
  tc = pt.table(outms,readonly=False)
  th = pt.table(tc.getkeyword('HISTORY'), readonly=False, ack=False)
  nr=th.nrows()
  th.addrows(1)
  tr=th.row()
  tr.put(nr,{'TIME': quantity('today').get('s').get_value(), 'OBSERVATION_ID':0,'MESSAGE': ' ', 'PRIORITY': ' ', 'ORIGIN': ' ','OBJECT_ID':0, 'APPLICATION':'mslin2circ','CLI_COMMAND':[''],'APP_PARAMS': ['']})


opt = optparse.OptionParser()
opt.add_option('-i','--inms',help='Input MS (format: ms:COLUMN, default column: DATA)',default='')
opt.add_option('-o','--outms',help='Output MS (format: ms:COLUMN, default ms: InputMS, default column: DATA)')
opt.add_option('-r','--reverse',action="store_true",default=False,help='Convert from circular to linear')
options, arguments = opt.parse_args()

if len( options.inms.split(':') ) == 2:
    incolumn = options.inms.split(':')[1]
else:
    incolumn = 'DATA'

if len( options.outms.split(':') ) == 2:
    outcolumn = options.outms.split(':')[1]
else:
    outcolumn = 'DATA'

inms = options.inms.split(':')[0]
outms = options.outms.split(':')[0]
checkfile(inms)
outms = setupiofiles(inms, outms, incolumn, outcolumn)

print "INFO: inms: "+inms+" (column: "+incolumn+")"
print "INFO: outms: "+outms+" (column: "+outcolumn+")"

if options.reverse == True:
   mscirc2lin(incolumn, outcolumn, outms)
else:
   mslin2circ(incolumn, outcolumn, outms)
updatehistory(outms)
