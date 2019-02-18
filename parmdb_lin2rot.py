#!/usr/bin/python

import matplotlib
matplotlib.use('GTK')
import lofar.expion.ionosphere
import lofar.parmdb
import sys
import os
import scipy
import time
import numpy as numpy
import math
import pyrap.tables
import scipy.signal

pdb = lofar.parmdb.parmdb('L71264_SAP000_SB000_uv.MS.dppp.dppp/instrument')
parms = pdb.getValuesGrid("*")
null, null, ants = numpy.array([ i.rsplit(':',2) for i in pdb.getNames() ]).transpose()
ants = set(ants)

cI = numpy.complex(0.,1.)

for ant in ants:
  print('updating ', ant)
  xx = parms['Gain:0:0:Real:'+ant]['values'] + cI * parms['Gain:0:0:Imag:'+ant]['values']
  yy = parms['Gain:1:1:Real:'+ant]['values'] + cI * parms['Gain:1:1:Imag:'+ant]['values']
  xy = parms['Gain:0:1:Real:'+ant]['values'] + cI * parms['Gain:0:1:Imag:'+ant]['values']
  yx = parms['Gain:1:0:Real:'+ant]['values'] + cI * parms['Gain:1:0:Imag:'+ant]['values'] 
  print(xx)
  print(yy)
  print(xy)
  print(yx)
  sys.exit()
  rr = xx - cI * xy + cI * yx + yy
  rl = xx + cI * xy + cI * yx - yy
  lr = xx - cI * xy - cI * yx - yy
  ll = xx + cI * xy - cI * yx - yy 
  parms['Gain:0:0:Real:'+ant]['values'][::] = rr.real
  parms['Gain:0:0:Imag:'+ant]['values'][::] = rr.imag
  parms['Gain:1:1:Real:'+ant]['values'][::] = ll.real
  parms['Gain:1:1:Imag:'+ant]['values'][::] = ll.imag
  parms['Gain:0:1:Real:'+ant]['values'][::] = rl.real
  parms['Gain:0:1:Imag:'+ant]['values'][::] = rl.imag
  parms['Gain:1:0:Real:'+ant]['values'][::] = lr.real
  parms['Gain:1:0:Imag:'+ant]['values'][::] = lr.imag

print("Writing...")
lofar.expion.parmdbmain.store_parms('test', parms, create_new = True)
