#!/usr/bin/python
# usage: plot_weights.py xxx.MS

import os, sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
#sys.path = ['/home/dijkema/opt/python-casacore/lib/python2.7/site-packages/']
from casacore.tables import taql

ms = sys.argv[1]

print "selecting on chan"
#t=taql('select ANTENNA, gmeans(DATA) as DATA, gmeans(WEIGHT_SPECTRUM) as WEIGHT from [[ SELECT ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM from '+ms+'], [SELECT ANTENNA2 AS ANTENNA, DATA, WEIGHT_SPECTRUM from '+ms+']] group by ANTENNA')
t1 = taql('select ANTENNA1, DATA, WEIGHT_SPECTRUM, TIME, FIELD_ID from '+ms+' where any(FLAG)==False')
t = taql('select ANTENNA, MEANS(GAGGR(DATA), 0) as DATA, MEANS(GAGGR(WEIGHT_SPECTRUM), 0) as WEIGHT from [SELECT ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM from $t1] group by ANTENNA')
weights = t.getcol('WEIGHT')[:,:,0]
data = np.abs(t.getcol('DATA')[:,:,0])

fig, ax1 = plt.subplots()
colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
#ax1.set_yscale('log')
for w, c in zip(weights, colors):
    ax1.scatter(xrange(len(w)),w,color=c,marker='.')
ax1.set_xlim(xmin=0,xmax=len(w))
ax1.set_ylim(ymin=np.min(weights),ymax=np.max(weights))
ax1.set_xlabel('Channel')
ax1.set_ylabel('Weights')
plt.legend()
plt.savefig('weightVSchan.png', bbox_inches='tight')

fig, ax1 = plt.subplots()
colors = iter(cm.rainbow(np.linspace(0, 1, len(data))))
#ax1.set_yscale('log')
for d, c in zip(data, colors):
    ax1.scatter(xrange(len(d)),d,color=c,marker='.')
ax1.set_xlim(xmin=0,xmax=len(d))
ax1.set_ylim(ymin=np.min(data),ymax=np.max(data))
ax1.set_xlabel('Channels')
ax1.set_ylabel('Data')
plt.legend()
plt.savefig('dataVSchan.png', bbox_inches='tight')


print "selecting on time"
t=taql('select TIME, MEANS(GAGGR(WEIGHT_SPECTRUM), 0) as WEIGHT, mscal.azel1()[1] as ELEV from $t1 group by TIME')
time=t.getcol('TIME')
elev=t.getcol('ELEV')
weights=t.getcol('WEIGHT')[:,0,0]

fig, ax1 = plt.subplots()
ax1.plot(time,weights,'k.')
ax2 = ax1.twinx()
ax2.plot(time,elev*180/np.pi,label='elevation')
ax1.set_xlabel('Time')
ax1.set_ylabel('Weights')
#ax1.set_yscale('log')
ax1.set_ylim(ymin=np.min(weights),ymax=np.max(weights))
ax2.set_ylabel('Elevation (deg)')
plt.legend()
plt.savefig('weightVStime.png', bbox_inches='tight')
