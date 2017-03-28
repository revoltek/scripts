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
t=taql('select ANTENNA, gmeans(DATA[FLAG]) as DATA, gmeans(WEIGHT_SPECTRUM[FLAG]) as WEIGHT from [[ SELECT ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM, FLAG from '+ms+'], [SELECT ANTENNA2 AS ANTENNA, DATA, WEIGHT_SPECTRUM, FLAG from '+ms+']] group by ANTENNA')
#weights = np.average(np.average(t.getcol('WEIGHT'),axis=0),axis=1)
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
t=taql('select TIME, gmean(WEIGHT_SPECTRUM[FLAG]) as WEIGHT, mscal.azel1()[1] as ELEV from '+ms+' group by TIME')
time=t.getcol('TIME')
elev=t.getcol('ELEV')
weights=t.getcol('WEIGHT')

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
