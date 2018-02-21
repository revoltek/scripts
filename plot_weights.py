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
t=taql('select TIME, ANTENNA, gmeans(DATA) as DATA, gmeans(WEIGHT_SPECTRUM) as WEIGHT from [[ SELECT ANTENNA1 AS ANTENNA, TIME, DATA, WEIGHT_SPECTRUM from '+ms+' where any(FLAG)==False], [SELECT ANTENNA2 AS ANTENNA, TIME, DATA, WEIGHT_SPECTRUM from '+ms+' where any(FLAG)==False]] group by ANTENNA')
#t1 = taql('select ANTENNA1, DATA, WEIGHT_SPECTRUM, TIME, FIELD_ID from '+ms+' where any(FLAG)==False')
#t = taql('select ANTENNA, MEANS(GAGGR(DATA), 0) as DATA, MEANS(GAGGR(WEIGHT_SPECTRUM), 0) as WEIGHT from [SELECT ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM from $t1] group by ANTENNA')
weights = t.getcol('WEIGHT')[:,:,3]
data = np.abs(t.getcol('DATA')[:,:,3])
ants_n = t.getcol('ANTENNA')
ants = taql('select NAME from '+ms+'/ANTENNA')
freq = taql('select CHAN_FREQ from '+ms+'/SPECTRAL_WINDOW')[0]['CHAN_FREQ']*1e-6

print "plot weightVSchan.png"
fig, ax1 = plt.subplots(figsize=(20,10))
colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
#ax1.set_yscale('log')
for w, c, a in zip(weights, colors, ants_n):
    ax1.scatter(freq,w,color=c,marker='.',label=ants[a]['NAME'])
ax1.set_xlim(xmin=min(freq),xmax=max(freq))
ax1.set_ylim(ymin=np.min(weights),ymax=np.max(weights))
ax1.set_xlabel('Freq [MHz]')
ax1.set_ylabel('Weights')
handles, labels = ax1.get_legend_handles_labels()
lgd = ax1.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.1,1))
plt.savefig('weightVSchan.png', bbox_inches='tight', bbox_extra_artists=(lgd,))

#print "plot datatVSchan.png"
#fig, ax1 = plt.subplots(figsize=(20,10))
#colors = iter(cm.rainbow(np.linspace(0, 1, len(data))))
##ax1.set_yscale('log')
#for d, c, a in zip(data, colors, ants_n):
#    ax1.scatter(xrange(len(d)),d,color=c,marker='.', label=ants[a]['NAME'])
#ax1.set_xlim(xmin=0,xmax=len(d))
#ax1.set_ylim(ymin=np.min(data),ymax=np.max(data))
#ax1.set_xlabel('Channels')
#ax1.set_ylabel('Data')
#handles, labels = ax1.get_legend_handles_labels()
#lgd = ax1.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.1,1))
#plt.savefig('dataVSchan.png', bbox_inches='tight', bbox_extra_artists=(lgd,))


print "selecting on time"
#t=taql('select TIME, MEANS(GAGGR(WEIGHT_SPECTRUM), 0) as WEIGHT, mscal.azel1()[1] as ELEV from $t1 group by TIME')
t=taql('select TIME, WEIGHT, means(gaggr(mscal.azel1()[1]),0) as ELEV from $t group by TIME')
time=t.getcol('TIME')
elev=t.getcol('ELEV')
weights=t.getcol('WEIGHT')[:,0,0]

print "plot weightVStime.png"
fig, ax1 = plt.subplots(figsize=(20,10))
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
