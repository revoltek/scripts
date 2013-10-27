#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#./whenweobserve.py datafile

import sys
import numpy as np
import matplotlib.pyplot as plt
import pyrap.measures as pm
import pyrap.quanta as qa
import lib_coordinates_mode as cm


def get_elev(t,ra,dec):
        pointing_s = me.direction('SUN')
        pointing_o = me.direction('j2000',qa.quantity(ra,'deg'),qa.quantity(dec,'deg'))
        tt = qa.quantity(t,'d')
        tt1 = me.epoch('utc',tt)
        me.doframe(tt1)
        obj_elev = me.measure(pointing_o,'azel')['m1']['value']*180/np.pi
        sun_elev = me.measure(pointing_s,'azel')['m1']['value']*180/np.pi
        return obj_elev, sun_elev

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.set_xlabel(r'Time [day]')
# toplabel in months
#ax.set_(r'')

types = np.dtype({'names':['name', 'ra', 'dec'], 'formats':['S100','S100','S100']})
data = np.loadtxt(sys.argv[1], comments='#', usecols=(0,1,2), unpack=False, dtype=types)
me = pm.measures()
position = me.position('wgs84',qa.quantity(3826600.961,'m'),qa.quantity(460953.402,'m'),qa.quantity(5064881.136,'m')) # superterp
me.doframe(position)

y=1
for obj in data:
        rah, ram, ras = obj['ra'].split(':')
        ra = cm.hmstora(float(rah), float(ram), float(ras)) # in deg
        decd, decm, decs = obj['dec'].split(':')
        dec = cm.dmstodec(float(decd), float(decm), float(decs)) # in deg
        for d in xrange(3):
                # count how many h the source is at elev>30 when sun is at elev<0
                observable = 0
                for h in xrange(24):
                        t = 2456293.5-2400000.5+d+h/24. # MJD for 1/1/2013 + the hours
                        obj_elev, sun_elev = get_elev(t,ra,dec)
                        print "RESULTS:", obj_elev, sun_elev
                        if obj_elev>30 and sun_elev<0: observable+=1
                # plot in greyscale from 0 to 12 h
                ax.plot(d,y,c=str(observable/24.))
        ax.text(1,y,obj['name'])
        y+=1

ax.set_ylim(ymin=0, ymax=y)
fig.savefig('whenweobserve.png')
fig.clf()

