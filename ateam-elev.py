#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin & Bas Van Der Tol
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

# Usage: ./ateam-dist.py RA (00h00m00s) DEC (+/-00d00m00s)
# Example: ./ateam-dist.py 12h30m49.4s +12d23m28s
# calculate the angular distance between an object and ateams/calibrators


from pylab import *
import pyrap.quanta as qa
import pyrap.tables as pt
import pyrap.measures as pm
import sys
import numpy

targets = [ {'name' : 'CasA', 'ra' : 6.123487680622104,  'dec' : 1.0265153995604648},
            {'name' : 'CygA', 'ra' : 5.233686575770755,  'dec' : 0.7109409582180791},
            {'name' : 'TauA', 'ra' : 1.4596748493730913, 'dec' : 0.38422502335921294},
            {'name' : 'VirA', 'ra' : 3.276086511413598,  'dec' : 0.21626589533567378},
            {'name' : 'SUN'},
            {'name' : 'JUPITER'}]

if len(sys.argv) == 2:
   msname = sys.argv[1]
else:
   print("Usage")
   print("   plot_Ateam_elevation.py <msname>")
   sys.exit()
   

# Create a measures object
me = pm.measures()

# Open the measurement set and the antenna and pointing table

ms = pt.table(msname)  

# Get the position of the first antenna and set it as reference frame
ant_table = pt.table(msname + '/ANTENNA')  
ant_no = 0
pos = ant_table.getcol('POSITION')
x = qa.quantity( pos[ant_no,0], 'm' )
y = qa.quantity( pos[ant_no,1], 'm' )
z = qa.quantity( pos[ant_no,2], 'm' )
position =  me.position( 'wgs84', x, y, z )
me.doframe( position )
#print position
ant_table.close()


# Get the first pointing of the first antenna
field_table = pt.table(msname + '/FIELD')
field_no = 0
direction = field_table.getcol('PHASE_DIR')
ra = direction[ ant_no, field_no, 0 ]
if ra<0: ra += 2*numpy.pi
dec = direction[ ant_no, field_no, 1 ]
targets.insert(0, {'name' : 'Pointing', 'ra' : ra, 'dec' : dec})
print(("Target ra/dec (deg):", targets[0]['ra']*180/numpy.pi, targets[0]['dec']*180/numpy.pi))
print(targets)
field_table.close()


# Get a ordered list of unique time stamps from the measurement set
time_table = pt.taql('select TIME from $1 orderby distinct TIME', tables = [ms])
time = time_table.getcol('TIME')
time1 = time/3600.0
time1 = time1 - floor(time1[0]/24)*24


clf()

ra_qa  = qa.quantity( targets[0]['ra'], 'rad' )
dec_qa = qa.quantity( targets[0]['dec'], 'rad' )
pointing =  me.direction('j2000', ra_qa, dec_qa)

for target in targets:
   
   t = qa.quantity(time[0], 's')
   t1 = me.epoch('utc', t)
   me.doframe(t1)

   if 'ra' in list(target.keys()):
      ra_qa  = qa.quantity( target['ra'], 'rad' )
      dec_qa = qa.quantity( target['dec'], 'rad' )
      direction =  me.direction('j2000', ra_qa, dec_qa)
      print((ra_qa, dec_qa))
   else :
      direction =  me.direction(target['name'])
      
   # Loop through all time stamps and calculate the elevation of the pointing
   el = []
   for t in time:
      t_qa = qa.quantity(t, 's')
      t1 = me.epoch('utc', t_qa)
      me.doframe(t1)
      a = me.measure(direction, 'azel')
      elevation = a['m1']
      el.append(elevation['value']/pi*180)
   el = numpy.array(el)

   plot(time1, el, label=target['name']+' ('+str(me.separation(pointing, direction))+')')
   
title('Pointing Elevation')
title('Elevation')
ylabel('Elevation (deg)');
xlabel('Time (h)');
l = legend(loc=5)
l.get_frame().set_alpha(0.5)
savefig('elevation.png')
show()
