#!/usr/bin/python
# -*- coding: utf-8 -*-



#
# Written by Bas van der Tol (vdtol@strw.leidenuniv.nl), March 2011.
#

from pylab import *
import pyrap.quanta as qa
import pyrap.tables as pt
import pyrap.measures as pm
import sys
import numpy

targets = [ {'name' : 'CasA', 'ra' : 6.123487680622104,  'dec' : 1.0265153995604648},
            {'name' : 'CygA', 'ra' : 5.233686575770755,  'dec' : 0.7109409582180791},
            {'name' : 'TauA', 'ra' : 1.4596748493730913, 'dec' : 0.38422502335921294},
            {'name' : 'HerA', 'ra' : 4.4119087330382163, 'dec' : 0.087135562905816893},
            {'name' : 'VirA', 'ra' : 3.276086511413598,  'dec' : 0.21626589533567378},
            {'name' : 'SUN'},
            {'name' : 'JUPITER'}]


if len(sys.argv) == 2:
   msname = sys.argv[1]
else:
   print "Usage"
   print "   plot_Ateam_elevation.py <msname>"
   

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
print position
ant_table.close()


# Get the first pointing of the first antenna
field_table = pt.table(msname + '/FIELD')
field_no = 0
direction = field_table.getcol('PHASE_DIR')
ra = direction[ ant_no, field_no, 0 ]
dec = direction[ ant_no, field_no, 1 ]
targets.insert(0, {'name' : 'Pointing', 'ra' : ra, 'dec' : dec})
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

separations = []

for target in targets:
   
   t = qa.quantity(time[0], 's')
   t1 = me.epoch('utc', t)
   me.doframe(t1)

   if 'ra' in target.keys():
      ra_qa  = qa.quantity( target['ra'], 'rad' )
      dec_qa = qa.quantity( target['dec'], 'rad' )
      direction =  me.direction('j2000', ra_qa, dec_qa)
   else :
      direction =  me.direction(target['name'])
      
   separations.append(me.separation(pointing, direction))
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

   plot(time1, el)
   
title('Pointing Elevation')
title('Elevation')
ylabel('Elevation (deg)');
#ylim(0,90)
xlabel('Time (h)');
#savefig('elevation.eps')
legend( [ target['name'] + '(' + separation.to_string() + ')' for target, separation in zip(targets, separations) ])
show()
