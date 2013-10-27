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

import numpy as np
import pylab as pl
import matplotlib, itertools
from coordinates_mode import *

def angsep(ra1deg, dec1deg, ra2deg, dec2deg):
    """Returns angular separation between two coordinates (all in degrees)"""
    import math

    ra1rad=ra1deg*math.pi/180.0
    dec1rad=dec1deg*math.pi/180.0
    ra2rad=ra2deg*math.pi/180.0
    dec2rad=dec2deg*math.pi/180.0

    # calculate scalar product for determination
    # of angular separation
    x=math.cos(ra1rad)*math.cos(dec1rad)*math.cos(ra2rad)*math.cos(dec2rad)
    y=math.sin(ra1rad)*math.cos(dec1rad)*math.sin(ra2rad)*math.cos(dec2rad)
    z=math.sin(dec1rad)*math.sin(dec2rad)

    if x+y+z >= 1: rad = 0
    else: rad=math.acos(x+y+z)

    # Angular separation
    deg=rad*180/math.pi
    return deg

def convert_sky_bbs_lsm(infilename,outfilename):

	types = np.dtype({'names':['Name', 'Type', 'Ra', 'Dec', 'I', 'Q', 'U', 'V', 'Maj', 'Min', 'PA', 'RefFreq', 'Spidx'],'formats':['S100','S100','S100','S100',np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,'S100']})
	data = np.loadtxt(infilename, comments='format', unpack=True, delimiter=', ', dtype=types)
	outfile=open(outfilename,'w')
	#name h m s d m s I Q U V spectral_index RM extent_X(rad) extent_Y(rad) pos_angle(rad) freq0
	outfile.write("## LSM file\n")
	outfile.write("### Name  | RA (hr,min,sec) | DEC (deg,min,sec) | I | Q | U | V | SI | RM | eX | eY | eP | freq0\n")
	for source in data:
		if source['Type'] == 'GAUSSIAN': source['Name']='G'+source['Name']
		elif source['Type'] == 'POINT': source['Name']='P'+source['Name']
		else:
			print "Unknown object type, ignoring "+source['Name']
			continue
		rahh, ramm, rass = source['Ra'].split(':')
		decdeg, decmm, decss = source['Dec'].split('.',2)
		# convert volume to peak
		i = source['I']/(1.133*source['Maj']*source['Min'])
		# BBS maj and min are FWHM in arcsec, in sagecal are radii in rad
		maj=((source['Maj']/3600.)*np.pi/180.)
		min=((source['Min']/3600.)*np.pi/180.)
		# in rad
		pa=((source['PA'])*np.pi/180.)
		outfile.write(source['Name']+' '+rahh+' '+ramm+' '+rass+' '+decdeg+' '+decmm+' '+decss+' '+str(i)+' '+str(source['Q'])+' '+str(source['U'])+' '+str(source['V'])+' ')
		outfile.write(source['Spidx'][1:-1].split(' ')[0]+' 0 '+str(maj)+' '+str(min)+' '+str(pa)+' '+str(source['RefFreq'])+'\n')
	outfile.close()

def create_cluster(infilename,outfilename,num_clusters,dist,do_cluster_png):
	types = np.dtype({'names':['Name', 'Rah', 'Ram', 'Ras', 'Decd', 'Decm', 'Decs', 'I', 'Q', 'U', 'V', 'SI', 'RM', 'eX', 'eY', 'eP', 'freq0'],'formats':['S100',np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float]})
	data = np.loadtxt(infilename, comments='#', unpack=True, delimiter=' ', dtype=types)
	outfile=open(outfilename,'w')

	fig=pl.figure()
	ax=fig.add_subplot(111)
	color=itertools.cycle(['r','b','g','c','m','y','k','#fb8e0d','#fb8eb9','#7c01f9'])

	# Select Q brightest sources Q=10)
	numQ = 0
	Q = []
	output = {}
	data.sort(order='I')
	for source in data[::-1]:
		if numQ == num_clusters: break
		ra = hmstora(source['Rah'],source['Ram'],source['Ras'])
		dec = dmstodec(source['Decd'],source['Decm'],source['Decs'])
		# ignore sources closer than $dist deg
		too_close = False
		for q in Q:
			if abs(angsep(ra,dec,q['Ra'],q['Dec'])) < dist: too_close = True
		if too_close: continue
		c = color.next()
		print "Adding centroid with flux: "+str(source['I'])+"(Ra="+str(ra)+" - Dec:"+str(dec)+"), color:"+c
		Q.append({'Name':source['Name'],'Ra':ra,'Dec':dec,'c':c})
		output[source['Name']]=[]
		numQ += 1

	# create clusters
	for source in data:
		# get angular distance from q and source
		s_ra = hmstora(source['Rah'],source['Ram'],source['Ras'])
		s_dec = dmstodec(source['Decd'],source['Decm'],source['Decs'])
		dist = 180
		for q in Q:
			newdist = abs(angsep(s_ra,s_dec,q['Ra'],q['Dec']))
			if newdist < dist:
				dist = newdist
				centroid = q['Name']
				centroid_ra = q['Ra']
				centroid_dec = q['Dec']
				c = q['c']
		output[centroid].append(source['Name'])
		# plot
		ax.plot(s_ra,s_dec,'o',color=c)
		ax.plot(centroid_ra,centroid_dec,'*',color=c, markersize=25)

	# print cluster file
	for q in Q:
		outfile.write(q['Name']+' ')
		for s in output[q['Name']]:
			outfile.write(s+' ')
		outfile.write('\n')
	outfile.close()
	ax.axis('equal')	
	pl.xlabel('RA (deg)')
	pl.ylabel('Dec (deg)')
	if do_cluster_png:
		print "Writing clusters.png"
		pl.savefig('clusters.png', format='png')
		
if __name__ == '__main__':
    import sys
    argc=len(sys.argv)
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <PyBDSM BBS catalog>')
    parser.add_option('--lsm', dest='file_lsm', help='Name of output lsm (sagecal) catalog [default: catalog.lsm]', metavar='VAL', default='catalog.lsm')
    parser.add_option('--clusters', dest='file_clusters', help='Name of the output cluster file [default catalogue.cluster]', metavar='VAL', default='catalog.clusters')
    parser.add_option('--numclusters', dest='num_clusters', help='Number of clusters [default: 10]', metavar='VAL', default=10)
    parser.add_option('--dist', dest='dist', help='Minimum distance from clusters brightest sources in deg [default: 0.3 deg]', metavar='VAL', default=0.3)
    parser.add_option('-c', action='store_true', dest='do_cluster_png', help='Print the png of the cluster disposition', default=False)
    (o, args) = parser.parse_args()

    if len(args) < 1: sys.exit("Missing BBS-format sky model.")
    convert_sky_bbs_lsm(args[0],o.file_lsm)
    create_cluster(o.file_lsm,o.file_clusters,int(o.num_clusters),float(o.dist),o.do_cluster_png)
    print "Done."
