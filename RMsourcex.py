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
#
# Version: 1.0
#
# RMsourceex.py will use a PyBDSM skymodel to look in a RM datacube for peaks above a certain
# noise treashold.
#
# Example: ./RMsourceex.py -c RMcube

import os, sys, re
import numpy as np
import pylab as pl

class Cube:
	"""
	A class to handle data cube
	"""
	def __init__(self, cube_data):
		self.cube_data = cube_data
		self.x_max = cube_data.shape[0]
		self.y_max = cube_data.shape[1]
		self.z_max = cube_data.shape[2]

	def pix2sky(self): pass

	def sky2pix(self): pass

	def getz(self,pix):
		"""
		Get the FD from a pixel
		Return:
		fd = the z vlues of the cube
		"""
		x = pix[0]
		y = pix[1]
		return self.cube_data[x][y]

	def findRMpeaks(self, pix, threshold):
		"""
	    	Look inside an RMcube at a certain pixel for peak above a threshold
	    	x = ra in pixel
		y = dec in pixel
		threshold = number of sigma above witch to look
		Return:
		detections = array of faraday depts pixels where detections are positive
		"""
		sigma = np.std(self.getz(pix))
		detections = []
		for i, phi in enumerate(self.getz(pix)):
   		 	if phi > threshold*sigma: detections.append(i)
 	   	return detections
		
class CubeFits(Cube):
	"""
	Class specific for fits files
	"""
	def __init__(self, cube_file):
		"""
		Extract the cube_data from the cube_file, set max and min values (in sky coordinates)
		"""
		import pyfits, pywcs
		# Put the cube in RA - DEC - RM order and save it
		Cube.__init__(self, np.transpose(pyfits.getdata(cube_file), (2, 1, 0)))
		self.wcs = pywcs.WCS(pyfits.getheader(cube_file))

		sky0 = self.pix2sky([0,0,0])
	    	skyN = self.pix2sky([self.x_max,self.y_max,self.z_max])
	    	self.ra_min = min(sky0[0],skyN[0])
		self.ra_max = max(sky0[0],skyN[0])
		self.ra_step = (self.ra_max-self.ra_min)/self.x_max
	    	self.dec_min = min(sky0[1],skyN[1])
	        self.dec_max = max(sky0[1],skyN[1])
		self.dec_step = (self.dec_max-self.dec_min)/self.y_max
		self.fd_min = min(sky0[2],skyN[2])
		self.fd_max = max(sky0[2],skyN[2])
		self.fd_step = (self.fd_max-self.fd_min)/self.z_max

	def pix2sky(self, pix):
		"""
		Convert pixel to sky coordinates
		pix = [ra,dec,fd] in pixel
		Return:
		sky = [ra,dec,fd] in sky coordinates
		"""
		return self.wcs.wcs_pix2sky([pix], 0)[0]

	def sky2pix(self, sky):
		"""
		Convert pixel to sky coordinates
		pix = [ra,dec,fd] in pixel
		Return:
		sky = [ra,dec,fd] in sky coordinates
		"""
		return self.wcs.wcs_sky2pix([sky], 0)[0]
		
	# not possible for RA and DEC separately since they are interconnected!
	def z2fd(self, z): return self.pix2sky([0,0,z])[2]
	def fd2z(self, fd): return round(self.sky2pix([0,0,fd])[2])
	def xy2radec(self, x, y): return self.pix2sky([x,y,0])[0:2]
	def radec2xy(self, ra, dec): return [round(x) for x in self.sky2pix([ra,dec,0])[0:2]]
	
class CubeMS(Cube):
	"""
	Class specific for MS files [TO BE IMPLEMENTED]
	"""
	def __init__(self, cube_file):
		sys.exit("ERROR: MS support still missing!")

def hmstora(rah,ram,ras):
	"""
	Convert RA in hours, minutes, seconds format to decimal
	degrees format.
	rah,ram,ras = RA values (h,m,s)
	Return:
	radegs = RA in decimal degrees
	"""
	hrs = (float(rah)+(float(ram)/60)+(float(ras)/3600.0)) % 24

	return 15*hrs

def dmstodec(decd,decm,decs):
	"""
	Convert Dec in degrees, minutes, seconds format to decimal
	degrees format.
	decd,decm,decs = list of Dec values (d,m,s)
	Return:
	decdegs = Dec in decimal degrees
	"""
	if decd < 0:
	    decm = -1*decm
	    decs = -1*decs

	decdegs = float(decd)+(float(decm)/60)+(float(decs)/3600.0)

	if abs(decdegs) > 90:
	    raise ValueError

	return decdegs

def createSkymodel(cube):
	"""
	Create a fake skymodel where every pixel is a source
	cube = the data cube object
	Return:
	skymodel_data = numpy array with all the pixels as sources
	"""
	dtypes = np.dtype({'names':['name','ra','dec'],'formats':['S100',np.float,np.float]})
	skymodel_data = np.ndarray((cube.x_max*cube.y_max),dtype=dtypes)
	for i in xrange(cube.x_max):
		sys.stdout.write("\rPopulating skymodel: %d%%" % ((i+1)*100./cube.x_max))
		sys.stdout.flush()
		for j in xrange(cube.y_max):
			skymodel_data['name'][j+i*cube.x_max] = str(i)+' - '+str(j)
			ra, dec  = cube.xy2radec(i,j)
			x, y = cube.radec2xy(ra,dec)
			skymodel_data['ra'][j+i*cube.x_max] = ra
			skymodel_data['dec'][j+i*cube.x_max] = dec
	return skymodel_data

def readPyBDSM(skymodel):
	"""
	Read PyBDSM skymodel file (BBS format)
	skymodel = the skymodel file generated by PyBDSM
	Return:
	skymodel_data = numpy array with all the sources in the skymodel
	"""
	names = []
	formats = []
	converters = {}
	f = open(skymodel, 'r')
	header = f.readline()
	header = re.sub(r'\s', '', header)
	header = re.sub(r'format=', '', header)
	headers = header.split(',')
	for i, h in enumerate(headers):
	    h = re.sub(r'=.*', '', h)
	    if h == "Name":
	        names.append('name')
	        formats.append('S100')
	    elif h == "Type":
	        names.append('type')
	        formats.append('S100')
	    elif h == "Ra":
	        names.append('ra')
	        formats.append(np.float)
	    	converters[i] = lambda x: hmstora(x.split(':')[0],x.split(':')[1],x.split(':')[2]) if x else 0.0
	    elif h == "Dec":
	        names.append('dec')
	        formats.append(np.float)
	   	converters[i] = lambda x: dmstodec(x.split('.')[0],x.split('.')[1],x.split('.')[2]) if x else 0.0
	    else:
	        names.append(h)
	        formats.append('S100')
	types = np.dtype({'names':names,'formats':formats})
	skymodel_data = np.loadtxt(skymodel, comments='format', unpack=True, dtype=types, delimiter=', ', converters=converters)
	return skymodel_data

def fixSkymodel(cube, skymodel_data):
	"""
	Remove sources present in the skymodel when outside the cube range
	cube = the RM cube data object
	skymodel_data = the skymodel
	Return:
	skymodel_data = an updated skymodel
	"""
	filt = []
	sou_in_field = []
	for i, source in enumerate(skymodel_data):
		ra = source['ra']
		dec = source['dec']
		# check if the source is inside the cube
		if ra > cube.ra_max or ra < cube.ra_min or dec > cube.dec_max or dec < cube.dec_min: filt.append(i)
	skymodel_data = np.delete(skymodel_data, filt)
	if (len(skymodel_data) == 0): sys.exit("No source found.")
	return skymodel_data

import optparse
opt = optparse.OptionParser(usage="%prog -c RMcube", version="%prog 1.0")
opt.add_option('-c', '--cube', help='RM cube to analyze in fots dormat (string)', type="string")
opt.add_option('-s', '--skymodel', help='Skymodel in BBS-format (string)', type="string", default='')
opt.add_option('-t', '--threshold', help='Number of sigma needed for a valid detection (double, default=5)', type="float", default=5.)
opt.add_option('-d', '--imagedir', help='Dir where to create plots (string, default="./")', type="string", default='./')
opt.add_option('-i', '--ignore', help='Comma separated list of Faraday depths to ignore for detection, x..y is interpreted as all the numbers between x and y (string, default='')', type="string", default='')
opt.add_option('-n', '--nospec', help='Do not produce spectra (bool, default=False)', action="store_true", default=False)
(options, args) = opt.parse_args()
cube_file = options.cube
skymodel = options.skymodel
threshold = options.threshold
imgdir = options.imagedir
ignore = options.ignore
nospec = options.nospec

# Options check
if cube_file == None:
	sys.exit("ERROR: RM cube needed")
if not os.access(imgdir, os.W_OK):
	sys.exit("ERROR: image directory not writable")

# Open RMcube (if it is a dir arrumes MSfile, otherwise a fits file)
print "Loading cube...",
sys.stdout.flush()
if os.path.isdir(cube_file) == True:
	cube = CubeMS(cube_file)
else:
	cube = CubeFits(cube_file)
print "done."
sys.stdout.flush()

# remove from faraday depths
ignore = []
if options.ignore != '':
	for i in options.ignore.split(','):
		if '..' in i:
			init, end = i.split('..')
			ignore = ignore + list(np.arange(float(init),float(end)+cube.fd_step,cube.fd_step))
		else:
			ignore = ignore + [float(i)]
	print "Ignoring: ", ignore

if skymodel == '':
	print "No skymodel, using the whole datacube."
	sys.stdout.flush()
	# create an artificial skymodel
	skymodel_data = createSkymodel(cube)
else:
	print "Loading skymodel...",
	sys.stdout.flush()
	skymodel_data = readPyBDSM(skymodel)
	# Remove sources outside the RM cube
	skymodel_data = fixSkymodel(cube, skymodel_data)
	print "done."
	sys.stdout.flush()

detections = dict((idx,[]) for idx in skymodel_data['name'])

# Look for peaks
detects_fd = []
detections_num = 0
for i, source in enumerate(skymodel_data):
	sys.stdout.write("\rSeeking peaks: %d%% (detections=%i)" % ((i+1)*100./len(skymodel_data), detections_num))
	sys.stdout.flush()
	detect_pix = cube.findRMpeaks(cube.radec2xy(source['ra'],source['dec']), threshold)
	detect_fd = [cube.z2fd(x) for x in detect_pix if cube.z2fd(x) not in ignore]
	detections_num += len(detect_fd)
	#print "DEBUG: ", detect_fd
	if detect_fd != []:
		detections[source['name']] = detect_fd
		detects_fd += list(detect_fd)
	        # Plot Faraday depth spectra of detected sources
		if nospec == False:
			fig=pl.figure()
			ax=fig.add_subplot(111)
			ax.plot(np.linspace(cube.fd_min,cube.fd_max,cube.z_max),cube.getz(cube.radec2xy(source['ra'],source['dec'])),'k-')
			ymin, ymax = ax.get_ylim()
			for dfd in detect_fd:
				ax.vlines(dfd, ymin, ymax,color='r', linestyles='solid')
			ax.plot(np.linspace(cube.fd_min,cube.fd_max,cube.z_max),cube.getz(cube.radec2xy(source['ra'],source['dec'])),'k-')
			x, y = cube.radec2xy(source['ra'],source['dec'])
			pl.xlabel('Faraday depth')
			pl.ylabel(r'F [Rad m$^{-2}$]')
			pl.title(source['name']+'-'+str(source['ra'])+' '+str(source['dec']))
			pl.savefig(imgdir+'/'+source['name']+'.png', format='png')
print "\n",

if detects_fd == []: sys.exit("No detections :(")

# Plot histogram of Faraday depths distribution
fig=pl.figure(figsize=(8, 8))
ax=fig.add_subplot(111)
ax.hist(detects_fd,bins=cube.z_max)
pl.xlabel('Faraday depth')
pl.ylabel('Number of positive detections')
pl.savefig(imgdir+'/hist.png', format='png')

# Plot spatial distribution of detections
fig=pl.figure(figsize=(8, 8))
ax=fig.add_subplot(111)
for source in skymodel_data:
	if detections[source['name']] != []: ax.plot(-source['ra'],source['dec'],'k.')
pl.xlabel('Ra')
pl.ylabel('Dec')
pl.savefig(imgdir+'/plane.png', format='png')

# Create ds9 region file
outfile=open(imgdir+'/detections.reg','w')
outfile.write('# Region file format: DS9 version 4.1\n')
outfile.write('# Filename: '+cube_file+' global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
outfile.write('fk5\n')
for source in skymodel_data:
	if detections[source['name']] != []:
		outfile.write('circle('+str(source['ra'])+','+str(source['dec'])+',10") # '+str(detections[source['name']])+'\n')
outfile.close()
