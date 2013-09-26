#!/usr/bin/env python
## fixed select typo on line 49  8/8/12 PC

import pyrap.images as pim
import pyrap.tables as pt
from numpy import *
import os, sys
import optparse


def main(opts,args):
	robust = opts.robust
	maxbl = opts.maxbl
	fov = opts.fov
	skymodel = opts.skymodel
	output = opts.output
	verbose = opts.verbose
	datacol = opts.datacol
	msname = args[0]
	if msname[0:2] == './':
		msname = msname[2:len(msname)]

	imname = "%s_rob%+.1f_maxbl%d_fov%d.img"%(msname.split('.')[0],robust,maxbl,fov)
	mask = "%s_rob%+.1f_maxbl%d_fov%d.mask"%(msname.split('.')[0],robust,maxbl,fov)
	
	ft = pt.table(msname+'/SPECTRAL_WINDOW')
	freq = ft.getcell('REF_FREQUENCY',0)
	cell = 3.e8/freq/maxbl*206265./3.
	if verbose:
		print "Frequency (MHz) =",freq/1.e6
		print "Cell size (arcsec) =",cell

	fovarcsec = fov*3600.
	npixrough = fovarcsec/cell
	npix = 2**(int(log(npixrough)/log(2.))+1)
	if verbose:
		print "NPIX (rough, rounded) =",npixrough,npix
	if npix > 2048:
		print 'WARNING !!!!! Large number of pixels ....'
	
	sqrbaseline = maxbl**2

	wplanes = min(int(fov*pi/180.*maxbl/(3.e8/freq)/2.)*2+1,257)
	if verbose:
		print "Number of wplanes =",wplanes

	os.system('rm -rf tempABCDxxxx*')
	print 'Determining noise, please wait...'
	os.system('awimager cellsize=%darcsec data=CORRECTED_DATA ms=%s nchan=10 niter=0 npix=256 operation=image robust=%f select="sumsqr(UVW[:2])<%.1e" stokes=IQUV wmax=%f wprojplanes=%d image=tempABCDxxxx > /dev/null'%(cell,msname,robust,sqrbaseline,maxbl,wplanes))
	image = pim.image('tempABCDxxxx')
	data = image.getdata()
	noise = std(data[0,3,100:150,100:150])
	os.system('rm -rf tempABCDxxxx*')
	threshold = 10.*noise

	os.system('rm -rf %s'%mask)
	if skymodel != 'none':
		print 'Making mask, please wait..'
		os.system('awimager cellsize=%darcsec data=CORRECTED_DATA ms=%s nchan=10 niter=0 npix=%d operation=empty robust=%f select="sumsqr(UVW[:2])<%.1e" stokes=IQUV wmax=%f wprojplanes=%d image=%s > /dev/null'%(cell,msname,npix,robust,sqrbaseline,maxbl,wplanes,mask))
		os.system('rm -rf temp.db')
		os.system('makesourcedb in=%s out=temp.db format=Name,Type,Ra,Dec,I,Q,U,V,ReferenceFrequency=\\\"60e6\\\",SpectralIndex=\\\"[0.0]\\\",MajorAxis,MinorAxis,Orientation > /dev/null'%skymodel)
		os.system('~csobey/MSSS/msss_mask.py %s temp.db > /dev/null'%(mask))
		os.system('rm -rf temp.db')

	parameters="""
cellsize=%darcsec
cyclefactor=1.
data=%s
mask=%s
ms=%s
nchan=10
niter=2500
npix=%d
operation=csclean
robust=%f
select="sumsqr(UVW[:2])<%.1e"
stokes=IQUV
threshold=%fJy
wmax=%f
wprojplanes=%d
image=%s
"""%(cell,datacol,mask,msname,npix,robust,sqrbaseline,threshold,maxbl,wplanes,imname)

	outfile = open(output,'w')
	print >>outfile,parameters
	outfile.close()

	print '\n\n\n Done -- now run "time awimager %s"'%output

usage = "%prog [options] msfile"
opt = optparse.OptionParser(usage)
opt.add_option('-v','--verbose',help='Verbose mode',default=False,action='store_true')
opt.add_option('-r','--robust',help='Robust [0]',default=0,type='float')
opt.add_option('-b','--maxbl',help='Max baseline (meters) [3000]',default=3000,type='float')
opt.add_option('-f','--fov',help='Field of view (degrees) [10]',default=10,type='float')
opt.add_option('-s','--skymodel',help='Sky model name [sky.model], use \'none\' to avoid using a mask for deconvolution',default='sky.model',type='string')
opt.add_option('-o','--output',help='Output parset name [aw.parset]',default='aw.parset',type='string')
opt.add_option('-d','--datacol',help='Data column to image [CORRECTED_DATA]',default='CORRECTED_DATA',type='string')
opts,args = opt.parse_args()
if len(args) != 1: 
	print 'Need required argument msname'
else:
	main(opts,args)

