#!/usr/bin/env python
# ./fakesources-hetdex.py mosaic.fits_resid.fits -m 1e-3 --logfile=test-blanked-results.txt -i 50 -n 6000 -p -I -1.6

from astropy.io import fits
import sys
from astropy import wcs
import numpy as np
import bdsf as bdsm
import argparse
import numexpr as ne
import re

import collections
import functools

class memoized(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)

@memoized
def gaussian(xsize,ysize,x0,y0,sx,sy,pa):
    X,Y = np.meshgrid(np.arange(0,xsize,1.0), np.arange(0,ysize,1.0))
    pa*=np.pi/180.0
    a=0.5*((np.cos(pa)/sx)**2.0+(np.sin(pa)/sy)**2.0)
    b=0.25*((-np.sin(2*pa)/sx**2.0)+(np.sin(2*pa)/sy**2.0))
    c=0.5*((np.sin(pa)/sx)**2.0+(np.cos(pa)/sy)**2.0)
    
    return ne.evaluate('exp(-(a*(X-x0)**2.0+2*b*(X-x0)*(Y-y0)+c*(Y-y0)**2.0))')

def restore_gaussian(image,norm,x,y,bmaj,bmin,bpa,guard,verbose=False):
    # deal with real image co-ords
    yd,xd=np.shape(image)
    if x is int:
        xp=x
        x=x+0.5
        yp=y
        y=y+0.5
    else:
        xp=int(np.trunc(x))
        yp=int(np.trunc(y))
    if verbose:
       print(x,y,xp,yp,norm)
    if xp<0 or yp<0 or xp>xd-1 or yp>yd-1:
       raise Exception('position out of range')
    xmin=xp-guard
    xmax=xp+guard
    ymin=yp-guard
    ymax=yp+guard
    if xmin<0:
        xmin=0
    if ymin<0:
        ymin=0
    if xmax>=xd:
        xmax=xd-1
    if ymax>=yd:
        ymax=yd-1
    x0=x-xmin
    y0=y-ymin
    image[ymin:ymax,xmin:xmax]+=norm*gaussian(xmax-xmin,ymax-ymin,x0,y0,bmaj,bmin,bpa)


# Add some fake sources to a FITS image and then see if they can be
# detected with pybdsm

# HETDEX version: use same arguments as Tim did

gfactor=2.0*np.sqrt(2.0*np.log(2.0))

parser = argparse.ArgumentParser(description='Completness testing')
parser.add_argument('file', metavar='FILE', nargs='+',
                   help='Root name of file to process')
parser.add_argument('-p','--plot',dest='plot',action='store_const',
                    const=True,default=False,
                    help='Produce diagnostic plots')
parser.add_argument('-i','--iterations',dest='iter',action='store',
                    type=int, default=1,
                    help='Number of iterations to go through')
parser.add_argument('-t','--tolerance',dest='toler',action='store',
                    type=float, default=25.0,
                    help='Tolerance for source matches in pixels')
parser.add_argument('-n','--number',dest='number',action='store',
                    type=int, default=100,
                    help='Number of fake sources to add')
parser.add_argument('-m','--minimum',dest='min',action='store',
                    type=float, default=3e-2,
                    help='Minimum flux to use')
parser.add_argument('-M','--maximum',dest='max',action='store',
                    type=float, default=10,
                    help='Minimum flux to use')
parser.add_argument('-I','--index',dest='index',action='store',
                    type=float, default=-0.6,
                    help='Power-law index of integrated number counts')
parser.add_argument('-s','--smear',dest='smear',action='store',
                    type=float, default=0,
                    help='Fraction of smearing (e.g. 0.2 for 20% reduction in peak flux)')
parser.add_argument('-l','--logfile',dest='logfile',action='store',
                    default=None,
                    help='File to log results to')
args = parser.parse_args()

fn=args.file[0]
fitsfile=fn
rmsmap='mosaic.fits_rms.fits'
meanmap='mosaic.fits_mean.fits'
fakefile=re.sub('.fits','.fake.fits',fitsfile)

if (args.plot):
    import matplotlib.pyplot as plt

if (args.logfile!=None):
    print('Logging to file',args.logfile)
    outfile=open(args.logfile,'w')

pnorm=args.min**args.index-args.max**args.index
sms=args.min**args.index

for c in range(0,args.iter):

    print('Doing iteration #',c)
    fp=fits.open(fitsfile)
    f=fp[0].data[0,0]
    prhd=fp[0].header
    bmaj=prhd.get('BMAJ')
    bmin=prhd.get('BMIN')
    bpa=prhd.get('BPA')
    w=wcs.WCS(prhd)
    cd1=-w.wcs.cdelt[0]
    cd2=w.wcs.cdelt[1]
    if abs(cd1-cd2) > 1E-14:
        print(cd1,cd2,cd1-cd2)
        raise Exception('Pixels are not square')
    (maxy,maxx)=f.shape

    print('File',fn)
    print('BMAJ',bmaj,'BMIN',bmin,'BPA',bpa)
    print('Pixel size',cd1)
    print('maxx',maxx,'maxy',maxy)
    bmaj/=cd1
    bmin/=cd2
    bmaj/=gfactor
    bmin/=gfactor
    print('Gaussian axes in pixels',bmaj,bmin)
    guard=int(bmaj*10)
    #Now add the sources

    sources=args.number
    xp=np.zeros(sources)
    yp=np.zeros(sources)
    fv=np.zeros(sources)
    rfv=np.zeros(sources)
    efv=np.zeros(sources)
    for i in range(0,sources):
        while True:
            x=np.random.random_sample()*maxx
            y=np.random.random_sample()*maxy
            if not(np.isnan(f[int(y),int(x)])):
                break
            
        fl=(sms-pnorm*np.random.random())**(1.0/args.index)
        print(i,x,y,fl)
        xp[i]=x
        yp[i]=y
        fv[i]=fl

        if args.smear > 0:
            restore_gaussian(f,fl*(1-args.smear),x,y,bmaj/np.sqrt(1-args.smear),bmin/np.sqrt(1-args.smear),bpa,guard)
        else:
            restore_gaussian(f,fl,x,y,bmaj,bmin,bpa,guard)
        
    fp.writeto(fakefile,clobber=True)
    fp.close()
    restfrq=54e6 # should work this out from the FITS headers eventually
    img=bdsm.process_image(fakefile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(160,50), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0,output_opts=True, output_all=True, atrous_do=True,atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None,frequency=restfrq, rmsmean_map_filename=[meanmap,rmsmap])

    found=[]
        
    for s in img.sources:
        (ra,dec)=s.posn_sky_centroid
        [[x,y,d,q]]=w.wcs_world2pix([[ra,dec,0,0]],0)
        #        print x,y,s.total_flux
        d=np.sqrt((xp-x)**2.0+(yp-y)**2.0)
        i=np.argmin(d)
        fd=abs(fv[i]-s.total_flux)/s.total_fluxE
        if (d[i]<args.toler) and (fd<10):
            rfv[i]=s.total_flux
            efv[i]=s.total_fluxE
            print(x,y,s.total_flux,i,d[i],fv[i])
            found.append(i)

    print('Total sources found',len(found))
    if args.plot:
        plt.hist(np.log10(fv),20,range=(np.log10(args.min),np.log10(args.max)))
        plt.hist(np.log10(fv[found]),20,range=(np.log10(args.min),np.log10(args.max)))
        plt.xlabel('Log10(Flux)')
        plt.show()

        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Real flux (Jy)')
        plt.ylabel('Recovered flux (Jy)')
        plt.errorbar(fv[found],rfv[found],yerr=efv[found],fmt='o')
        plt.plot([fv[found].min()*0.9,fv[found].max()*1.1],[fv[found].min()*0.9,fv[found].max()*1.1])
        plt.show()

    if (args.logfile!=None):
        for i in range(0,sources):
            outfile.write('%i %i %g %g %g\n' % (c,i,fv[i],rfv[i],efv[i]))

