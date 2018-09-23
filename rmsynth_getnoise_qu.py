#!/usr/bin/env python

import pyfits
from numpy import *
import sys,os
import glob,argparse

def myfit(x,y, fn):
        # Find the width of a Gaussian distribution by computing the second moment of
        # the data (y) on a given axis (x)

        w = sqrt(abs(sum(x**2*y)/sum(y)))

        # Or something like this. Probably a Gauss plus a power-law with a lower cutoff is necessary
        #func = lambda x, a, b: np.exp(0.5*x**2/a**2) + x**(-1.*b)

        #[popt, pvar] = curve_fit(func, x, y)
        # Need a newer version of scipy for this to work...

        #a = popt[0]
        #b = popt[1]
        if fn!='':
          import pylab
          pylab.ioff()
          pylab.semilogy(x, y, '+', label="Data")
          #pylab.semilogy(x, np.exp(-0.5*x**2/a**2) + x**(-1.*b), label="Noise fit")
          pylab.semilogy(x, exp(-0.5*x**2/w**2), label="Noise fit")
          pylab.title('Normalized pixel distribution')
          pylab.ylabel('Rel. Num. of Pixels')
          pylab.xlabel('Pixel brightness (Jy/beam)')
          pylab.legend(loc=3)
          pylab.savefig(fn)
          #pylab.show()
          pylab.close()
        return w

def calcnoise(fitsfile):
  h = pyfits.open(fitsfile)
  d = squeeze(h[0].data)
  imsize = d.shape
  assert len(imsize)==2
  nx = imsize[-1]
  ny = imsize[-2]
  Id = d[ny/3:2*ny/3,nx/3:2*nx/3].flatten()
  if len(Id[isnan(Id)])==len(Id):
    return -1.
  else:
    rms = std(Id)
    mval = mean(Id)
    Id = Id[logical_and(Id<mval+3.*rms,Id>mval-3.*rms)]
    #print mval,rms,len(Id)

    #hrange = (-1,1)
    Ih = histogram(Id, bins=100)#, range=hrange) # 0 = values, 1 = left bin edges
    if max(Ih[0])==0.:
      return -1.
    Ix = Ih[1][:-1] + 0.5*(Ih[1][1] - Ih[1][0])
    Iv = Ih[0]/float(max(Ih[0]))
    Inoise = myfit(Ix, Iv, '')

    return Inoise

def main(args):
  qfitslist = sorted(glob.glob(args.qfitslist))
  ufitslist = sorted(glob.glob(args.ufitslist))
  #ifitslist = sorted(glob.glob(args.ifitslist)) # added by Shane

  old_stdout = sys.stdout
  log_file = open("QU_noise.txt","w")
  sys.stdout = log_file

  assert len(qfitslist)==len(ufitslist)
  qnoisevals = zeros(len(qfitslist))
  unoisevals = zeros(len(ufitslist))
  for i, fitsfile in enumerate(qfitslist):
    qnoisevals[i] = calcnoise(fitsfile)
  for i, fitsfile in enumerate(ufitslist):
    unoisevals[i] = calcnoise(fitsfile)
  qmeannoise = median(qnoisevals[abs(qnoisevals)<1.])
  qstdnoise = std(qnoisevals[abs(qnoisevals)<1.])
  #print 'Q median, std:',qmeannoise,qstdnoise
  print qmeannoise,qstdnoise
  umeannoise = median(unoisevals[abs(unoisevals)<1.])
  ustdnoise = std(unoisevals[abs(unoisevals)<1.])
  #print 'U median, std:',umeannoise,ustdnoise
  print umeannoise,ustdnoise
  qbadones = logical_or(qnoisevals>(qmeannoise+args.cliplev*qstdnoise),qnoisevals==-1.)
  ubadones = logical_or(unoisevals>(umeannoise+args.cliplev*ustdnoise),unoisevals==-1.)
  #print sum(asarray(qbadones,dtype=int)),'of',len(qfitslist),'are bad (Q)'
  #print sum(asarray(ubadones,dtype=int)),'of',len(ufitslist),'are bad (U)'
  totalbad = logical_or(qbadones,ubadones)
  #print sum(asarray(totalbad,dtype=int)),'of',len(qfitslist),'are bad in Q -or- U'
  print sum(asarray(totalbad,dtype=int))
  #print fitslist[asarray(badones,dtype=int)]
  if not args.delete: print 'Nothing will be deleted, but these are the files that would be with the -d option activated:'
  for i,f in enumerate(qfitslist):
    if totalbad[i]:
	#print qfitslist[i],ufitslist[i]
	if args.delete:
		os.remove(qfitslist[i])
		os.remove(ufitslist[i])
		#os.remove(ifitslist[i]) # added by Shane
  if not args.delete: print 'Nothing in the above list was deleted, use -d to take that action'

  sys.stdout = old_stdout
  log_file.close()

ap = argparse.ArgumentParser()
ap.add_argument('qfitslist',help='Wildcard list of Q fits files')
ap.add_argument('ufitslist',help='Wildcard list of U fits files')
#ap.add_argument('ifitslist',help='Wildcard list of I fits files')
ap.add_argument('--delete','-d',help='Delete bad channel maps? [default False]',default=False,action='store_true')
ap.add_argument('--cliplev','-c',help='Clip level in sigma, make this number lower to be more aggressive [default 5]',default=5.,type=float)
args=ap.parse_args()
main(args)

