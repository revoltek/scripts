#!/usr/bin/env python
# Inject fake point sources in a residual mosaic and check which ones PyBDSF recovers.
#
# single mosaic:
#   survey_completeness.py HPxxxx-mosaicI.fits -o HPxxxx-mosaicI.injection.txt -r HPxxxx-mosaicI.rms.fits -e HPxxxx-mosaicI.mean.fits
# all the residual mosaics in a directory (rms and mean maps are taken from the catalogue directory):
#   survey_completeness.py --batch ~/storage/surveytgts/mosaics/healpix-res/ -r ~/storage/surveytgts/catalogues/ -e ~/storage/surveytgts/catalogues/
#
# -i is the total number of cycles wanted: a run on an existing output file adds only the
# missing cycles. The results are analysed by survey_completeness_analysis.py


from astropy.io import fits
import os, sys, subprocess
from astropy import wcs
import numpy as np
import bdsf as bdsm
import argparse
import numexpr as ne
import re

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

gfactor=2.0*np.sqrt(2.0*np.log(2.0))

parser = argparse.ArgumentParser(description='Completness testing')
parser.add_argument('file', metavar='FILE', nargs='*',
                   help='Root name of file to process (a residual image)')
parser.add_argument('-b','--batch',dest='batch',action='store',
                    default=None, metavar='DIR',
                    help='Auto-discover *mosaicI.fits files in DIR and process each one; -r and -e must then be the directories with the rms and mean maps')
parser.add_argument('-p','--plot',dest='plot',action='store_const',
                    const=True,default=False,
                    help='Produce diagnostic plots')
parser.add_argument('-i','--iterations',dest='iter',action='store',
                    type=int, default=5,
                    help='Total number of iterations wanted; the ones already in the output file are not repeated (default: 5)')
parser.add_argument('-t','--tolerance',dest='toler',action='store',
                    type=float, default=None,
                    help='Tolerance for source matches in pixels (default: use --tolerance-beam)')
parser.add_argument('-T','--tolerance-beam',dest='toler_beam',action='store',
                    type=float, default=1.0,
                    help='Tolerance for source matches in units of the beam FWHM (default: 1)')
parser.add_argument('--redo',dest='redo',action='store_const',
                    const=True,default=False,
                    help='Start from scratch even if the output file exists')
parser.add_argument('--seed',dest='seed',action='store',
                    type=int, default=None,
                    help='Seed of the random number generator (default: not set)')
parser.add_argument('-n','--number',dest='number',action='store',
                    type=int, default=200,
                    help='Number of fake sources to add')
parser.add_argument('-m','--minimum',dest='min',action='store',
                    type=float, default=3e-3,
                    help='Minimum flux to use')
parser.add_argument('-M','--maximum',dest='max',action='store',
                    type=float, default=10,
                    help='Maximum flux to use')
parser.add_argument('-I','--index',dest='index',action='store',
                    type=float, default=-0.6,
                    help='Power-law index of integrated number counts')
parser.add_argument('-s','--smear',dest='smear',action='store',
                    type=float, default=0,
                    help='Fraction of smearing (e.g. 0.2 for 20% reduction in peak flux)')
parser.add_argument('-o','--outfile',dest='outfile',action='store',
                    default=None,
                    help='File to log results to')
parser.add_argument('-r','--rmsmap',dest='rmsmap',action='store',
                    default=None,
                    help='BDSF RMS map (important to re-use the originally created one as the "empty" image has different rms properties)')
parser.add_argument('-e','--meanmap',dest='meanmap',action='store',
                    default=None,
                    help='BDSF mean map')
args = parser.parse_args()

if args.batch:
    import glob as _glob
    residuals = sorted(_glob.glob(os.path.join(args.batch, '*mosaicI.fits')))
    if not residuals:
        print('No *mosaicI.fits files found in %s.' % args.batch)
        sys.exit(1)
    if args.rmsmap is None or args.meanmap is None:
        parser.error('--batch requires -r and -e (directories with the rms and mean maps)')
    for residual in residuals:
        print('Batch: working on %s' % residual)
        basename = os.path.basename(residual).replace('.fits', '')
        rms     = os.path.join(args.rmsmap, basename + '.rms.fits')
        mean    = os.path.join(args.meanmap, basename + '.mean.fits')
        outfile = os.path.join(args.batch, basename + '.injection.txt')
        cmd = [sys.argv[0], residual,
               '-m', str(args.min), '-M', str(args.max),
               '-i', str(args.iter), '-n', str(args.number),
               '-I', str(args.index), '-T', str(args.toler_beam),
               '-o', outfile, '-r', rms, '-e', mean]
        if args.toler is not None:
            cmd += ['-t', str(args.toler)]
        if args.smear:
            cmd += ['-s', str(args.smear)]
        if args.redo:
            cmd += ['--redo']
        if args.seed is not None:
            cmd += ['--seed', str(args.seed)]
        print('Running command: %s' % ' '.join(cmd))
        subprocess.run(cmd, check=True)
    sys.exit(0)

if not args.file:
    parser.error('Provide at least one FILE or use --batch')

fn=args.file[0]
fitsfile=fn
fakefile=re.sub('.fits','.fake.fits',fitsfile)

if (args.plot):
    import matplotlib.pyplot as plt

if args.seed is not None:
    np.random.seed(args.seed)

# matching tolerance in pixels, from the beam unless given explicitly
with fits.open(fitsfile) as fp:
    fwhm_pix=fp[0].header['BMAJ']/abs(fp[0].header['CDELT2'])
toler=args.toler if args.toler is not None else args.toler_beam*fwhm_pix
print('Matching tolerance: %.2f pixels (beam FWHM: %.2f pixels)' % (toler, fwhm_pix))

# the first line of the output file records how the matches were made; cycles are
# appended to an existing file only if it was made in the same way
tag='# tolerance_pix=%.3f min=%g max=%g index=%g number=%i smear=%g' % (toler, args.min, args.max, args.index, args.number, args.smear)
first_cycle=0
if args.outfile:
    print('Logging to file', args.outfile)
    if os.path.exists(args.outfile) and not args.redo:
        with open(args.outfile) as old:
            lines=old.read().splitlines()
        if not lines or lines[0]!=tag:
            print('ERROR: %s was made with different parameters (or by an older version of this script); use --redo to overwrite it.' % args.outfile)
            sys.exit(1)
        first_cycle=1+max([int(l.split()[0]) for l in lines[2:]], default=-1)
        if first_cycle>=args.iter:
            print('%s already has %i cycles, nothing to do.' % (args.outfile, first_cycle))
            sys.exit(0)
        print('%s already has %i cycles, adding %i.' % (args.outfile, first_cycle, args.iter-first_cycle))
        outfile = open(args.outfile, 'a')
    else:
        outfile = open(args.outfile, 'w')
        outfile.write(tag+'\n')
        outfile.write('cycle idx ra dec fv rfv efv\n')

pnorm=args.min**args.index-args.max**args.index
sms=args.min**args.index

for c in range(first_cycle,args.iter):

    print('Doing iteration #',c)
    fp=fits.open(fitsfile)
    f=fp[0].data
    prhd=fp[0].header
    bmaj=prhd.get('BMAJ')
    bmin=prhd.get('BMIN')
    bpa=prhd.get('BPA')
    restfrq=prhd.get('RESTFRQ')
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
    ras=np.zeros(sources)
    decs=np.zeros(sources)
    for i in range(sources):
        while True:
            x=np.random.random_sample()*maxx
            y=np.random.random_sample()*maxy
            if not(np.isnan(f[int(y),int(x)])):
                break
            
        fl=(sms-pnorm*np.random.random())**(1.0/args.index)
        #print(i,x,y,fl)
        xp[i]=x
        yp[i]=y
        fv[i]=fl
        [[ra,dec]]=w.wcs_pix2world([[x,y]],0)
        ras[i] = ra
        decs[i] = dec

        if args.smear > 0:
            restore_gaussian(f,fl*(1-args.smear),x,y,bmaj/np.sqrt(1-args.smear),bmin/np.sqrt(1-args.smear),bpa,guard)
        else:
            restore_gaussian(f,fl,x,y,bmaj,bmin,bpa,guard)
        
    fp.writeto(fakefile, overwrite=True)
    fp.close()
    #restfrq=54e6 # should work this out from the FITS headers eventually
    if args.rmsmap and args.meanmap:
        fake_dir = os.path.dirname(os.path.abspath(fakefile))
        rmsmap_rel  = os.path.relpath(os.path.abspath(args.rmsmap),  fake_dir)
        meanmap_rel = os.path.relpath(os.path.abspath(args.meanmap), fake_dir)
        img=bdsm.process_image(fakefile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=50, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0,output_opts=True, output_all=True, atrous_do=False,atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None,frequency=restfrq, rmsmean_map_filename=[meanmap_rel,rmsmap_rel])
    else:
        img=bdsm.process_image(fakefile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=50, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0,output_opts=True, output_all=True, atrous_do=False,atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None,frequency=restfrq)

    found=[]
        
    for s in img.sources:
        (ra,dec)=s.posn_sky_centroid
        [[x,y]]=w.wcs_world2pix([[ra,dec]],0)
        #print x,y,s.total_flux
        d=np.sqrt((xp-x)**2.0+(yp-y)**2.0)
        i=np.argmin(d)
        fd=abs(fv[i]-s.total_flux)/s.total_fluxE
        if (d[i]<toler) and (fd<10):
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

    if args.outfile:
        for i in range(sources):
            outfile.write('%i %i %.6f %.6f %g %g %g\n' % (c,i, ras[i], decs[i], fv[i], rfv[i], efv[i]))
        outfile.flush() # a completed cycle is safe even if a later one is interrupted

if args.outfile:
    outfile.close()

