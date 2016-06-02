#!/usr/bin/env python

import pyregion
import scipy.stats
from astropy.io import fits
from astropy import wcs
import numpy as np
import sys
import warnings
from linearfit import *

def flatten(f,channel=0,freqaxis=0):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        return f[0].header,f[0].data

    w = wcs.WCS(f[0].header)
    wn=wcs.WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    slice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            slice.append(np.s_[:],)
        elif i==freqaxis:
            slice.append(channel)
        else:
            slice.append(0)
        
# slice=(0,)*(naxis-2)+(np.s_[:],)*2
    return header,f[0].data[slice]

class RadioError(Exception):
    """Base class for exceptions in this module."""
    pass

class radiomap:
    """ Process a fits file as though it were a radio map, calculating beam areas etc """
    def __init__(self, filename, verbose=False):
        self.filename = filename
        self.fitsfile=fits.open(filename)
        # Catch warnings to avoid datfix errors
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gfactor=2.0*np.sqrt(2.0*np.log(2.0))
            self.f=self.fitsfile[0]
            self.prhd=self.fitsfile[0].header

            # Get units and resolution
            self.units=self.prhd.get('BUNIT')
            if self.units is None:
                self.units=self.prhd.get('UNIT')
            if self.units!='JY/BEAM' and self.units!='Jy/beam':
                print 'Warning: units are',self.units,'but code expects JY/BEAM'
            self.bmaj=self.prhd.get('BMAJ')
            self.bmin=self.prhd.get('BMIN')
            if self.bmaj is None:
                # Try RESOL1 and RESOL2
                self.bmaj=self.prhd.get('RESOL1')
                self.bmin=self.prhd.get('RESOL2')
            if self.bmaj is None:
                if verbose:
                    print 'Can\'t find BMAJ in headers, checking history'
                try:
                    history=self.prhd['HISTORY']
                except KeyError:
                    history=None
                if history is not None:
                    for line in history:
                        if 'HISTORY' in line:
                            continue # stops it finding nested history
                        if 'CLEAN BMAJ' in line:
                            bits=line.split()
                            self.bmaj=float(bits[3])
                            self.bmin=float(bits[5])
                                
            if self.bmaj is None:
                raise RadioError('No beam information found')

            w=wcs.WCS(self.prhd)
            cd1=-w.wcs.cdelt[0]
            cd2=w.wcs.cdelt[1]
            if ((cd1-cd2)/cd1)>1.0001 and ((self.bmaj-self.bmin)/self.bmin)>1.0001:
                raise RadioError('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

            self.bmaj/=cd1
            self.bmin/=cd2
            if verbose:
                print 'beam is',self.bmaj,'by',self.bmin,'pixels'

            self.area=2.0*np.pi*(self.bmaj*self.bmin)/(gfactor*gfactor)
            if verbose:
                print 'beam area is',self.area,'pixels'

            # Remove any PC... keywords we may have, they confuse the pyregion WCS
            for i in range(1,5):
                for j in range(1,5):
                    self.quiet_remove('PC0%i_0%i' % (i,j))
                
            # Now check what sort of a map we have
            naxis=len(self.fitsfile[0].data.shape)
            if verbose: print 'We have',naxis,'axes'
            self.cube=False
            if naxis<2 or naxis>4:
                raise RadioError('Too many or too few axes to proceed (%i)' % naxis)
            if naxis>2:
                # a cube, what sort?
                frequency=0
                self.cube=True
                freqaxis=-1
                stokesaxis=-1
                for i in range(3,naxis+1):
                    ctype=self.prhd.get('CTYPE%i' % i)
                    if 'FREQ' in ctype:
                        freqaxis=i
                    elif 'STOKES' in ctype:
                        stokesaxis=i
                    elif 'VOPT' in ctype:
                        pass
                    else:
                        print 'Warning: unknown CTYPE %i = %s' % (i,ctype)
                if verbose:
                    print 'This is a cube with freq axis %i and Stokes axis %i' % (freqaxis, stokesaxis)
                if stokesaxis>0:
                    nstokes=self.prhd.get('NAXIS%i' % stokesaxis)
                    if nstokes>1:
                        raise RadioError('Multiple Stokes parameters present, not handled')
                if freqaxis>0:
                    nchans=self.prhd.get('NAXIS%i' % freqaxis)
                    if verbose:
                        print 'There are %i channel(s)' % nchans
                    self.nchans=nchans
            else:
                self.nchans=1
                    

            # Various possibilities for the frequency. It's possible
            # that a bad (zero) value will be present, so keep
            # checking if one is found.

            if not(self.cube) or freqaxis<0:
                # frequency, if present, must be in another keyword
                frequency=self.prhd.get('RESTFRQ')
                if frequency is None or frequency==0:
                    frequency=self.prhd.get('RESTFREQ')
                if frequency is None or frequency==0:
                    frequency=self.prhd.get('FREQ')
                if frequency is None or frequency==0:
                    # It seems some maps present with a FREQ ctype
                    # even if they don't have the appropriate axes!
                    # The mind boggles.
                    for i in range(5):
                        type_s=self.prhd.get('CTYPE%i' % i)
                        if type_s is not None and type_s[0:4]=='FREQ':
                            frequency=self.prhd.get('CRVAL%i' % i)
                self.frq=[frequency]
                # now if there _are_ extra headers, get rid of them so pyregion WCS can work
                for i in range(3,5):
                    for k in ['CTYPE','CRVAL','CDELT','CRPIX','CROTA','CUNIT']:
                        self.quiet_remove(k+'%i' %i)
                self.headers=[self.prhd]
                self.d=[self.fitsfile[0].data]
            else:
                # if this is a cube, frequency/ies should be in freq header
                basefreq=self.prhd.get('CRVAL%i' % freqaxis)
                deltafreq=self.prhd.get('CDELT%i' % freqaxis)
                self.frq=[basefreq+deltafreq*i for i in range(nchans)]
                self.d=[]
                self.headers=[]
                for i in range(nchans):
                    header,data=flatten(self.fitsfile,freqaxis=freqaxis,channel=i)
                    self.d.append(data)
                    self.headers.append(header)
            for i,f in enumerate(self.frq):
                if f is None:
                    print('Warning, can\'t get frequency %i -- set to zero' % i)
                    self.frq[i]=0
            if verbose:
                print 'Frequencies are',self.frq,'Hz'
            #self.fitsfile.close()

    def quiet_remove(self,keyname):
        if self.prhd.get(keyname,None) is not None:
            self.prhd.remove(keyname)

                
#            self.fhead,self.d=flatten(fitsfile)


class applyregion:
    """ apply a region from pyregion to a radiomap """
    def __init__(self,rm,region,offsource=None,mask=None,wht=None,robustrms=3):
        """
        provides:
        rms -- the rms in the aperture
        robustrms -- the rms for pixels below robustrms * the normal rms (it should cut sources)
        flux -- the flux of the aperture
        mean -- the mean in the apertur (if wht then is weighted)
        error -- error on the flux given the rms in offsource
        """
        self.rms=[]
        self.flux=[]
        self.error=[]
        self.mean=[]
        self.robustrms=[]
        self.mean_error=[]

        for i,d in enumerate(rm.d):
            mask_r=region.get_mask(hdu=rm.f,shape=np.shape(d))
            pixels=np.sum(mask_r)
            if mask is not None:
                # save mask fits as debug
                #rm.fitsfile[0].data = np.logical_and(mask[i],mask_r).astype(np.float)
                #rm.fitsfile.writeto('debugmask.fits')
                #sys.exit(1)
                data = np.extract(np.logical_and(mask[i],mask_r),d)
            else:
                data = np.extract(mask_r,d)

            self.rms.append(scipy.stats.nanstd(data))
            self.robustrms.append(scipy.stats.nanstd(data[np.where(data < robustrms * self.rms[-1])]))
            self.flux.append(data[np.logical_not(np.isnan(data))].sum()/rm.area)

            if wht:
                mask=region.get_mask(hdu=wht.f,shape=np.shape(wht.d[0]))
                pixels=np.sum(mask)
                data_error=np.extract(mask,wht.d)
                #print data[np.logical_not(np.isnan(data))], data_error[np.logical_not(np.isnan(data))]
                # https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html
                self.mean.append( np.average( data[np.logical_not(np.isnan(data))], weights=1./data_error[np.logical_not(np.isnan(data_error))]**2 ) )
                self.mean_error.append( np.sqrt(1./np.sum(1./data_error[np.logical_not(np.isnan(data_error))]**2)) )
            else:
                self.mean.append(scipy.stats.nanmean(data))

            # calc noise
            if offsource is not None:
                self.error.append(offsource[i]*np.sqrt(pixels/rm.area))

def printflux(rms,region,bgr=None,wht=None,fluxerr=None,minrms=0,nsigma=0,label=''):
    """
    rms -- list of radio map objects
    region -- region to work on
    bgr -- background region for noise/subtraction (only flux and spidx)
    wht -- radiomap of the weights in sigma (only mean)
    fluxerr -- percentage of flux error for spidx maps (only spidx)
    """
    for rm in rms:
        if bgr:
            bg_ir = pyregion.open(bgr).as_imagecoord(rm.headers[0])
            bg = applyregion(rm,bg_ir)
            noise = bg.rms
        else:
            noise = None

        fg=applyregion(rm,region,offsource=noise)

        for i in range(rm.nchans):
            freq=rm.frq[i]
        
            if noise is not None:
                print rm.filename,label,'%8.4g %10.6g %10.6g' % (freq,fg.flux[i],fg.error[i])
            else:
                print rm.filename,label,'%8.4g %10.6g' % (freq,fg.flux[i])

def printmean(rms,region,bgr=None,wht=None,fluxerr=None,nsigma=0,label=''):
    for rm in rms:
        fg=applyregion(rm,region,wht=wht)

        for i in range(rm.nchans):
            freq=rm.frq[i]
            if wht:
                print rm.filename,label,'%8.4g %10.6g %10.6g' % (freq,fg.mean[i],fg.mean_error[i])
            else:
                print rm.filename,label,'%8.4g %10.6g' % (freq,fg.mean[i])

def printspidx(rms,region,bgr=None,wht=None,fluxerr=0,nsigma=0,label=''):
    freqs = []
    fluxes = []
    errors = []
    noises = []

    mask = (np.zeros_like(rms[0].d) == 0)
    for rm in rms:
        if bgr:
            bg_ir=pyregion.open(bgr).as_imagecoord(rm.headers[0])
            bg=applyregion(rm,bg_ir)
            noises.append(bg.rms)
            # likely brakes with channelled images
            if nsigma > 0: mask = np.logical_and(mask, np.array(rm.d) > (np.array(bg.rms)*nsigma) )
        else:
            noises.append(None)

    for i, rm in enumerate(rms):
        fg=applyregion(rm,region,offsource=noises[i],mask=mask)
        freqs += rm.frq
        fluxes += fg.flux
        if bgr: errors += fg.error

    # lin reg
    if bgr:
        yerr = 0.434*np.sqrt(np.array(errors)**2+(np.array(fluxerr)*np.array(fluxes)/100)**2)/np.array(fluxes)
    else:
        yerr = None

    (a, b, sa, sb) = linear_fit_bootstrap(x=np.log10(freqs), y=np.log10(fluxes), yerr=yerr)
    print label,'%8.4g %8.4g' % (a, sa)


def radioflux(files,fgr,bgr=None,individual=False,action='Flux',wht=None,fluxerr=0,nsigma=0,verbose=False):
    """Determine the flux in a region file for a set of files. This is the
    default action for the code called on the command line, but
    may be useful to other code as well.

    Keyword arguments:
    files -- list of files (mandatory)
    fdr -- foreground region name (mandatory)
    bgr -- background region name (optional)
    individual -- separate region into individual sub-regions
    action -- what to do once fluxes are measured: allows a user-defined action
              which must be a drop-in replacement for printflux
    wht -- weight file in sigmas
    fluxerr -- flux error in % for spidxmap
    """
    action = {'flux':printflux, 'mean':printmean, 'spidx':printspidx}[action]

    if wht:
        wht = radiomap(wht,verbose=verbose)

    rms = [] # radio maps
    for filename in files:
        rms.append(radiomap(filename,verbose=verbose))

    fg_ir=pyregion.open(fgr).as_imagecoord(rms[0].headers[0])

    if individual:
        for n,reg in enumerate(fg_ir):
            fg=pyregion.ShapeList([reg])
            r=action(rms,fg,bgr=bgr,wht=wht,fluxerr=fluxerr,nsigma=nsigma,label=n+1)
    else:
        r=action(rms,fg_ir,bgr=bgr,wht=wht,fluxerr=fluxerr,nsigma=nsigma)

    return r
        
if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Measure fluxes from FITS files.')
    parser.add_argument('files', metavar='FILE', nargs='+',
                        help='FITS files to process')
    parser.add_argument('-f','--foreground', dest='fgr', action='store',default='ds9.reg',help='Foreground region file to use')
    parser.add_argument('-b','--background', dest='bgr', action='store',default='',help='Background region file to use')
    parser.add_argument('-i','--individual', dest='indiv', action='store_true',default=False,help='Break composite region file into individual regions')
    parser.add_argument('-e','--fluxerr', dest='fluxerr', action='store',default=0, type=float, help='Flux error in % for spidx maps')
    parser.add_argument('-s','--sigma', dest='nsigma', action='store',default=0, type=float, help='Try to cut all the images above a certain sigma. Only pixel over that sigma in ALL the images are considered. Valid only for spidx.')
    parser.add_argument('-a','--action', dest='action', action='store',default='flux',help='Action to perform: flux, mean, spidx')
    parser.add_argument('-w','--weights', dest='weights', action='store',default=None,help='If action=mean then weight the mean with the values in this "sigma" map [e.g. fitsfile=spidx.fits, -w spidx-rms.fits gives a flux-weighted spidx]')
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true',default=False,help='Be verbose')

    args = parser.parse_args()

    radioflux(args.files,args.fgr,args.bgr,args.indiv,args.action,args.weights,args.fluxerr,args.nsigma,verbose=args.verbose)
