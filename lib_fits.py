#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
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

from astropy.wcs import WCS as pywcs
from astropy.io import fits as pyfits
import numpy as np
import os, sys, logging, re

def flatten(filename, channel=0, freqaxis=0, hdu=0):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    f = pyfits.open(filename)

    naxis=f[hdu].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        #pass
        return f[hdu].header,f[hdu].data

    w = pywcs(f[hdu].header)
    wn = pywcs(naxis=2)

    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]

    header = wn.to_header()
    header["NAXIS"]=2
    header["NAXIS1"]=f[hdu].header['NAXIS1']
    header["NAXIS2"]=f[hdu].header['NAXIS2']
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[hdu].header.get(k)
        if r:
            header[k]=r

    dataslice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            dataslice.append(np.s_[:],)
        elif i==freqaxis:
            dataslice.append(channel)
        else:
            dataslice.append(0)

    # add freq
    header["FREQ"] = find_freq(f[hdu].header)

    # add beam if present
    try:
        header["BMAJ"]=f[hdu].header['BMAJ']
        header["BMIN"]=f[hdu].header['BMIN']
        header["BPA"]=f[hdu].header['BPA']
    except:
        pass

    # slice=(0,)*(naxis-2)+(np.s_[:],)*2
    return header, f[hdu].data[tuple(dataslice)]


def correct_beam_header(header):
    """ 
    Find the primary beam headers following AIPS convenction
    """
    if ('BMAJ' in header) and ('BMIN' in header) and ('PA' in header): return header
    elif 'HISTORY' in header:
        for hist in header['HISTORY']:
            if 'AIPS   CLEAN BMAJ' in hist:
                # remove every letter from the string
                bmaj, bmin, pa = re.sub(' +', ' ', re.sub('[A-Z ]*=','',hist)).strip().split(' ')
                header['BMAJ'] = float(bmaj)
                header['BMIN'] = float(bmin)
                header['BPA'] = float(pa)
    return header

def find_freq(header):
    """
    Find frequency value in most common places of a fits header
    """
    if not header.get('RESTFRQ') is None and not header.get('RESTFRQ') == 0:
        return header.get('RESTFRQ')
    elif not header.get('FREQ') is None and not header.get('FREQ') == 0:
        return header.get('FREQ')
    else:
        for i in range(5):
            type_s = header.get('CTYPE%i' % i)
            if type_s is not None and type_s[0:4] == 'FREQ':
                return header.get('CRVAL%i' % i)

    return None # no freq information found
 

class Image(object):

    def __init__(self, imagefile):
        """
        imagefile: name of the fits file
        """

        self.imagefile = imagefile
        header = pyfits.open(imagefile)[0].header
        header = correct_beam_header(header)

        try:
            beam = [header['BMAJ'], header['BMIN'], header['BPA']]
        except:
            logging.warning('%s: No beam information found.' % self.imagefile)
            sys.exit(1)
        logging.debug('%s: Beam: %.1f" %.1f" (pa %.1f deg)' % \
                (self.imagefile, beam[0]*3600., beam[1]*3600., beam[2]))

        self.freq = find_freq(header)
        if self.freq is None:
            logging.error('%s: No frequency information found.' % self.imagefile)
            sys.exit(1)
        logging.debug('%s: Frequency: %.0f MHz' % (self.imagefile, self.freq/1e6))

        self.noise = None
        self.img_hdr, self.img_data = flatten(self.imagefile)
        self.set_beam(beam)
        self.set_freq(self.freq)
        self.ra = self.img_hdr['CRVAL1']
        self.dec = self.img_hdr['CRVAL2']

    def write(self, filename=None):
        if filename is None:
            filename = self.imagefile
        pyfits.writeto(filename, self.img_data, self.img_hdr, overwrite=True)

    def set_beam(self, beam):
        self.img_hdr['BMAJ'] = beam[0]
        self.img_hdr['BMIN'] = beam[1]
        self.img_hdr['BPA'] = beam[2]

    def set_freq(self, freq):
        self.img_hdr['RESTFREQ'] = freq
        self.img_hdr['FREQ'] = freq

    def get_beam(self):
        return [self.img_hdr['BMAJ'], self.img_hdr['BMIN'], self.img_hdr['BPA']]

    def get_wcs(self):
        return pywcs(self.img_hdr)

    def apply_region(self, regionfile, blankvalue=np.nan, invert=False):
        """
        Blank inside mask
        invert: blank outside region
        """
        import pyregion
        if not os.path.exists(regionfile):
            logging.error('%s: Region file not found.' % regionfile)
            sys.exit(1)

        logging.debug('%s: Apply region %s' % (self.imagefile, regionfile))
        r = pyregion.open(regionfile)
        mask = r.get_mask(header=self.img_hdr, shape=self.img_data.shape)
        if invert: self.img_data[~mask] = blankvalue
        else: self.img_data[mask] = blankvalue


    def apply_mask(self, mask, blankvalue=np.nan, invert=False):
        """
        Blank inside mask
        invert: blank outside mask
        """
        logging.debug('%s: Apply mask' % self.imagefile)
        if invert: self.img_data[~mask] = blankvalue
        else: self.img_data[mask] = blankvalue


#    def calc_noise(self, niter=1000, eps=None, sigma=5):
#        """
#        Return the rms of all the pixels in an image
#        niter : robust rms estimation
#        eps : convergency criterion, if None is 1% of initial rms
#        """
#        if eps == None: eps = np.nanstd(self.img_data)*1e-3
#        data = self.img_data[ ~np.isnan(self.img_data) ] # remove nans
#        if len(data) == 0: return 0
#        oldrms = 1.
#        for i in range(niter):
#            rms = np.nanstd(data)
#            if np.abs(oldrms-rms)/rms < eps:
#                self.noise = rms
#                #print('%s: Noise: %.3f mJy/b' % (self.imagefile, self.noise*1e3))
#                logging.debug('%s: Noise: %.3f mJy/b' % (self.imagefile, self.noise*1e3))
#                return rms
#
#            data = data[np.abs(data)<sigma*rms]
#            oldrms = rms
#        raise Exception('Noise estimation failed to converge.')

    def calc_noise(self, sigma=7):
        """
        Return the rms of all the pixels in an image
        """
        from astropy.stats import median_absolute_deviation
        data = self.img_data[ ~np.isnan(self.img_data) ] # remove nans
        if len(data) == 0: return 0
        mad = median_absolute_deviation(data)
        mad_old = 0.
        while mad != mad_old:
             #print('MAD: %f uJy"' % (mad*1e6))
             mad_old = mad
             mad = median_absolute_deviation(data[np.where(np.abs(data) < sigma * mad)])
    
        rms = np.nanstd( data[np.where(np.abs(data) < sigma * mad)] )
        if np.isnan(rms):
            raise("cal_noise didn't converge")
        logging.debug('Noise: %.3f mJy/b' % (rms*1e3))
        self.noise = rms
        return rms

    def convolve(self, target_beam):
        """
        Convolve *to* this rsolution
        beam = [bmaj, bmin, bpa]
        """
        from lib_beamdeconv import deconvolve_ell, EllipticalGaussian2DKernel
        from astropy import convolution

        # if difference between beam is negligible <1%, skip - it mostly happens when beams are exactly the same
        beam = self.get_beam()
        if (np.abs((target_beam[0]/beam[0])-1) < 1e-2) and (np.abs((target_beam[1]/beam[1])-1) < 1e-2) and (np.abs(target_beam[2] - beam[2]) < 1):
            logging.debug('%s: do not convolve. Same beam.' % self.imagefile)
            return
        # first find beam to convolve with
        convolve_beam = deconvolve_ell(target_beam[0], target_beam[1], target_beam[2], beam[0], beam[1], beam[2])
        if convolve_beam[0] is None:
            logging.error('Cannot deconvolve this beam.')
            sys.exit(1)
        logging.debug('%s: Convolve beam: %.3f" %.3f" (pa %.1f deg)' \
                % (self.imagefile, convolve_beam[0]*3600, convolve_beam[1]*3600, convolve_beam[2]))
        # do convolution on data
        bmaj, bmin, bpa = convolve_beam
        assert abs(self.img_hdr['CDELT1']) == abs(self.img_hdr['CDELT2'])
        pixsize = abs(self.img_hdr['CDELT1'])
        fwhm2sigma = 1./np.sqrt(8.*np.log(2.))
        gauss_kern = EllipticalGaussian2DKernel((bmaj*fwhm2sigma)/pixsize, (bmin*fwhm2sigma)/pixsize, (90+bpa)*np.pi/180.) # bmaj and bmin are in pixels
        self.img_data = convolution.convolve(self.img_data, gauss_kern, boundary=None)
        self.img_data *= (target_beam[0]*target_beam[1])/(beam[0]*beam[1]) # since we are in Jt/b we need to renormalise
        self.set_beam(target_beam) # update beam

    def apply_shift(self, dra, ddec):
        """
        Shift header by dra/ddec
        dra, ddec in degree
        """
        # correct the dra shift for np.cos(DEC*np.pi/180.) -- only in the log as the reference val doesn't need it!
        logging.info('%s: Shift %.2f %.2f (arcsec)' % (self.imagefile, dra*3600*np.cos(self.dec*np.pi/180.), ddec*3600))
        dec = self.img_hdr['CRVAL2']
        self.img_hdr['CRVAL1'] += dra
        self.img_hdr['CRVAL2'] += ddec

