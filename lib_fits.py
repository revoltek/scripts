from astropy.wcs import WCS as pywcs
from astropy.io import fits as pyfits
import numpy as np
import os, sys, logging

def flatten(filename, channel=0, freqaxis=0):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    f = pyfits.open(filename)

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        return f[0].header,f[0].data

    w = pywcs(f[0].header)
    wn = pywcs(naxis=2)

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


def correct_beam_header(header):
    """ 
    Find the primary beam headers following AIPS convenction
    """
    import pyfits, re
    if ('BMAJ' in header) and ('BMIN' in header) and ('PA' in header): return header
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
    if not header.get('RESTFREQ') is None and not header.get('RESTFREQ') == 0:
        return header.get('RESTFREQ')
    elif not header.get('FREQ') is None and not header.get('FREQ') == 0:
        return header.get('FREQ')
    else:
        for i in xrange(5):
            type_s = header.get('CTYPE%i' % i)
            if type_s is not None and type_s[0:4] == 'FREQ':
                return header.get('CRVAL%i' % i)

    return None # no freq information found
 

class Image(object):

    def __init__(self, imagefile):

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

    def set_beam(self, beam):
        self.img_hdr['BMAJ'] = beam[0]
        self.img_hdr['BMIN'] = beam[1]
        self.img_hdr['BPA'] = beam[2]

    def get_beam(self):
        return [self.img_hdr['BMAJ'], self.img_hdr['BMIN'], self.img_hdr['BPA']]

    def get_wcs(self):
        return pywcs(self.img_hdr)

    def apply_region(self, regionfile, blankvalue=np.nan, invert=False):
        """
        Blank inside mask
        invert: blank outside mask
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

    def calc_noise(self, niter=50, eps=1e-5):
        """
        Return the rms of all the pixels in an image
        niter : robust rms estimation
        eps : convergency
        """
        data = self.img_data[ ~np.isnan(self.img_data) ] # remove nans
        oldrms = 1.
        for i in range(niter):
            rms = np.nanstd(data)
            if np.abs(oldrms-rms)/rms < eps:
                self.noise = rms
                logging.debug('%s: Noise: %.3f mJy/b' % (self.imagefile, self.noise*1e3))
                return

            data = data[np.abs(data)<5*rms]
            oldrms = rms
        raise Exception('Failed to converge')
