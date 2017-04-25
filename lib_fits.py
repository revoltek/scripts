from astropy.wcs import WCS as pywcs
from astropy.io import fits as pyfits
import numpy as np

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

class Image(object):

    def __init__(self, imagefile):

        self.fitsfile = imagefile
        self.freq = 
        self.img_hdr, self.img_data = flatten(self.imagefile)
        self.beam = [img_hdr['BMAJ'], img_hdr['BMIN'], img_hdr['BPA']]
	self.noise = None
        
    def get_wcs(self):
        return pywcs(self.img_hdr)

    def applyRegion(self, regionfile, blankvalue=np.nan, invert=False)
        """
        Blank inside mask
        invert: blank outside mask
        """
        logging.debug('%s - apply region' % (self.fistfile))
        r = pyregion.open(regionfile)
        mask = r.get_mask(header=self.img_hdr, shape=self.img_data.shape)
        if invert: self.img_data[~mask] = blankvalue
        else: self.img_data[mask] = blankvalue

    def calc_noise(self, niter=20, eps=1e-5):
        """
        Return the rms of all the pixels in an image
        niter : robust rms estimation
        eps : convergency
        """
        with pyfits.open(self.imagefile) as fits:
            data = fits[0].data
            oldrms = 1.
            for i in range(niter):
                rms = np.nanstd(data)
                if np.abs(oldrms-rms)/rms < eps:
                    self.noise = rms
                    logging.debug('Noise for %s: %f' % (self.imagefile, self.noise))
                    return
                data = data[np.abs(data)<5*rms]
                oldrms = rms
            raise Exception('Failed to converge')
