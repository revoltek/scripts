#!/usr/bin/env python3
#i -*- coding: utf-8 -*-
#
# Copyright (C) 2022 - Francesco de Gasperin, Henrik Edler
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
import os, sys, logging, re, copy

from astropy.wcs import WCS as pywcs
from astropy.io import fits as pyfits
from astropy.cosmology import FlatLambdaCDM
from astropy.nddata import Cutout2D
from astropy.coordinates import match_coordinates_sky, SkyCoord, FK5
from astropy.convolution import Gaussian2DKernel
import pyregion
import astropy.units as u

def flatten(filename, channel=0, stokes=0):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    f = pyfits.open(filename)
    header_init = f[0].header

    naxis = header_init['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        pass
        #return f[0].header,f[0].data

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
    header["NAXIS1"]=f[0].header['NAXIS1']
    header["NAXIS2"]=f[0].header['NAXIS2']
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    dataslice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            dataslice.append(np.s_[:],)
        elif header_init['CTYPE%i'%i] == 'FREQ':
            #print("channel:", channel)
            dataslice.append(channel)
        elif header_init['CTYPE%i'%i] == 'STOKES':
            #print("stokes:", stokes)
            dataslice.append(stokes)
        else:
            dataslice.append(0)

    # add freq
    freq = find_freq(f[0].header)
    if freq is not None:
        header["FREQ"] = freq

    # add beam if present
    try:
        header["BMAJ"]=f[0].header['BMAJ']
        header["BMIN"]=f[0].header['BMIN']
        header["BPA"]=f[0].header['BPA']
    except:
        pass

    # slice=(0,)*(naxis-2)+(np.s_[:],)*2
    return header, f[0].data[tuple(dataslice)]


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
    elif not header.get('RFALPHA') is None and not header.get('RFALPHA') == 0:
        return header.get('RFALPHA')
    else:
        for i in range(5):
            type_s = header.get('CTYPE%i' % i)
            if type_s is not None and type_s[0:4] == 'FREQ':
                return header.get('CRVAL%i' % i)

    return None # no freq information found
 
class AllImages():

    def __init__(self, filenames, channel=0, stokes=0):
        if len(filenames) == 0:
            logging.error('Cannot find images!')
            raise ValueError()

        self.filenames = filenames
        self.images = []
        img_list, freqs = [], []
        for filename in filenames:
            img_list.append(Image(filename))
            freqs.append(img_list[-1].freq)
        self.images = [img_list[i] for i in np.argsort(freqs)]
        self.freqs = np.sort(freqs)


    def __len__(self):
        return len(self.images)

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        try:
            nextimage = self.images[self.index]
        except IndexError:
            raise StopIteration
        self.index += 1
        return nextimage

    def __getitem__(self, x):
        return self.images[x]

    def align_catalogue(self):
        [image.make_catalogue() for image in self]
        # ref cat - use lowest scaled noise as reference.
        noise = [image.calc_noise()*(image.get_freq()/(54.e6))**0.8 for image in self]
        ref_idx = np.argmin(noise)
        ref_cat = self[ref_idx].cat
        logging.info(f'Reference cat: {self[ref_idx].imagefile}')
        # keep only point sources
        target_beam = self.common_beam(circbeam=True)
        for i, image in enumerate(self):
            if i == ref_idx:
                # skip ref_cat
                continue
            # cross match
            idx_match, sep, _ = match_coordinates_sky(SkyCoord(ref_cat['RA'], ref_cat['DEC']), \
                                                      SkyCoord(image.cat['RA'], image.cat['DEC']))
            idx_matched_ref = np.arange(0, len(ref_cat))[sep < target_beam[0] * u.degree]
            idx_matched_img = idx_match[sep < target_beam[0] * u.degree]

            # find & apply shift
            if len(idx_match) < 5:
                logging.warning('%s: Not enough matches found, assume no shift.' % image.imagefile)
                continue

            dra = ref_cat['RA'][idx_matched_ref] - image.cat['RA'][idx_matched_img]
            dra[dra > 180] -= 360
            dra[dra < -180] += 360
            ddec = ref_cat['DEC'][idx_matched_ref] - image.cat['DEC'][idx_matched_img]
            flux = ref_cat['Peak_flux'][idx_matched_ref]
            image.apply_shift(np.average(dra, weights=flux), np.average(ddec, weights=flux))
            logging.info(f'{image.imagefile} Applying shift: {3600*np.average(dra, weights=flux):.5f}", {3600*np.average(ddec, weights=flux):.5f}".')
            # debug output
            matches_lsm = image.cat_lsm.copy()
            matches_lsm.select(idx_match, aggregate=True)
            matches_lsm.write(image.imagefile + '-matched.reg', format='ds9', clobber=True)
        # debug output
        ref_cat_lsm = self.images[ref_idx].cat_lsm
        ref_cat_lsm.write(image.imagefile + '-refcat.reg', format='ds9', clobber=True)

    def center_at(self, ra, dec):
        """
        Re-align all images to a common center
        Parameters
        ----------
        ra: float, Right ascension in deg
        dec: float, declination in deg
        """
        for image in self.images:
            image.apply_recenter_cutout(ra, dec)

    def common_beam(self, circbeam=True):
        """
        Return parameters of the smallest common beam
        Parameters
        ----------
        circbeam: bool, optional. Default True - force beam circular

        Returns
        -------
        bmaj, bmin, bpa: array-like. Common beam in deg

        """
        all_beams = [image.get_beam() for image in self.images]
        if all_beams.count(all_beams[0]) == len(all_beams):
            # all beams are already exactly the same
            return all_beams[0]

        if circbeam:
            maxmaj = np.max([image.get_beam()[0] for image in self.images])
            maxmin = np.max([image.get_beam()[1] for image in self.images])
            if maxmaj == maxmin:
                target_beam = [maxmaj, maxmaj, 0.]
            else:
                # add 1% to prevent crash in convolution when making only half of the beam larger
                target_beam = [maxmaj * 1.01, maxmaj * 1.01, 0.]
        else:
            from radio_beam import Beams
            my_beams = Beams([image.get_beam()[0] for image in self.images] * u.deg,
                             [image.get_beam()[1] for image in self.images] * u.deg,
                             [image.get_beam()[2] for image in self.images] * u.deg)
            common_beam = my_beams.common_beam()
            target_beam = [common_beam.major.to_value('deg'), common_beam.minor.to_value('deg'), common_beam.pa.to_value('deg')]
        return target_beam

    def convolve_to(self, beam=None, circbeam=False):
        """
        Convolve all images to a common beam. By default, convolve to smallest common beam.

        Parameters
        ----------
        beam: list, optional. Default = None
            Beam parameters [b_major, b_minor, b_pa] in [asec, asec, deg]. None: find smallest common beam
        circbeam: bool, optional. Default = False
            Force circular beam
        """

        if beam is None:
            target_beam = self.common_beam(circbeam=circbeam)
        else:
            target_beam = [beam[0] / 3600., beam[1] / 3600., beam[2]]
        logging.info('Final beam: %.1f" %.1f" (pa %.1f deg)' \
                     % (target_beam[0] * 3600., target_beam[1] * 3600., target_beam[2]))

        for image in self.images:
            image.convolve(target_beam)

    def regrid_common(self, size=None, region=None, pixscale=None, radec=None, square=False, action='regrid'):
        """
        Move all images to a common grid
        Parameters
        ----------
        size: float or array-like of size 2, optional. Default = None
            Size of the new grid in degree. If not a list of size two, is assumed to be square.
            If not provided, automatically determines the largest size that fits all images.
        region: ds9 region used to restrict the image to just cover it
        pixscale: float, optional. Default = derive from the beam of first image
            Size of a square pixel in arcseconds
        radec: RA [deg] and Dec [deg] where to chenter the final image, otherwise use first image
        square: bool, optional. Default = True
            If False, do not force square image.
        action: regrid, header, regrid_header
            The function can perform the regrid or just return the common header or both
        """
        
        rwcs = pywcs(naxis=2)
        rwcs.wcs.ctype = self.images[0].get_wcs().wcs.ctype
        if pixscale:
            cdelt = pixscale / 3600.
        else:
            cdelt = self.images[0].get_beam()[1] / 4.  # 1/4 of minor axes (deg)
        logging.info('Pixel scale: %f"' % (cdelt * 3600.))
        rwcs.wcs.cdelt = [-cdelt, cdelt]
        if radec:
            mra = radec[0]
            mdec = radec[1]
        else:
            midpix = np.array(self.images[0].img_data.shape)/2
            mra, mdec = self.images[0].get_wcs().all_pix2world(midpix[1], midpix[0], 0, ra_dec_order=True)
        rwcs.wcs.crval = [mra, mdec]

        if region:
            r = pyregion.open(region)
            mask = r.get_mask(header=self.images[0].img_hdr, shape=self.images[0].img_data.shape)
            intermediate = pyfits.PrimaryHDU(mask.astype(float), self.images[0].img_hdr)
            intermediate.writeto('__mask.fits', overwrite=True)
            w = self.images[0].get_wcs()
            y, x = mask.nonzero()
            ra_max, dec_max = w.all_pix2world(np.max(x), np.max(y), 0, ra_dec_order=True)
            ra_min, dec_min = w.all_pix2world(np.min(x), np.min(y), 0, ra_dec_order=True)
            # 2.2 to go a bit larger than the diameter
            size = [2.2*np.max( [ np.max([np.abs(dec_max-mdec),np.abs(dec_min-mdec)]), np.max([np.abs(ra_max-mra),np.abs(ra_min-mra)]) ] )]
            os.system('rm __mask.fits')
            #print(ra_min,ra_max,dec_min,dec_max,'size:',size)

        # Calculate sizes of all images to find smallest size that fits all images
        sizes = np.empty((len(self.images), 2))
        for i, image in enumerate(self.images):
            sizes[i] = np.array(image.img_data.shape) * image.get_degperpixel()
        if size:
            size = np.array(size)
            if np.any(np.min(sizes, axis=1) < np.min(size)):
                logging.warning(f'Requested size {size} is larger than smallest image size {np.min(sizes, axis=1)} in at least one dimension. This will result in NaN values in the regridded images.')
        else:
            size = np.min(sizes, axis=0)
            if square:
                size = np.min(size)
        try: # if we have square, make it 2 dim
            if len(size) == 1:
                size = np.repeat(size,2)
        except TypeError:
            size = np.repeat(size, 2)

        # fits file are ordered the opposite way: (dec,ra)
        ysize = int(np.rint(np.array([size[0]]) / cdelt))
        xsize = int(np.rint(np.array([size[-1]]) / cdelt))
        if xsize % 2 != 0: xsize += 1
        if ysize % 2 != 0: ysize += 1
        rwcs.wcs.crpix = [xsize / 2, ysize / 2]

        regrid_hdr = rwcs.to_header()
        regrid_hdr['NAXIS'] = 2
        regrid_hdr['NAXIS1'] = xsize
        regrid_hdr['NAXIS2'] = ysize
        regrid_hdr['EQUINOX'] = 2000.0
        regrid_hdr['RADESYSA'] = 'J2000'

        logging.info(f'Regridded image size: {size} deg ({ysize:.0f},{xsize:.0f} pixels))')
        if action == 'regrid' or action == 'regrid_header':
            for image in self.images:
                image.regrid(regrid_hdr)
        if action == 'header' or action == 'regrid_header':
            return regrid_hdr

    def suffix_exists(self, suffix):
        """ Check if suffix exists for all images"""
        return np.all([os.path.exists(name.replace('.fits', f'-{suffix}.fits')) for name in self.filenames])

    def write(self, suffix, inflate=False):
        """ Write all (changed) images to imagename-suffix.fits"""
        for image in self.images:
            image.write(image.imagefile.replace('.fits', f'-{suffix}.fits'), inflate=inflate)


class Image(object):

    def __init__(self, imagefile, channel=0, stokes=0):
        """
        imagefile: name of the fits file
        """

        logging.info(f"Open {imagefile}")
        self.imagefile = imagefile
        header = pyfits.open(imagefile)[0].header
        header = correct_beam_header(header)
        self.img_hdr_orig = header
        
        try:
            self.img_hdr_orig["EQUINOX"]
        except KeyError:
            try:
                self.img_hdr_orig["EQUINOX"] = self.img_hdr_orig["EPOCH"]
            except KeyError:
                self.img_hdr_orig["EQUINOX"] = 2000.0 

        if self.img_hdr_orig["EQUINOX"] != 2000:
            logging.warning(f'Equinox is not 2000, but {self.img_hdr_orig["EQUINOX"]}. transforming to J2000')
            ra_b1950 = self.img_hdr_orig["CRVAL1"] # RA in degrees
            dec_b1950 = self.img_hdr_orig["CRVAL2"] # Dec in degrees
            
            coord_b1950 = SkyCoord(ra=ra_b1950, dec=dec_b1950, unit='deg', frame=FK5, equinox="B1950")
            coord_j2000 = coord_b1950.transform_to(FK5(equinox="J2000"))
            
            self.img_hdr_orig["EQUINOX"] = 2000.0
            self.img_hdr_orig["CRVAL1"] = coord_j2000.ra.deg
            self.img_hdr_orig["CRVAL2"] = coord_j2000.dec.deg
            
        try:
            beam = [header['BMAJ'], header['BMIN'], header['BPA']]
        except:
            logging.warning('%s: No beam information found.' % self.imagefile)
            beam = [0,0,0]

        logging.debug('%s: Beam: %.1f" %.1f" (pa %.1f deg)' % \
                (self.imagefile, beam[0]*3600., beam[1]*3600., beam[2]))

        freq = find_freq(header)
        if freq is None:
            logging.warning('%s: No frequency information found.' % self.imagefile)
            # sys.exit(1)
        else:
            logging.debug('%s: Frequency: %.0f MHz' % (self.imagefile, freq/1e6))
        try:
            self.mhz = np.round(freq / 1e6).astype(int) # handy for suffixes etc.
        except TypeError:
            self.mhz = None

        self.noise = None
        self.img_hdr, self.img_data = flatten(self.imagefile, channel=channel, stokes=stokes)
        self.img_data_orig = self.img_data.copy()
        
        try:
            self.img_hdr["EQUINOX"]
        except KeyError:
            try:
                self.img_hdr["EQUINOX"] = self.img_hdr["EPOCH"]
            except KeyError:
                self.img_hdr["EQUINOX"] = 2000.0 

        if self.img_hdr["EQUINOX"] != 2000:
            ra_b1950 = self.img_hdr["CRVAL1"] # RA in degrees
            dec_b1950 = self.img_hdr["CRVAL2"] # Dec in degrees
            
            coord_b1950 = SkyCoord(ra=ra_b1950, dec=dec_b1950, unit='deg', frame=FK5, equinox="B1950")
            coord_j2000 = coord_b1950.transform_to(FK5(equinox="J2000"))
            
            self.img_hdr["EQUINOX"] = 2000.0
            self.img_hdr["CRVAL1"] = coord_j2000.ra.deg
            self.img_hdr["CRVAL2"] = coord_j2000.dec.deg
        
        self.img_hdu = pyfits.ImageHDU(data=self.img_data, header=self.img_hdr)
        self.set_beam(beam)
        self.set_freq(freq)
        self.ra = self.img_hdr['CRVAL1']
        self.dec = self.img_hdr['CRVAL2']
        self.get_degperpixel() # sets self.degperpixel (faster to call, no WCS call)



    def write(self, filename=None, inflate=False):
        """
        Write to fits-file
        Parameters
        ----------
        filename: str, filename
        inflate: bool, optional. Default=False
                If False, write as flat 2D-fits file. If true, inflate to 4D.
        """
        if filename is None:
            filename = self.imagefile
        if inflate:
            # Inflate a fits file so that it becomes a 4D image. Return new header and data
            hdr_inf = pyfits.Header()
            hdr_inf['SIMPLE'  ] = self.img_hdr_orig['SIMPLE']
            hdr_inf['BITPIX'  ] = self.img_hdr_orig['BITPIX']
            hdr_inf['NAXIS'   ] = 4
            hdr_inf['NAXIS1'  ] = self.img_hdr['NAXIS1']
            hdr_inf['NAXIS2'  ] = self.img_hdr['NAXIS2']
            hdr_inf['NAXIS3'  ] = 1
            hdr_inf['NAXIS4'  ] = 1
            hdr_inf['EXTEND'  ] = 'T'
            hdr_inf['BUNIT'   ] = 'JY/BEAM'
            hdr_inf['RADESYS' ] = 'FK5'
            hdr_inf['EQUINOX' ] = 2000.
            hdr_inf['BMAJ'    ] = self.img_hdr['BMAJ']
            hdr_inf['BMIN'    ] = self.img_hdr['BMIN']
            hdr_inf['BPA'     ] = self.img_hdr['BPA']
            hdr_inf['EQUINOX' ] = self.img_hdr_orig['EQUINOX']
            hdr_inf['BTYPE'   ] = 'INTENSITY'
            hdr_inf['TELESCOP'] = self.img_hdr_orig['TELESCOP']
            hdr_inf['OBJECT'  ] = self.img_hdr_orig['OBJECT']
            hdr_inf['CTYPE1'  ] = self.img_hdr['CTYPE1']
            hdr_inf['CRPIX1'  ] = self.img_hdr['CRPIX1']
            hdr_inf['CRVAL1'  ] = self.img_hdr['CRVAL1']
            hdr_inf['CDELT1'  ] = self.img_hdr['CDELT1']
            hdr_inf['CUNIT1'  ] = self.img_hdr['CUNIT1']
            hdr_inf['CTYPE2'  ] = self.img_hdr['CTYPE2']
            hdr_inf['CRPIX2'  ] = self.img_hdr['CRPIX2']
            hdr_inf['CRVAL2'  ] = self.img_hdr['CRVAL2']
            hdr_inf['CDELT2'  ] = self.img_hdr['CDELT2']
            hdr_inf['CUNIT2'  ] = self.img_hdr['CUNIT2']
            hdr_inf['CTYPE3'  ] = 'FREQ'
            hdr_inf['CRPIX3'  ] = 1.
            hdr_inf['CRVAL3'  ] = self.get_freq()
            hdr_inf['CDELT3'  ] = 10000000.
            hdr_inf['CUNIT3'  ] = 'Hz'
            hdr_inf['CTYPE4'  ] = 'STOKES'
            hdr_inf['CRPIX4'  ] = 1.
            hdr_inf['CRVAL4'  ] = 1.
            hdr_inf['CDELT4'  ] = 1.
            hdr_inf['CUNIT4'  ] = ' '
            pyfits.writeto(filename, self.img_data[np.newaxis,np.newaxis], hdr_inf, overwrite=True, output_verify='fix')
        else:
            pyfits.writeto(filename, self.img_data, self.img_hdr, overwrite=True, output_verify='fix')

    def set_beam(self, beam):
        self.img_hdr['BMAJ'] = beam[0]
        self.img_hdr['BMIN'] = beam[1]
        self.img_hdr['BPA'] = beam[2]

    def get_beam(self):
        return [self.img_hdr['BMAJ'], self.img_hdr['BMIN'], self.img_hdr['BPA']]

    def set_freq(self, freq):
        if freq:
            self.img_hdr['RESTFREQ'] = freq
            self.img_hdr['FREQ'] = freq
            self.freq = freq

    def get_freq(self):
        try:
            return self.img_hdr['FREQ']
        except:
            return None

    def get_beam_area(self, unit='arcsec'):
        """
        Return area of psf.
        Parameters
        ----------
        unit: string, optional. Default = arcsec.
            Units in which to return the area. Either arcsec or pixel
        Returns
        -------
        beam area: float
        """
        b = self.get_beam()
        beam_area_squaredeg = (2*np.pi * b[0] * b[1]) / (8. * np.log(2.))
        if unit in ['arcsec', 'asec']:
            return beam_area_squaredeg * 3600**2
        elif unit in ['pix','pixel']:
            return beam_area_squaredeg / self.degperpixel**2


    def get_wcs(self):
        return pywcs(self.img_hdr)

    def apply_region(self, regionfile, blankvalue=np.nan, invert=False):
        """
        Blank inside mask
        invert: blank outside region
        """
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

    def blank_noisy(self, nsigma):
        """
        Set to nan pixels below nsigma*noise
        """
        nans_before = np.sum(np.isnan(self.img_data))
        self.img_data[np.isnan(self.img_data)] = 0  # temporary set nans to 0 to prevent error in "<"
        self.img_data[np.where(self.img_data <= nsigma * self.noise)] = np.nan
        nans_after = np.sum(np.isnan(self.img_data))
        logging.debug('%s: Blanked pixels %i -> %i' % (self.imagefile, nans_before, nans_after))

    def calc_noise(self, niter=1000, eps=None, sigma=5, bg_reg=None, force_recalc=False):
        """
        Return the rms of all the pixels in an image
        niter : robust rms estimation
        eps : convergency criterion, if None is 1% of initial rms
        bg_reg : If ds9 region file provided, use this as background region
        force_recalc : recalculate noise even if already set
        """
        if not self.noise is None and force_recalc == False:
            print('WARNING: Noise already set, and force_recalc=False.')
            return self.noise

        if bg_reg is not None:
            if not os.path.exists(bg_reg):
                logging.error('%s: Region file not found.' % bg_reg)
                sys.exit(1)

            logging.debug('%s: Apply background region %s' % (self.imagefile, bg_reg))
            r = pyregion.open(bg_reg)
            mask = r.get_mask(header=self.img_hdr, shape=self.img_data.shape)
            self.noise = np.nanstd(self.img_data[mask])
            logging.info('%s: Noise: %.3f mJy/b' % (self.imagefile, self.noise * 1e3))
        else:
            from astropy.stats import median_absolute_deviation
            if eps == None: eps = 1e-3
            data = self.img_data[ ~np.isnan(self.img_data) & (self.img_data != 0) ] # remove nans and 0s
            initial_len = len(data)
            if initial_len == 0: return 0
            mad_old = 0.
            for i in range(niter):
                 mad = median_absolute_deviation(data)
                 logging.debug('%s: MAD noise: %f uJy on %f%% data' % (self.imagefile, mad*1e6, 100*len(data)/initial_len))
                 if np.isnan(mad): break
                 if np.abs(mad_old-mad)/mad < eps:
                     rms = np.nanstd( data )
                     self.noise = rms
                     #print('%s: Noise: %.3f mJy/b' % (self.imagefile, self.noise*1e3))
                     logging.debug('%s: Noise: %.3f mJy/b (data len: %i -> %i - %.2f%%)' % (self.imagefile, self.noise*1e3, initial_len, len(data), 100*len(data)/initial_len))
                     return rms
    
                 data = data[np.abs(data) < (sigma*mad)]
                 mad_old = mad
            raise Exception('Noise estimation failed to converge.')


    def convolve(self, target_beam, stokes=True):
        """
        Convolve *to* this rsolution
        beam = [bmaj, bmin, bpa]
        """
        from lib_beamdeconv import deconvolve_ell, EllipticalGaussian2DKernel
        from astropy import convolution
        
        # if difference between beam is negligible <1%, skip - it mostly happens when beams are exactly the same
        beam = self.get_beam()
        try:
            if (np.abs((target_beam[0]/beam[0])-1) < 1e-2) and (np.abs((target_beam[1]/beam[1])-1) < 1e-2) and (np.abs(target_beam[2] - beam[2]) < 1):
                logging.debug('%s: do not convolve. Same beam.' % self.imagefile)
                return
            elif (target_beam[0] < beam[0]) or target_beam[1] < beam[1]:
                raise ValueError('%s: target beam is smaller than current beam. Cannot convolve!.' % self.imagefile)
        except ZeroDivisionError: pass  # catch case where we have delta scale beam (model image)
        # first find beam to convolve with
        convolve_beam = deconvolve_ell(target_beam[0], target_beam[1], target_beam[2], beam[0], beam[1], beam[2])
        if convolve_beam[0] is None:
            logging.error('Cannot deconvolve this beam.')
            sys.exit(1)
        logging.debug('%s: Convolve beam: %.3f" %.3f" (pa %.1f deg)' \
                % (self.imagefile, convolve_beam[0]*3600, convolve_beam[1]*3600, convolve_beam[2]))
        # do convolution on data
        bmaj, bmin, bpa = convolve_beam
        #print(self.imagefile,self.img_hdr['CDELT1'], self.img_hdr['CDELT2'])
        
        assert abs(self.img_hdr['CDELT1']) - abs(self.img_hdr['CDELT2']) < 1e-6
        pixsize = abs(self.img_hdr['CDELT1'])
        fwhm2sigma = 1./np.sqrt(8.*np.log(2.))
        gauss_kern = EllipticalGaussian2DKernel((bmaj*fwhm2sigma)/pixsize, (bmin*fwhm2sigma)/pixsize, (90+bpa)*np.pi/180.) # bmaj and bmin are in pixels
        self.img_data = convolution.convolve(self.img_data, gauss_kern, boundary=None, preserve_nan=True)
        if stokes: # if not stokes image (e.g. spectral index, do not renormalize)
            self.img_data *= (target_beam[0]*target_beam[1])/(beam[0]*beam[1]) # since we are in Jy/b we need to renormalise

        self.set_beam(target_beam) # update beam

    def regrid(self, regrid_hdr):
        """ Regrid image to new header """
        from reproject import reproject_interp, reproject_exact
        reproj = reproject_exact
        # store some info so to reconstruct headers
        beam = self.get_beam()
        freq = self.get_freq()
        logging.debug('%s: regridding' % (self.imagefile))
        self.img_data, __footprint = reproj((self.img_data, self.img_hdr), regrid_hdr, parallel=True)
        # update headers
        self.img_hdr = copy.copy(regrid_hdr)
        self.set_freq(freq)
        self.set_beam(beam)
        self.get_degperpixel() # update

    def apply_shift(self, dra, ddec):
        """
        Shift header by dra/ddec
        dra, ddec in degree
        """
        # correct the dra shift for np.cos(DEC*np.pi/180.) -- only in the log as the reference val doesn't need it!
        logging.info('Shift %.2f %.2f arcsec (%s)' % (dra*3600*np.cos(self.dec*np.pi/180.), ddec*3600, self.imagefile))
        dec = self.img_hdr['CRVAL2']
        self.img_hdr['CRVAL1'] += dra
        self.img_hdr['CRVAL2'] += ddec

    def pixel_covariance(self, pix1, pix2):
        """
        Get the covariance matrix
        Source Finding in the Era of the SKA (Precursors): AEGEAN2.0 -- Hancock, Trott, Hurley-Walker
        Parameters
        ----------
        pix1
        pix2

        Returns
        -------

        """
        b = self.get_beam()
        fwhm2sigma = 1./np.sqrt(8.*np.log(2.))
        b_sig_pix_ma = b[0] * fwhm2sigma / self.degperpixel
        b_sig_pix_min = b[1] * fwhm2sigma / self.degperpixel
        theta = np.deg2rad(b[2])
        dx, dy = np.array(pix1 - pix2)

        uncorrelated_variance = self.noise**2 # TODO which is the best mode to find noise?

        # NOTE: It might very well be that there is an error in the definition of the angle here! So maybe +pi/2 or -theta would be correct.
        if np.linalg.norm([dx, dy]) > 3*b_sig_pix_ma:
            Cij = 0.
        else:
            # According to AEGEAN2.0 there is a factor of 1/4 in the exponent. But that is actually the wrong way around, since we need the square of the 2DGaussian and not the squareroot?
            Cij = np.exp(-((dx*np.sin(theta)+dy*np.cos(theta))/(b_sig_pix_ma))**2
                         -((dx*np.cos(theta)-dy*np.sin(theta))/(b_sig_pix_min))**2)
        Cij *= uncorrelated_variance
        #print(Cij)
        return Cij

    def make_catalogue(self):
        """
        Create catalogue for image alignmnt
        """
        import bdsf
        import lsmtool as lsm
        from astropy.table import Table

        img_cat = self.imagefile+'.cat'
        if not os.path.exists(img_cat):
            bdsf_img = bdsf.process_image(self.imagefile, rms_box=(100,30), \
                                          thresh_pix=5, thresh_isl=3, atrous_do=False, \
                                          adaptive_rms_box=True, adaptive_thresh=100, rms_box_bright=(30,10), quiet=True)
            bdsf_img.write_catalog(outfile=img_cat, catalog_type='srl', format='fits', clobber=True)
            bdsf_img.write_catalog(outfile=img_cat.replace('.cat', '.skymodel'), catalog_type='gaul', format='bbs', bbs_patches='source', clobber=True, srcroot='src')
        else:
            logging.warning('%s already exists, using it.' % img_cat)

        cat = Table.read(img_cat)
        # remove extended sources
        extended_src = (cat['Peak_flux'] / cat['Total_flux']) < 0.1 # ~extended source
        extended_src[cat['S_Code'] == 'M'] = True # multiple-gaussian source
        extended_src[cat['S_Code'] == 'C'] = True # one gaussian + other sources island
        # remove same sources from skymodel
        cat_lsm = lsm.load(img_cat.replace('.cat', '.skymodel'))
        # TODO this spams the logging
        for srcid in cat[extended_src]['Source_id']:
            cat_lsm.remove(f'Patch == src_patch_s{srcid}')
        cat.remove_rows(np.argwhere(extended_src))
        self.cat = cat
        self.cat_lsm = cat_lsm
        logging.debug('%s: Number of sources detected: %i; removed %i extended sources.' % (self.imagefile, len(self.cat), sum(extended_src)) )

    def get_degperpixel(self):
        """
        Return the number of degrees per image pixel. This assumes SQUARE pixels!
        Returns
        -------
        degperpixel: float
        """
        wcs = self.get_wcs()
        self.degperpixel = np.abs(wcs.all_pix2world(0, 0, 0)[1] - wcs.all_pix2world(0, 1, 0)[1])
        return self.degperpixel

    def get_degperkpc(self, z):
        """
        How many degrees are there per kpc? Assume H0=70km/S/Mpcm O_m = 0.3

        Parameters
        ----------
        z: Source redshift

        Returns
        -------
        degperkpc: float
        """
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        return cosmo.arcsec_per_kpc_proper(z).value / 3600.


    def get_pixelperkpc(self, z):
        """
        Return the number of pixel per kpc. This assumes SQUARE pixels!
        Returns
        -------
        pixelperkpc: float
        """
        return self.get_degperkpc(z) / self.get_degperpixel()

