#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2020 - Francesco de Gasperin
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


import os, sys, logging, re
import numpy as np
from astropy.io import fits
import pyregion
from lib_fits import Image

class RadioImage(Image):

    def __init__(self, imagefile):
        """
        imagefile: name of the fits file
        """
        Image.__init__(self,imagefile)
        self.hdu = fits.PrimaryHDU(self.img_data, self.img_hdr)
        self.region_noise = None

        # calculate beam area in pixels
        cd1 = abs(self.img_hdr['CDELT1'])
        cd2 = abs(self.img_hdr['CDELT2'])
        bmaj = self.img_hdr['BMAJ']
        bmin = self.img_hdr['BMIN']
        if ((cd1-cd2)/cd1)>1.0001 and ((bmaj-bmin)/bmin)>1.0001:
                raise RadioError('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

        gfactor = 2.0*np.sqrt(2.0*np.log(2.0))
        self.barea = 2.0*np.pi*(bmaj/cd1*bmin/cd2)/(gfactor*gfactor) # in pixels
        logging.info('Beam area is',self.barea,'pixels')

        self.masks = []


    def set_region(self, regionfile, individual=False):
        self.masks = []
        region = pyregion.open(regionfile).as_imagecoord(self.img_hdr)
        if individual:
            for region_split in region:
                self.masks.append( pyregion.ShapeList([region_split]).get_mask(hdu=self.hdu,shape=np.shape(self.img_data)) )
        else:
            self.masks.append( region.get_mask(hdu=self.hdu,shape=np.shape(self.img_data)) )


    def set_region_noise(self, regionfile):
        self.region_noise = regionfile

    def get_flux(self, flux_scale_err = 0, nsigma = 0, with_upper_limits = False, upper_limit_sigma = 3):
        """
        flux_scale_err = percentage of error to add to the flux scale (e.g. 0.05 to have a 5% error)
        nsigma: use only pixels above this sigma
        with_upper_limits: if no detection, set the value at upper_limit_sigma sigma. It also returns a bool array with True for limits
        upper_limit_sigma: numer of sigmas to consider a flux a limit (default: 3)
        """

        # set self.noise
        if self.region_noise is None:
            self.calc_noise()
        else:
            self.calc_noise(bg_reg = self.region_noise)

        fluxes = []
        errors = []
        for mask in self.masks:
            mask = np.logical_and(mask, ~np.isnan(self.img_data))
            flux = np.nansum( self.img_data[mask] ) / self.barea
            error_rms = self.noise * np.sqrt( np.count_nonzero(mask) / self.barea )
            error_flux = flux_scale_err*flux
            error = np.sqrt(error_rms**2+error_flux**2)
            fluxes.append(flux)
            errors.append(error)

        fluxes = np.array(fluxes)
        errors = np.array(errors)

        if with_upper_limits:
            upper_limits = np.zeros_like(fluxes, dtype=bool)
            is_limit = np.where(fluxes < (upper_limit_sigma * errors))
            upper_limits[ is_limit ] = True
            fluxes[ is_limit ] = upper_limit_sigma*errors[ is_limit ]
            return fluxes, errors, upper_limits
        else:
            return fluxes, errors
