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

class radioImage(Image):

    def __init__(self, imagefile):
        """
        imagefile: name of the fits file
        """
        Image.__init__(self,imagefile)
        self.region = None
        self.region_noise = None
        self.hdu = fits.PrimaryHDU(self.img_data, self.img_hdr)

        # calculate beam area in pixels
        cd1 = abs(self.img_hdr['CDELT1'])
        cd2 = abs(self.img_hdr['CDELT2'])
        bmaj = self.img_hdr['BMAJ']
        bmin = self.img_hdr['BMAJ']
        if ((cd1-cd2)/cd1)>1.0001 and ((bmaj-bmin)/bmin)>1.0001:
                raise RadioError('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

        gfactor = 2.0*np.sqrt(2.0*np.log(2.0))
        self.barea = 2.0*np.pi*(bmaj/cd1*bmin/cd2)/(gfactor*gfactor) # in pixels
        logging.info('Beam area is',self.barea,'pixels')

        self.masks = []


    def set_region(self, regionfile, individual=False):
        region = pyregion.open(regionfile).as_imagecoord(self.img_hdr)
        if individual:
            for region_split in region:
                self.masks.append( pyregion.ShapeList([region_split]).get_mask(hdu=self.hdu,shape=np.shape(self.img_data)) )
        else:
            self.masks.append( region.get_mask(hdu=self.hdu,shape=np.shape(self.img_data)) )


    def set_region_noise(self, regionfile):
        region_noise = pyregion.open(regionfile).as_imagecoord(self.img_hdr)
        self.mask_noise = region_noise.get_mask(hdu=self.hdu,shape=np.shape(self.img_data))

    def get_flux(self, nsigma = 0):
        """
        nsigma: use only pixels above this sigma
        """

        # set self.noise
        if self.region_noise is None:
            self.calc_noise()
        else:
            self.noise = np.nanstd( self.img_data[self.mask_noise] )


        fluxes = []
        errors = []
        for mask in self.masks:
            mask = np.logical_and(mask, ~np.isnan(self.img_data))
            flux = np.nansum( self.img_data[mask] ) / self.barea
            error = self.noise * np.sqrt( np.count_nonzero(mask) / self.barea )
            fluxes.append(flux)
            errors.append(error)

        return fluxes, errors
