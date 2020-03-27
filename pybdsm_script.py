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

import bdsf
import sys

f = sys.argv[1]
try:
    freq = sys.argv[2]
except:
    freq = None

img = bdsf.process_image(f, advanced_opts=True, detection_image=f.replace('-pb',''), interactive=False, thresh_pix=5., thresh_isl=4., \
        adaptive_rms_box=True, rms_box_bright=(70,10), adaptive_thresh=6., frequency=freq)
img.write_catalog(format='fits',catalog_type='srl', clobber=True)
#img.write_catalog(format='ds9',catalog_type='srl', clobber=True)
#img.write_catalog(format='ds9',catalog_type='gaul', clobber=True)
#img.export_image(outfile=f+'_gaus_resid.fits', img_type='gaus_resid',clobber=True)
#img.export_image(outfile=f+'_gaus_model.fits', img_type='gaus_model',clobber=True)
img.export_image(outfile=f+'_rms.fits', img_type='rms',clobber=True)
#img.export_image(outfile=f+'_mean.fits', img_type='mean',clobber=True)
