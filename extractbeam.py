#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013 - Francesco de Gasperin
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
#
#usage: extractbeam.py imagename

import pyrap.images as pim
from pyrap import quanta
import sys

for img in sys.argv[1:]:
    this_pim = pim.image(img)
    info_dict = this_pim.info()['imageinfo']['restoringbeam']
    # get beam info
    bpar_ma = quanta.quantity(info_dict['major']).get_value('arcsec')
    bpar_mi = quanta.quantity(info_dict['minor']).get_value('arcsec')
    bpar_pa = quanta.quantity(info_dict['positionangle']).get_value('deg')
    print('\n{0} - Beam: maj {1:0.3f} (arcsec), min {2:2.3f} (arcsec), pa {3:0.2f} (deg)'.format(img, bpar_ma, bpar_mi,bpar_pa))

