#!/usr/bin/env python3
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

# reg2fits.py reg img.fits

from lib_fits import Image
import os, sys

im = Image(sys.argv[2])
regionfile = sys.argv[1]

im.apply_region(regionfile, blankvalue=1, invert=False)
im.apply_region(regionfile, blankvalue=0, invert=True)

im.write(sys.argv[1].replace('reg','fits'))
