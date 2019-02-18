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
# Usage: casapy  --nogui --nologger -c ~/bin/scripts/casa_reg2mask.py image region
# output is region.mask

import sys
import numpy as np
img = sys.argv[5]
region = sys.argv[6]

if os.path.exists(region+'.mask'): os.system('rm -r '+region+'.mask')
os.system('cp -r '+img+' '+region+'.mask')
ia.open(region+'.mask')
reg = rg.fromtextfile(filename=region, shape=ia.shape(), csys=ia.coordsys().torecord())
ia.close()
ia.open(region+'.mask')
ia.set(pixels='0')
ia.set(pixels='1', region=reg)
print("Output: "+region+".mask")
ia.close()

