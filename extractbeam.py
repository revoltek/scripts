#!/usr/bin/python
#usage: extractbeam.py imagename

import pyrap.images as pim
from pyrap import quanta
import sys

this_pim = pim.image(sys.argv[1])
info_dict = this_pim.info()['imageinfo']['restoringbeam']
# get beam info
bpar_ma = quanta.quantity(info_dict['major']).get_value('arcsec')
bpar_mi = quanta.quantity(info_dict['minor']).get_value('arcsec')
bpar_pa = quanta.quantity(info_dict['positionangle']).get_value('deg')
print '\nmean Beam: {0:0.3f} maj (arcsec), {1:2.3f} min (arcsec), {2:0.2f} pa (deg)'.format(bpar_ma, bpar_mi,bpar_pa)

