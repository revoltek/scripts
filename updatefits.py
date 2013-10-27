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

# Usage: updatefits.py -setbeam max min pa -setkeyword name=value fitsfile

import pyfits,sys,optparse

opt = optparse.OptionParser(usage="%prog [-setbeam max,min,pa] [-setkeyword keyword=value] fitsfile", version="%prog 0.1")
opt.add_option('-b', '--setbeam', help='Set beam minaxis maxaxis and position angle to three comma-separated numbers (arcsec,arcsec,degree) [ex: 123,123,90]')
opt.add_option('-k', '--setkeyword', help='Set a keyword to a specific value')
opt.add_option('-s', '--show', action="store_true", help='Print headers', default=False)
(options, img) = opt.parse_args()
setbeam = options.setbeam
setkeyword = options.setkeyword
show = options.show
sys.stdout.flush()

try:
	hdulist = pyfits.open(img[0], mode='update')
except:
	print "ERROR: problems opening file "+img[0]
	sys.exit(1)

if ( show is True ):
	print hdulist[0].header

if ( not setkeyword is None ):
	try: keyword, value = setkeyword.split('=')
	except:
		print "ERROR: the format for \"--setkeyword\" is keyword=value"
		sys.exit(1)
	prihdr = hdulist[0].header
	print "Setting",keyword,"=",value
	prihdr.update(keyword, value)

if ( not setbeam is None ):
	try: bmaj,bmin,pa = setbeam.split(',')
	except:
		print "ERROR: the format for \"--setbeam\" is max,min,pa (arcsec,arcsec,deg)"
		sys.exit(1)
	prihdr = hdulist[0].header
	print "Setting beam to ",bmaj,bmin,pa
	bmaj = float(bmaj)/3600.
	bmin = float(bmin)/3600.
	pa = float(pa)
	prihdr.update('BMAJ', bmaj)
	prihdr.update('BMIN', bmin)
	prihdr.update('BPA', pa)

hdulist.flush()
print "Done!"
