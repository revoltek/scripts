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

# Usage: casapy  --nogui --nologger -c ~/bin/scripts/image_armonizer.py image1 image2 image3 ...

import sys
import numpy as np
images = sys.argv[5:]

do_beam = True
do_regrid = True
newincr = 5 # arcsec
do_cut = True
region_file = 'cut.crtf'
to_fits = False
clean = True

todelete = []

print "##########################################"
print "# Converting:"
for i, img in enumerate(images):
    # if it is a fits file transform it into an MS image
    if img[-5:] == '.fits':
        print "Convert to MS: "+img
        importfits(fitsimage=img, overwrite=True, imagename=img.replace('.fits','.image'))
        images[i] = img.replace('.fits','.image')

print "##########################################"
print "# Compute Beam:"
bmaxmaj = 0
for i, img in enumerate(images):
    bmaj = qa.convert(imhead(imagename=img,mode='get',hdkey='beammajor'),'arcsec')
    bmin = qa.convert(imhead(imagename=img,mode='get',hdkey='beamminor'),'arcsec')
    bpa = imhead(imagename=img,mode='get',hdkey='beampa')
    print img, ": Major:", bmaj['value'], bmaj['unit'], "- Minor:", bmin['value'], bmin['unit']
    if bmaxmaj < bmaj['value']:
        bmaxmaj = bmaj['value']+0.01
        bmaxmaj_img = img
print "Max beam is: ", bmaxmaj

if do_beam:
    print "##########################################"
    print "# Do Beam:"
    # armonize beams to the biggest
    for i, img in enumerate(images):
        print "Convolving (to", bmaxmaj, "arcsec):", img
        imsmooth(imagename=img, kernel='gauss', major=str(bmaxmaj)+'arcsec', minor=str(bmaxmaj)+'arcsec', pa='0deg', targetres=True, overwrite=True, outfile=img+'-conv'+str(bmaxmaj))
        images[i] = img+'-conv'+str(bmaxmaj)
        todelete.append(img)


if do_cut:
    print "##########################################"
    print "# Do Cut:"
    print "Region:", region_file
    for i, img in enumerate(images):
        print "Cutting:", img
        ia.open(img)   # casa image class
        #ia.summary()   # will print out some info
        #box = rg.box([10,10], [50,50])   # casa regionmanager class
        im2 = ia.subimage(outfile=img+'-cut', region=region_file, overwrite=True)
        images[i] = img+'-cut'
        todelete.append(img)
        ia.close()

if do_regrid:
    print "##########################################"
    print "# Do Regrid:"
    # regrid to the first image size and pixel-size 1/5 of the beam if not set
    print "Setting pixel to", newincr, "arcsec"
    if not type(newincr) is int:
        newincr = bmaxmaj/5.
    newincr = qa.convert({'unit':'arcsec', 'value':newincr},'rad')['value']

    for i, img in enumerate(images):
        print "Regridding:", img

        # pre-check: set ref pix to 0,0 and consequently the refval
        ia.open(img)
        cs = ia.coordsys()
        if i == 0:
            newrefv = cs.referencevalue()['numeric']
            newrefv[0] = ia.toworld([0,0])['numeric'][0]
            newrefv[1] = ia.toworld([0,0])['numeric'][1]
#            print "ref val", cs.referencevalue()

            newrefp = cs.referencepixel()['numeric']
            newrefp[0] = 0.
            newrefp[1] = 0.
#            print "ref pix", cs.referencepixel()

            # scale the shape according to the new increment per pixel
            newshape = ia.shape()
            rescaling = cs.increment()['numeric'][1]/newincr
            newshape[0] = ia.shape()[0]*rescaling
            newshape[1] = ia.shape()[1]*rescaling

        refv = cs.referencevalue()['numeric']
        refv[0] = newrefv[0]
        refv[1] = newrefv[1]
        cs.setreferencevalue(value=refv)

        refp = cs.referencepixel()['numeric']
        refp[0] = newrefp[0]
        refp[1] = newrefp[1]
        cs.setreferencepixel(value=refp)

        cs.setincrement([-newincr,newincr])
        cs.setprojection(type='SIN', parameters=[0,0])
        shape = ia.shape()
        shape[0] = newshape[0]
        shape[1] = newshape[1]
        imrr = ia.regrid(outfile=img+'-regridded', csys=cs.torecord(), shape=shape, method='cubic', overwrite=True)

        images[i] = img+'-regridded'
        todelete.append(img)
        cs.done()
        ia.close()
        # if missing beam: use ia.setrestoringbeam()

if do_cut:
    print "##########################################"
    print "# Do Cut:"
    print "Region:", region_file
    for i, img in enumerate(images):
        print "Cutting:", img
        ia.open(img)   # casa image class
        #ia.summary()   # will print out some info
        #box = rg.box([10,10], [50,50])   # casa regionmanager class
        im2 = ia.subimage(outfile=img+'-cut', region=region_file, overwrite=True)
        images[i] = img+'-cut'
        todelete.append(img)
        ia.close()

if to_fits:
    print "##########################################"
    print "# Do To_FITS:"
    for i, img in enumerate(images):
        os.system('image2fits in='+img+' out='+img+'.fits')
        todelete.append(img)

if clean:
    print "##########################################"
    print "# Cleaning up..."
    os.system('rm *log')
    for img in todelete:
        os.system('rm -r '+img)

print "Done."