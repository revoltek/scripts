#!/usr/bin/python
# casapy  --nogui --nologger -c ~/bin/scripts/image_armonizer.py image1 image2 image3 ...

import sys
import numpy as np
images = sys.argv[5:]

do_beam = True
do_regrid = True
do_cut = True
region_file = 'cut.rgn'


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
    print "# De Beam:"
    # armonize beams to the biggest
    for i, img in enumerate(images):
        print "Convolving (to", bmaxmaj, "arcsec):", img
        imsmooth(imagename=img, kernel='gauss', major=str(bmaxmaj)+'arcsec', minor=str(bmaxmaj)+'arcsec', pa='0deg', targetres=True, overwrite=True, outfile=img+'-conv'+str(bmaxmaj))
        images[i] = img+'-conv'+str(bmaxmaj)

if do_regrid:
    print "##########################################"
    print "# De Regrid:"
    # regrid to the first image size and pixel-size 1/5 of the beam
    newincr = qa.convert({'unit':'arcsec', 'value':bmaxmaj/5.},'rad')['value']
    print "Setting pixel to", newincr*180/np.pi/3600., "srcsec"
    for i, img in enumerate(images):
        print "Regridding:", img

        # pre-check: set ref pix to 0,0 and consequently the refval
        ia.open(img)
        cs = ia.coordsys()
        if i == 0:
            newrefv = cs.referencevalue()['numeric']
            newrefv[0] = ia.toworld([0,0])['numeric'][0]
            newrefv[1] = ia.toworld([0,0])['numeric'][1]
            cs.setreferencevalue(value=newrefv)
#            print "ref val", cs.referencevalue()

            newrefp = cs.referencepixel()['numeric']
            newrefp[0] = 0.
            newrefp[1] = 0.
            cs.setreferencepixel(value=newrefp)
#            print "ref pix", cs.referencepixel()

            cs.setincrement([-newincr,newincr])

            imrr = ia.regrid(outfile=img+'-regridded', csys=cs.torecord(), method='cubic', overwrite=True)

        else:
            refv = cs.referencevalue()['numeric']
            refv[0] = newrefv[0]
            refv[1] = newrefv[1]
            cs.setreferencevalue(value=refv)

            refp = cs.referencepixel()['numeric']
            refp[0] = newrefp[0]
            refp[1] = newrefp[1]
            cs.setreferencepixel(value=refp)

            cs.setincrement([-newincr,newincr])
            thisshape = ia.shape() #np.ones(len(cs.axiscoordinatetypes()))
            #thisshape[0] = shape[0]
            #thisshape[1] = shape[1]
            imrr = ia.regrid(outfile=img+'-regridded', csys=cs.torecord(), shape=thisshape, method='cubic', overwrite=True)

        images[i] = img+'-regridded'
        cs.done()
        ia.close()

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
