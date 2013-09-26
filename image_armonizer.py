#!/usr/bin/python
# casapy  --nogui --nologger -c ~/bin/image_armonizer.py image1 image2 image3 ...

import sys
import numpy as np
images = sys.argv[5:]

do_beam = True
do_regrid = False

# load images
bmaxmaj = 0
for img in images:
	bmaj = qa.convert(imhead(imagename=img,mode='get',hdkey='beammajor'),'arcsec')
	bmin = qa.convert(imhead(imagename=img,mode='get',hdkey='beamminor'),'arcsec')
	bpa = imhead(imagename=img,mode='get',hdkey='beampa')
	if bmaxmaj < bmaj['value']:
		print "Setting bmaj to ", bmaj
		bmaxmaj = bmaj['value']
		bmaxmaj_img = img
	print "Load image:", img, "(Major:", bmaj['value'], bmaj['unit'], "- Minor:", bmin['value'], bmin['unit'], ")"

if do_beam:
    # armonize beams to the biggest
    for i, img in enumerate(images):
	    print "Convolving", img, "(to", bmaxmaj, "arcsec)"
	    imsmooth(imagename=img, kernel='gauss', major=str(bmaxmaj)+'arcsec', minor=str(bmaxmaj)+'arcsec', pa='0deg', targetres=True, overwrite=True, outfile=img+'-convolved'+str(bmaxmaj))
	    images[i] = img+'-convolved'+str(bmaxmaj)

# find smallest image
#for img in images:
#	print "Scanning", img
#	ia.open(img)
#	cs = ia.coordsys()
#	axisra = cs.findcoordinate(type='direction')[1][0]
#	axisdec = cs.findcoordinate(type='direction')[1][1]
#	ra = ia.toworld([0,0])[axisra]
#	dec = ia.toworld([0,0])[axisdec]
#	ra = ia.toworld([ia.shape()[axisra],ia.shape()[axisdec]])[axisra]
#	dec = ia.toworld([ia.shape()[axisra],ia.shape()[axisdec]])[axisdec]
#	if ra < minra

if do_regrid:
    # regrid to the first image coord sys
    newincr = qa.convert({'unit':'arcsec', 'value':bmaxmaj/3.},'rad')['value']
    for i, img in enumerate(images):
    	print "Working on", img
    	print "Setting pixel to", newincr, "rad"
    	# pre-check: set ref pix to 0,0 and consequently the refval
    	ia.open(img)
    	cs = ia.coordsys()
    	refv = cs.referencevalue()
    	refv['numeric'][0] = ia.toworld([0,0])[0]
    	refv['numeric'][1] = ia.toworld([0,0])[1]
    	refp = cs.referencepixel()
    	refp['numeric'][0] = 0
    	refp['numeric'][1] = 0
    	imrr = ia.regrid(outfile=img+'-preregridded', csys=cs.torecord(), method='cubic')
    	# regridding
    	ia.open(img+'-preregridded')
    	cs = ia.coordsys()
    	cs.setincrement([-newincr,newincr])
    	#print "ref val", cs.referencevalue()
    	#print "ref pix", cs.referencepixel()
    	if i == 0:
    		imrr = ia.regrid(outfile=img+'-regridded', csys=cs.torecord(), method='cubic')
    		cs.done()
    		ia.close()
    		ia.open(img+'-regridded')
    		cs = ia.coordsys()
    		shape = ia.shape()
    		newrefv = cs.referencevalue()
    		newrefp = cs.referencepixel()
    	else:
    		refv = cs.referencevalue()
    		refv['numeric'][0] = newrefv['numeric'][0]
    		refv['numeric'][1] = newrefv['numeric'][1]
    		cs.setreferencevalue(refv)
    		refp = cs.referencepixel()
    		refp['numeric'][0] = newrefp['numeric'][0]
    		refp['numeric'][1] = newrefp['numeric'][1]
    		cs.setreferencepixel(refp)
    		thisshape = np.ones(len(cs.axiscoordinatetypes()))
    		thisshape[0] = shape[0]
    		thisshape[1] = shape[1]
    		imrr = ia.regrid(outfile=img+'-regridded', csys=cs.torecord(), shape=thisshape, method='cubic')
    	cs.done()
    	ia.close()
