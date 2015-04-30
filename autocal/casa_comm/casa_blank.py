#!/usr/bin/python
# casapy --nogui --log2term --nologger -c  casa_blank.py picklefile
# blank all the pixel outside this region
# NOTE: overwrite the image! And the blanked image is the smallest around the given region.

import os, sys, pickle

def casa_blank(imgs = [], region = '', inverse=False, setTo = 0):
    """
    img = image to blank
    region = blank all pixels inside this region to the value "setTo"
    setTo = value to set pixels
    inverse = instead of blanking pixel return an image that covers only the "region" region
    out: the new, blanked image will overwrite the old one
    """
    if type(imgs) is str: imgs = [imgs]
    for img in imgs:
        if inverse:
            default('immath')
            immath(imagename = img, mode = 'evalexpr', expr = 'IM0', region = region, outfile = img+'-out')
            os.system('rm -r '+img)
            os.system('mv '+img+'-out '+img)
        else:
            ia.open(img)
            reg = rg.fromtextfile(filename=region, shape=ia.shape(), csys=ia.coordsys().torecord())
            ia.set(pixels=setTo, region=reg)
            ia.close()

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_blank(**params)
