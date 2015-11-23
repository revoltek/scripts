#!/usr/bin/python
# casapy --nogui --log2term --nologger -c  casa_blank.py picklefile
# blank all the pixel outside this region
# NOTE: overwrite the image! And the blanked image is the smallest around the given region.

import os, sys, pickle

def casa_fits2ms(imgs = [], del_fits = False):
    """
    img = image to convert
    """
    if type(imgs) is str: imgs = [imgs]
    for img in imgs:
        # define output name
        if img[-5:] == '.fits': imagename = img.replace('.fits','.ms')
        else: imagename = img+'.ms'

        default('importfits')
        importfits(fitsimage=img, imagename=imagename)

        if del_fits: os.system('rm '+img)

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_fits2ms(**params)
