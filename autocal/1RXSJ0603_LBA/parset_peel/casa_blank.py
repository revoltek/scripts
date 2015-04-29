#!/usr/bin/python
# casapy --nogui --log2term --nologger -c  casa_blank.py picklefile
# blank all the pixel outside this region
# NOTE: overwrite the image! And the blanked image is the smallest around the given region.

import os, sys, pickle

def casa_blank(img = [], region = ''):
    """
    img = image to blank
    region = blank all pixel outside this region
    out: the new, blanked image will overwrite the old one
    """
    default('immath')
    immath(imagename = img, mode = 'evalexpr', expr = 'IM0', region = region, outfile = img+'-out')
    os.system('rm -r '+img)
    os.system('mv '+img+'-out '+img)

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_blank(**params)
