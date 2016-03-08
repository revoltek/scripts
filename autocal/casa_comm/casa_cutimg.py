#!/usr/bin/python
# casapy --nogui --log2term --nologger -c  casa_blank.py picklefile
# cut an image to a specified region

import os, sys, pickle

def casa_cutimg(img = '', region = '', out = ''):
    """
    img = image to cut
    region = region to restrict image to
    out = output file name
    """
    if out == '': out = img+'-cut'
    default('immath')
    immath(imagename = img, mode = 'evalexpr', expr = 'IM0', region = region, outfile = out)

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_cutimg(**params)
