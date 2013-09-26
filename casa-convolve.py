#!/usr/bin/python
import sys
import os
out = sys.stdout
args = sys.argv[5:]
img = args[0]
res = args[1]

# CONVOLVE
imsmooth(imagename=img, kernel='gauss', major=res+'arcsec', minor=res+'arcsec', pa='0deg', targetres=True,outfile=img+'-convolved'+res)
