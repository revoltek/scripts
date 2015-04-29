#!/usr/bin/python
# casapy --nogui --log2term --nologger -c this_script.py

import os, sys

args = sys.argv[6:]

active_ms = args[0]
imagename = args[1]
if len(args) == 3:
    mask = args[2]
else:
    mask = ''

clean(vis=active_ms,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=128,niter=5000,gain=0.1,psfmode="clark",\
        imagermode="csclean",multiscale=[0,3,9,18],interactive=False,\
        mask=mask,imsize=[256],cell=['4arcsec'],weighting="briggs",robust=0,usescratch=False,cyclefactor=5,cyclespeedup=-1,nterms=2)
