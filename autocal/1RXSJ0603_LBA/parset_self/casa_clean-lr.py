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

clean(vis=active_ms,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=25000,gain=0.1,psfmode="clark",\
        imagermode="csclean",multiscale=[],interactive=False,\
        mask=mask,imsize=[2500],cell=['15arcsec'],weighting="briggs",robust=0,usescratch=False,cyclefactor=2,cyclespeedup=-1,nterms=2,\
        uvtaper=True,outertaper='60arcsec')
