#!/usr/bin/python
# casapy --nogui --log2term --nologger -c this_script.py

import os, sys

args = sys.argv[6:]

active_ms = args[0]
clean(vis=active_ms,imagename="img/"+active_ms.replace('.MS', '-wide'),mode="mfs",gridmode="widefield",wprojplanes=512,niter=5000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=[0,3,9,27],interactive=False,mask='/cep3home/fdg/scripts/autocal/1RXSJ0603_LBA/tooth_mask.crtf',imsize=[4096],cell=['5arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=2)
