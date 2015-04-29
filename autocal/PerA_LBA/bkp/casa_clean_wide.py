#!/usr/bin/python
# casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_clean_wide.py ms imgname mask
import os, sys

args = sys.argv[6:]
active_ms = args[0]
i = args[1]
if len(args) == 3: usermask = args[2]
else: usermask = ''

default('clean')
clean(vis=active_ms,imagename='img/widefield-'+str(i),mode="mfs",gridmode="widefield",wprojplanes=512,niter=5000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask=usermask,imsize=[4096],cell=['5arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=2)
