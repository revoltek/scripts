#!/usr/bin/python
# casapy --nogui --log2term --nologger -c this_script.py picklefile

import os, sys, pickle

def casa_clean(msfile='', imagename='', imsize=1024, mask='', niter=1000, multiscale=[], wproj=512, cell='3arcsec', uvtaper=False, outertaper=[]):
    default('clean')
    clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=wproj,niter=niter,gain=0.1,psfmode="clark",\
        imagermode="csclean",multiscale=multiscale,interactive=False, mask=mask,imsize=[int(imsize)],cell=cell,weighting="briggs",robust=0,\
        usescratch=False,cyclefactor=5,cyclespeedup=-1,nterms=2,uvtaper=uvtaper,outertaper=outertaper)

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_clean(**params)
