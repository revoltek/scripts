#!/usr/bin/python
# casapy --nogui --log2term --nologger -c this_script.py picklefile

import sys, pickle

def casa_clean(msfile='', imagename='', imtype='normal', mask=''):
    """
    imtype: 'lr' for low resolution, large version of the image (larger FoV to get possible new extended emission)
            'wide' for wide FoV image
            'normal' for full res M87 image
    """

    if imtype == 'normal': 
        default('clean')
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=1000,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[0,3,9,27],interactive=False,mask=mask,imsize=[1024],cell=['5arcsec'],weighting="briggs",robust=0.,cyclefactor=8,cyclespeedup=-1,nterms=2)

    elif imtype == 'lr':
        default('clean')
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=50,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask=mask,imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=8,cyclespeedup=-1,nterms=2, antenna='CS*') #uvtaper=T,outertaper='60arcsec')

    elif imtype == 'wide':
        default('clean')
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=500,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,imsize=[1024],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=2,uvtaper=T,outertaper='120arcsec', mask=mask)

    else:
        print "Wrong imtype."

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_clean(**params)
