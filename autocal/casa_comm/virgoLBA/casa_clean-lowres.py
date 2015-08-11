#!/usr/bin/python
# casapy --nogui --log2term --nologger -c this_script.py picklefile

import sys, pickle

def casa_clean(msfile='', imagename='', imtype='normal', mask=''):
    """
    imtype: 'lr' for low resolution, large version of the image (larger FoV to get possible new extended emission)
            'wide' for wide FoV image
            'widemasked' for wide FoV image but with less cycles (useful when given a mask)
            'normal' for full res M87 image
            'dirty' quick dirty image for testing purposes
    """

    if imtype == 'normal': 
        default('clean')
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=50,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask='/home/fdg/scripts/autocal/VirA_LBA/m87.crtf',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=8,cyclespeedup=-1,nterms=3)
        default('clean')
        scales=[0, 5, 10, 20, 40]
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=5000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask=mask,imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=3)

        clean(vis=msfile,imagename=imagename+'-r-1',mode="mfs",gridmode="widefield",wprojplanes=512,niter=50,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask='/home/fdg/scripts/autocal/VirA_LBA/m87.crtf',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=-1,cyclefactor=8,cyclespeedup=-1,nterms=3)
        default('clean')
        scales=[0, 5, 10, 20, 40]
        clean(vis=msfile,imagename=imagename+'-r-1',mode="mfs",gridmode="widefield",wprojplanes=512,niter=5000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask=mask,imsize=[512],cell=['10arcsec'],weighting="briggs",robust=-1,cyclefactor=5,cyclespeedup=-1,nterms=3)

        clean(vis=msfile,imagename=imagename+'-r1',mode="mfs",gridmode="widefield",wprojplanes=512,niter=50,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask='/home/fdg/scripts/autocal/VirA_LBA/m87.crtf',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=1,cyclefactor=8,cyclespeedup=-1,nterms=3)
        default('clean')
        scales=[0, 5, 10, 20, 40]
        clean(vis=msfile,imagename=imagename+'-r1',mode="mfs",gridmode="widefield",wprojplanes=512,niter=5000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask=mask,imsize=[512],cell=['10arcsec'],weighting="briggs",robust=1,cyclefactor=5,cyclespeedup=-1,nterms=3)

    elif imtype == 'wide':
        default('clean')
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=50,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask='',imsize=[1024],cell=['15arcsec'],weighting="briggs",robust=0,cyclefactor=8,cyclespeedup=-1,nterms=3)
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=5000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=[0,3,6,12,24,48],interactive=False,mask=mask,imsize=[1024],cell=['15arcsec'],weighting="briggs",robust=0,cyclefactor=8,cyclespeedup=-1,nterms=3)

    else:
        print "Wrong imtype."

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_clean(**params)
