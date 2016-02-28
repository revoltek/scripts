#!/usr/bin/python
# casapy --nogui --log2term --nologger -c this_script.py picklefile

import os, sys, pickle

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
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=300,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask='/home/fdg/scripts/autocal/VirgoLBA/m87.crtf',imsize=[1024],cell=['3arcsec'],weighting="briggs",robust=-0.5,cyclefactor=8,cyclespeedup=-1,nterms=3)
        os.system('cp -r '+imagename+'.residual.tt0 '+imagename+'-bkp1')

        default('clean')
        scales=[0, 5, 10, 20]
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=250,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask='/home/fdg/scripts/autocal/VirgoLBA/m87.crtf',imsize=[1024],cell=['3arcsec'],weighting="briggs",robust=-0.5,cyclefactor=8,cyclespeedup=-1,nterms=3)
        os.system('cp -r '+imagename+'.residual.tt0 '+imagename+'-bkp2')

        default('clean')
        scales=[0, 5, 10, 20, 40, 80, 160, 320]
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=15000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask='/home/fdg/scripts/autocal/VirgoLBA/m87.crtf',imsize=[512],cell=['4arcsec'],weighting="briggs",robust=-0.5,cyclefactor=5,cyclespeedup=-1,nterms=3)

    elif imtype == 'dirty':
        default('clean')
        clean(vis=msfile,imagename=imagename,mode="mfs",niter=1,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,imsize=[512],cell=['3arcsec'],weighting="briggs",robust=-0.5,cyclefactor=8,cyclespeedup=-1,nterms=1)

    elif imtype == 'lr':
        default('clean')
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=50,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask='/home/fdg/scripts/autocal/VirgoLBA/m87.crtf',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=8,cyclespeedup=-1,nterms=3,uvtaper=T,outertaper='60arcsec')

        default('clean')
        scales=[0, 5, 10, 20, 40]
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=5000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask='/home/fdg/scripts/autocal/VirgoLBA/m87.crtf',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=3,uvtaper=T,outertaper='60arcsec')

    elif imtype == 'wide':
        default('clean')
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=10000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,imsize=[2500],cell=['10arcsec'],weighting="briggs",robust=0,nterms=2,uvtaper=True,outertaper='30arcsec', mask=mask)

    elif imtype == 'widemasked':
        default('clean')
        clean(vis=msfile,imagename=imagename,mode="mfs",gridmode="widefield",wprojplanes=512,niter=5000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,imsize=[2500],cell=['10arcsec'],weighting="briggs",robust=0,nterms=2,uvtaper=True,outertaper='30arcsec', mask=mask)

    else:
        print "Wrong imtype."

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_clean(**params)
