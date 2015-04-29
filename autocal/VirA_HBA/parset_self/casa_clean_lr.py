#!/usr/bin/python
# MPIP: casaargs = msfile
import sys

args = sys.argv[6:]
i = args[-1]
args.remove(i)
active_ms = args[-1]

default('clean')
clean(vis=active_ms,imagename='img/lr-'+str(i),mode="mfs",gridmode="widefield",wprojplanes=1024,niter=50,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask='',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=8,cyclespeedup=-1,nterms=3,uvtaper=T,outertaper='60arcsec')
default('clean')
scales=[0, 5, 10, 20]
clean(vis=active_ms,imagename='img/lr-'+str(i),mode="mfs",gridmode="widefield",wprojplanes=1024,niter=150,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask='',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=8,cyclespeedup=-1,nterms=3,uvtaper=T,outertaper='60arcsec')
default('clean')
scales=[0, 5, 10, 20, 40, 80, 160]
clean(vis=active_ms,imagename='img/lr-'+str(i),mode="mfs",gridmode="widefield",wprojplanes=1024,niter=5000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask='',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=3,uvtaper=T,outertaper='60arcsec')


ft(vis=active_ms, model=['img/lr-'+str(i)+'.model.tt0', 'img/lr-'+str(i)+'.model.tt1','img/lr-'+str(i)+'.model.tt2'], nterms=3, usescratch=True)
uvsub(vis=active_ms)
clean(vis=active_ms,imagename='img/lr-uvsub-'+str(i),antenna='',mode="mfs",gridmode="widefield",wprojplanes=1024,niter=500,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,imsize=[1024],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=3,uvtaper=T,outertaper='120arcsec')
