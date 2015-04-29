#!/usr/bin/python
# MPIP: casaargs = msfile
import sys

args = sys.argv[6:]
i = args[-1]
args.remove(i)
active_ms = args

default('clean')
clean(vis=active_ms,imagename='img/normal_clean-'+str(i),mode="mfs",gridmode="widefield",wprojplanes=512,niter=250,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask='/home/fdg/scripts/autocal/VirA_LBA/m87.crtf',imsize=[512],cell=['3arcsec'],weighting="briggs",robust=-1,cyclefactor=8,cyclespeedup=-1,nterms=2)

default('clean')
scales=[0, 5, 10, 20]
clean(vis=active_ms,imagename='img/normal_clean-'+str(i),mode="mfs",gridmode="widefield",wprojplanes=512,niter=250,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask='/home/fdg/scripts/autocal/VirA_LBA/m87.crtf',imsize=[512],cell=['3arcsec'],weighting="briggs",robust=-1,cyclefactor=8,cyclespeedup=-1,nterms=2)

default('clean')
scales=[0, 5, 10, 20, 40, 80, 160, 320]
clean(vis=active_ms,imagename='img/normal_clean-'+str(i),mode="mfs",gridmode="widefield",wprojplanes=512,niter=10000,gain=0.1,psfmode="clark",imagermode="csclean",multiscale=scales,interactive=False,mask='/home/fdg/scripts/autocal/VirA_LBA/m87.crtf',imsize=[512],cell=['3arcsec'],weighting="briggs",robust=-1,cyclefactor=5,cyclespeedup=-1,nterms=2)
