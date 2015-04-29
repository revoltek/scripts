#!/usr/bin/python
# MPIP: casaargs = msfile
import sys

args = sys.argv[6:]
i = args[-1]
args.remove(i)
active_ms = args

clean(vis=active_ms,imagename='img/low-res-'+str(i),mode="mfs",gridmode="widefield",wprojplanes=128,niter=150,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[0],interactive=False,mask='/cep3home/fdg/scripts/autocal/VirA_LBA/m87.crtf',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=2,uvtaper=T,outertaper='60arcsec')
clean(vis=active_ms,imagename='img/low-res-'+str(i),mode="mfs",gridmode="widefield",wprojplanes=128,niter=20000,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[0,5,15,30],interactive=False,mask='/cep3home/fdg/scripts/autocal/VirA_LBA/m87.crtf',imsize=[512],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=2,uvtaper=T,outertaper='60arcsec')

ft(vis=active_ms, model=['img/low-res-'+str(i)+'.model.tt0', 'img/low-res-'+str(i)+'.model.tt1'], nterms=2, usescratch=True)
uvsub(vis=active_ms)
clean(vis=active_ms,imagename='img/low-res-uvsub-'+str(i),antenna='CS*',mode="mfs",gridmode="widefield",wprojplanes=128,niter=500,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[0],interactive=False,imsize=[1024],cell=['10arcsec'],weighting="briggs",robust=0,cyclefactor=5,cyclespeedup=-1,nterms=2,uvtaper=T,outertaper='100arcsec')
