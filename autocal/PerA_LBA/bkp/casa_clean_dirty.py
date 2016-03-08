#!/usr/bin/python
# MPIP: casaargs = msfile
import sys,os

args = sys.argv[6:]
i = args[-1]
args.remove(i)
active_ms = args[-1]


clean(vis=active_ms,imagename='img/'+os.path.basename(active_ms)+'_'+str(i),mode="mfs",niter=1,gain=0.05,psfmode="clark",imagermode="csclean",multiscale=[],interactive=False,mask='/cep3home/fdg/scripts/autocal/VirA_LBA/m87.crtf',imsize=[512],cell=['3arcsec'],weighting="briggs",robust=0,cyclefactor=10,cyclespeedup=-1,nterms=1)
