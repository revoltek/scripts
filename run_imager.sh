#!/bin/bash 

source ~tasse/jaws26/init.sh
awimager ms=${1} image=${2}  weight=briggs npix=8192 robust=0. cellsize=2.5arcsec data=CORRECTED_DATA padding=1.3 niter=75000 stokes=I threshold=0.0 operation=mfclark timewindow=300 wmax=200000 cyclefactor=1.5 gain=0.1 PBCut=0.03 UseLIG=1 UseWSplit=1 ApplyBeamCode=1 ApplyElement=0 TWElement=20 MakeDirtyCorr=0 SingleGridMode=1 FindNWplanes=1 ChanBlockSize=1
