#!/bin/bash

ms=$1
img=img/$2
mask=$3

awimager ms=$ms \
    image=$img \
    weight=briggs robust=0 \
    npix=4096 cellsize=7arcsec \
    padding=1. \
    data=CORRECTED_DATA \
    stokes=I \
    niter=100000 \
    operation=mfclark \
    oversample=5 \
    wmax=20000 \
    cyclefactor=1.5 \
    gain=0.1 \
    timewindow=300 \
    ApplyElement=0 \
    UVmin=0.1 \
    mask=$mask
