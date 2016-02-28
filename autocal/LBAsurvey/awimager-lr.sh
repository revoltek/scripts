#!/bin/bash

ms=$1
img=$2
mask=$3

awimager ms=$ms \
    image=$img \
    weight=briggs robust=0 \
    npix=700 cellsize=30arcsec \
    padding=1. \
    data=CORRECTED_DATA stokes=I \
    niter=25000 \
    operation=mfclark \
    oversample=5 \
    wmax=20000 \
    cyclefactor=1.5 \
    gain=0.1 \
    timewindow=300 \
    ApplyElement=0 \
    UVmin=0.1 \
    mask=$mask \
    filter=120arcsec,120arcsec,0deg
