#!/bin/bash

ms=$1
img=$2
mask=$3

source /opt/cep/tools/citt/lofarinit.csh 

awimager clean verbose=1 output.imagename=$img data.ms=$ms clean.niter=100000 clean.maskimage=$mask image.npix=5000 image.cellsize=5arcsec image.nterms=2 weight.type=robust weight.robust=0
