#!/bin/bash

#singularity build ~/storage/MeerKATpol.simg docker://tpasini/pol_meerkat:latest
singularity shell --pid --writable-tmpfs --cleanenv -B/homes/fdg,/local/work/fdg,/iranet/groups/ulu/fdg/,/iranet/groups/lofar/containers/ ~/storage/MeerKATpol.simg
