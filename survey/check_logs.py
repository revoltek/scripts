#!/usr/bin/env python

import os, sys, glob

for tgt in sorted(glob.glob('*-pointings/p*')):
    os.chdir('/home/fdg/data/LBAsurvey/%s' % tgt)
    print ("-------- Working on: %s" % tgt)
    cmd = 'cat logs/ddfacetM-c1.log | grep "iter" | tail -n1'
    os.system(cmd)
