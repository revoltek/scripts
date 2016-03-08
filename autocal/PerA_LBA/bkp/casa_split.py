#!/usr/bin/python
# MPIP: casaargs = msfile
import sys

args = sys.argv[6:]
msfile = args[0]

split(vis=msfile, outputvis=msfile.replace('.MS','.split.MS'))
