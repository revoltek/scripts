#!/usr/bin/python
# MPIP: casaargs = msfile
import sys

args = sys.argv[6:]
msfile = args[0]

default('uvsub')
uvsub(vis=msfile)
