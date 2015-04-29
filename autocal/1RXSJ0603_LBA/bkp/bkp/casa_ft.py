#!/usr/bin/python
# MPIP: casaargs = msfile
import sys

args = sys.argv[6:]
msfile = args[0]
mod = args[1]
tty = int(args[2])

default('ftw')
if tty == 1:
    ftw(vis=msfile, model=mod, nterms=1, wprojplanes=512, usescratch=True)
if tty == 2:
    ftw(vis=msfile, model=[mod+'.tt0', mod+'.tt1'], nterms=2, wprojplanes=512, usescratch=True)
if tty == 3:
    ftw(vis=msfile, model=[mod+'.tt0', mod+'.tt1', mod+'tt2'], nterms=3, wprojplanes=512, usescratch=True)
