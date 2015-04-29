#!/usr/bin/python
# MPIP: casaargs = msfile
import sys

args = sys.argv[6:]
msfile = args[0]
mod = args[1]
tty = int(args[2])

default('ft')
if tty == 1:
    ft(vis=msfile, model=mod, nterms=1, usescratch=True)
if tty == 2:
    ft(vis=msfile, model=[mod+'.tt0', mod+'.tt1'], nterms=2, usescratch=True)
if tty == 3:
    ft(vis=msfile, model=[mod+'.tt0', mod+'.tt1', mod+'tt2'], nterms=3, usescratch=True)
