#!/usr/bin/python
# MPIP: casaargs = msfile
import sys

args = sys.argv[6:]
msfile = args[0]

ft(vis=msfile, model=['/cep1home/fdg/scripts/autocal/VirA_LBA/selfap10.model.tt0', '/cep1home/fdg/scripts/autocal/VirA_LBA/selfap10.model.tt1'], nterms=2, usescratch=True)
