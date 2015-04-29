#!/usr/bin/python
#casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_ft.py $ms $model
import os, sys

args = sys.argv[6:]
msfile = args[0]
mod = args[1]
wpp = int(args[2])

def update_freq_range(modelfile):
    """Allow a model file to cover a wide range of freq,
    larger than the one it was originally created for
    """
    ia.open(modelfile)
    csys = ia.coordsys()
    q = csys.increment()
    q['numeric'][3] = 1e10 # a random large value
    csys.setincrement(q)
    ia.setcoordsys(csys.torecord())
    ia.close()


tty = 1
if os.path.exists(mod+'.tt1'): tty = 2
if os.path.exists(mod+'.tt1') and os.path.exists(mod+'.tt2'): tty = 3

default('ft')
if tty == 1:
    update_freq_range(mod)
    ftw(vis=msfile, model=mod, nterms=1, usescratch=True, wprojplanes=wpp)
if tty == 2:
    update_freq_range(mod+'.tt0')
    update_freq_range(mod+'.tt1')
    ftw(vis=msfile, model=[mod+'.tt0', mod+'.tt1'], nterms=2, usescratch=True, wprojplanes=wpp)
if tty == 3:
    update_freq_range(mod+'.tt0')
    update_freq_range(mod+'.tt1')
    update_freq_range(mod+'.tt2')
    ftw(vis=msfile, model=[mod+'.tt0', mod+'.tt1', mod+'.tt2'], nterms=3, usescratch=True, wprojplanes=wpp)
