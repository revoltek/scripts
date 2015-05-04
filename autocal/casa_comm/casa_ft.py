#!/usr/bin/python
#casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_ft.py picklefile
import os, sys, pickle

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

def casa_ft(msfile='', model='', incr=False, wproj=None):
    """
    msfile = MS name
    model = model name, no tt* extension
    incr = bool, if to append to existing MODEL_DATA
    wproj = number of wprojplanes, or None for normal FFT
    """
    tty = 1
    if os.path.exists(model+'.tt1'): tty = 2
    if os.path.exists(model+'.tt1') and os.path.exists(model+'.tt2'): tty = 3

    if tty == 1:
        update_freq_range(model)
        if wproj != None:
            default('ftw')
            ftw(vis=msfile, model=model, nterms=1, usescratch=True, incremental=incr, wprojplanes=wproj)
        else:
            default('ft')
            ft(vis=msfile, model=model, nterms=1, usescratch=True, incremental=incr)
    if tty == 2:
        update_freq_range(model+'.tt0')
        update_freq_range(model+'.tt1')
        if wproj != None:
            default('ftw')
            ftw(vis=msfile, model=[model+'.tt0', model+'.tt1'], nterms=2, usescratch=True, incremental=incr, wprojplanes=wproj)
        else:
            default('ft')
            ft(vis=msfile, model=[model+'.tt0', model+'.tt1'], nterms=2, usescratch=True, incremental=incr)
    if tty == 3:
        update_freq_range(model+'.tt0')
        update_freq_range(model+'.tt1')
        update_freq_range(model+'.tt2')
        if wproj != None:
            default('ftw')
            ftw(vis=msfile, model=[model+'.tt0', model+'.tt1', model+'.tt2'], nterms=3, usescratch=True, incremental=incr, wprojplanes=wproj)
        else:
            default('ft')
            ft(vis=msfile, model=[model+'.tt0', model+'.tt1', model+'.tt2'], nterms=3, usescratch=True, incremental=incr)

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_ft(**params)
