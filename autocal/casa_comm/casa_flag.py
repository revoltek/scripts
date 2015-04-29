#!/usr/bin/python
# casapy --nogui --log2term --nologger -c this_script.py picklefile

import sys, pickle

def casa_flag(msfile=''):
    """
    Flag the CORRECTED_DATA of this msfile
    """
    default('flagdata')
    statsflags = flagdata(vis=msfile, mode='summary', spwchan=False, spwcorr=False, basecnt=False, action='calculate', flagbackup=False, savepars=False, async=False)
    clearstat()
    print "INFO: Initial flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"

    default('flagdata')
    flagdata(vis=msfile, mode='tfcrop', datacolumn='corrected', action='apply')

    default('flagdata')
    statsflags = flagdata(vis=msfile, mode='summary', spwchan=False, spwcorr=False, basecnt=False, action='calculate', flagbackup=False, savepars=False, async=False)
    clearstat()
    print "INFO: Final flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"

params = pickle.load( open( sys.argv[6], "rb" ) )
os.system('rm '+sys.argv[6])
casa_flag(**params)
