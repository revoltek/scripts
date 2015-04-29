#!/usr/bin/python
# run minimal tfcrop
import sys

args = sys.argv[6:]
active_ms = args[0]

default('flagdata')
statsflags = flagdata(vis=active_ms, mode='summary', spwchan=False, spwcorr=False, basecnt=False, action='calculate', flagbackup=False, savepars=False, async=False)
clearstat()
print "INFO: Initial flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"

default('flagdata')
flagdata(vis=active_ms, mode='tfcrop', datacolumn='corrected', action='apply')

default('flagdata')
statsflags = flagdata(vis=active_ms, mode='summary', spwchan=False, spwcorr=False, basecnt=False, action='calculate', flagbackup=False, savepars=False, async=False)
clearstat()
print "INFO: Initial flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%"
