#!/usr/bin/python

# rename all SB in current dir scaline their SB number by this factor
f = -244

import os, sys, re
import glob

for filename in sorted(glob.glob('*MS')):
    SBnum = str(re.findall(r'\d+', filename)[-1])
    newfilename = filename.replace('SB'+SBnum, 'SB%03d' % (int(SBnum)+f) )
    print 'mv '+filename+' '+newfilename
    os.system('mv '+filename+' '+newfilename)
