#!/usr/bin/python
# concat SBs in an observations with gaps

prefix = 'cal'
sbs = [ '%03d' % int(x) for x in range(0,243)]

import os, sys, re
import glob
from lib_pipeline import *

s = Scheduler(qsub=False, max_threads=10, dry=False)

os.system('mkdir concat')

for sb in sbs:
    SBtoConcat = []

    for obs in sorted(glob.glob('L*')):
        # format L340894_SB078_uv.dppp.MS
        SBtoConcat.append(obs+'/'+obs+'_SB'+sb+'_uv.dppp.MS')

    newfilename = 'concat/'+prefix+'_SB'+sb+'.MS'
    s.add('concat_timehack.py -o '+newfilename+' '+' '.join(SBtoConcat)+'' , \
        log=sb+'.log', cmd_type='python')

s.run(check=True)
