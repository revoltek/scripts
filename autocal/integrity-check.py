#!/usr/bin/env python
# check if each the calibrator SB has the correspective tgt SB and viceversa
# integrity-check.py dir1 dir2

import sys, os, glob, re
from lib_pipeline import *

set_logger()

dir1 = sys.argv[1]
dir2 = sys.argv[2]

calname = glob.glob(dir1+'/*MS')[0].split('/')[1].split('_')[0]
print 'Dir1: ', calname
tgtname = glob.glob(dir2+'/*MS')[0].split('/')[1].split('_')[0]
print 'Dir2: ', tgtname

for cal in sorted(glob.glob(dir1+'/*MS')):
    if not os.path.exists(cal.replace(dir1,dir2).replace(calname,tgtname)):
        logging.warning("Cannot find " + cal.replace(dir1,dir2).replace(calname,tgtname)+' remove '+cal)
for tgt in sorted(glob.glob(dir2+'/*MS')):
    if not os.path.exists(tgt.replace(dir2,dir1).replace(tgtname,calname)):
        logging.warning("Cannot find " + tgt.replace(dir2,dir1).replace(tgtname,calname)+' remove '+tgt)

logging.info("Done.")
