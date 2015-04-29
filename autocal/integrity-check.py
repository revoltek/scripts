#!/usr/bin/env python
# check if each the calibrator SB has the correspective tgt SB and viceversa

import sys, os, glob, re
from lib_pipeline import *

set_logger()

calname = glob.glob('cals/*MS')[0].split('/')[1].split('_')[0]
print 'Calname', calname
tgtname = glob.glob('tgts/*MS')[0].split('/')[1].split('_')[0]
print 'Tgtname', tgtname

for cal in sorted(glob.glob('cals/*MS')):
    if not os.path.exists(cal.replace('cals','tgts').replace(calname,tgtname)):
        logging.warning("Cannot find " + cal.replace('cals','tgts').replace(calname,tgtname)+' remove '+cal)
for tgt in sorted(glob.glob('tgts/*MS')):
    if not os.path.exists(tgt.replace('tgts','cals').replace(tgtname,calname)):
        logging.warning("Cannot find " + tgt.replace('tgts','cals').replace(tgtname,calname)+' remove '+tgt)

logging.info("Done.")
