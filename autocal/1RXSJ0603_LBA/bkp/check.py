#!/usr/bin/env python

import sys, os, glob, re
from pipeline_lib import *

set_logger()

################################################
# Integrity check
for cal in sorted(glob.glob('cals/*MS')):
    if not os.path.exists(cal.replace('cals','tgts').replace('3C196','tooth')):
        logging.warning("Cannot find " + cal.replace('cals','tgts').replace('3C196','tooth')+' remove '+cal)
for tgt in sorted(glob.glob('tgts/*MS')):
    if not os.path.exists(tgt.replace('tgts','cals').replace('tooth','3C196')):
        logging.warning("Cannot find " + tgt.replace('tgts','cals').replace('tooth','3C196')+' remove '+tgt)

logging.info("Done.")
