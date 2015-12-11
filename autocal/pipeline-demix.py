#!/usr/bin/python
# demix of a set of SBs from a given dir, output is on local dir

parset_dir = '/home/fdg/scripts/autocal/VirA_LBA/parset_demix/'
#origmss_dir = '/data/scratch/fdg/virgoLBAis/tgts-bkp/'
origmss_dir = '/data/scratch/fdg/virgoLBAis/demix-test/orig/'

###################################################

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from lib_pipeline import *

set_logger()
check_rm('logs')
s = Scheduler(dry=False)
mss = sorted(glob.glob(origmss_dir+'/*MS'))

##############################################
# Initial processing (2/2013->2/2014)
#logging.warning('Fix beam table')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)

##############################################
# Demix
logging.info('Demixing...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP_demix.parset msin='+ms+' msout='+os.path.basename(ms)+' demixer.instrumentmodel='+os.path.basename(ms)+'/instrument_demix', log=os.path.basename(ms)+'_demix.log', cmd_type='NDPPP')
s.run(check=True)

logging.info("Done.")
