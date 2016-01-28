#!/usr/bin/python
# demix of a set of SBs from a given dir, output is in the local dir

#parset_dir = '/home/fdg/scripts/autocal/VirgoLBA/parset_demix/'
parset_dir = '/home/fdg/scripts/autocal/CygLBA/parset_demix/'
origmss_dir = '/data/scratch/fdg/CygnusLBAis/tgts1-full/'

###################################################

import sys, os, glob
import numpy as np
from lib_pipeline import *

set_logger()
check_rm('logs')
s = Scheduler(dry=False, max_threads = 5) # set here max number of threads here
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
    if os.path.exists(os.path.basename(ms)): continue
    s.add('NDPPP '+parset_dir+'/NDPPP_demix.parset msin='+ms+' msout='+os.path.basename(ms)+' demixer.instrumentmodel='+os.path.basename(ms)+'/instrument_demix', log=os.path.basename(ms)+'_demix.log', cmd_type='NDPPP')
s.run(check=True)

logging.info("Done.")
