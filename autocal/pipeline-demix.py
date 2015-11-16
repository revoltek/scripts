#!/usr/bin/python
# demix of a set of SBs from a given dir, output is on local dir

skymodel = '/home/fdg/scripts/model/3C295-allfield.skymodel'
parset_dir = '/home/fdg/scripts/autocal/LBAmode/parset_demix/'
origmss_dir = '/lofar2/fdg/virgoHBA/...'

###################################################

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from lib_pipeline import *

set_logger()
s = Scheduler(dry=False)
mss = sorted(glob.glob(origmss_dir+'/*MS'))

#################################################
# Clean
logging.info('Cleaning...')
check_rm('*log')

##############################################
# Initial processing (2/2013->2/2014)
#logging.warning('Fix beam table')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)

##############################################
# Demix
logging.info('Calibrating with skymodel: '+skymodel)
for ms in mss:
    s.add('NDPPP msin='+ms+' msout='+os.path.basename(ms)+' '+parset_dir+'/NDPPP_demix.parset '+skymodel, log=ms+'_demix.log', cmd_type='NDPPP')
s.run(check=True)

logging.info("Done.")
