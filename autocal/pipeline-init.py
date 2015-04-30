#!/usr/bin/python
# initial calibration of the calibrator, sol flag

#skymodel = '/home/fdg/model/3C196-allfield.skymodel' # tooth
#parset_dir = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_init' # tooth
#skymodel = '/home/fdg/model/3C295-allfield.skymodel' # virgo LBA
#parset_dir = '/home/fdg/scripts/autocal/VirA_LBA/parset_init' # virgo LBA
skymodel = '/home/fdg/model/3C196-allfield.skymodel' # perseus LBA
parset_dir = '/home/fdg/scripts/autocal/PerA_LBA/parset_init' # perseus LBA

###################################################

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from lib_pipeline import *

set_logger()
s = Scheduler(qsub=False, max_threads=20, dry=False)
mss = sorted(glob.glob('*MS'))

#################################################
# Clean
logging.info('Cleaning...')
check_rm('*log')

##############################################
# Initial processing
logging.info('Fix beam table')
for ms in mss:
    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
s.run(check=False)

##############################################
# [PARALLEL] initial calibrator
logging.info('Calibrating with skymodel: '+skymodel)
for ms in mss:
    s.add('calibrate-stand-alone -f '+ms+' '+parset_dir+'/bbs-cal_field.parset '+skymodel, log=ms+'_cal.log', cmd_type='BBS')
s.run(check=True)

##############################################
# Clock check and flagging
check_rm('globaldb')
os.system('mkdir globaldb')

logging.info('Running LoSoTo...')
for i, ms in enumerate(mss):
    num = re.findall(r'\d+', ms)[-1]
    logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')

check_rm('plot')
os.makedirs('plot')
check_rm('cal.h5')
s.add('H5parm_importer.py -v cal.h5 globaldb', log='losoto.log', cmd_type='python')
s.run(check=False)
s.add('losoto.py -v cal.h5 '+parset_dir+'/losoto.parset', log='losoto.log', log_append=True, cmd_type='python')
s.run(check=False)
s.add('H5parm_exporter.py -v cal.h5 globaldb', log='losoto.log', log_append=True, cmd_type='python')
s.run(check=True)

logging.info("Done.")
