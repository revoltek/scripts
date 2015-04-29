#!/usr/bin/python
# initial calibration of the calibrator, sol flag

skymodel = '/home/fdg/model/3C196-allfield.skymodel'
max_threads = 20

###################################################

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from lib_pipeline import *

set_logger()

#################################################
# Clean
logging.info('Cleaning...')
check_rm('*log')

##############################################
# Initial processing
logging.info('Fix beam table')
cmds=[]
for ms in sorted(glob.glob('*MS')):
    cmds.append('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms+' > '+ms+'_fixbeam.log 2>&1')
#thread_cmd(cmds, max_threads)

##############################################
# [PARALLEL] initial calibrator
logging.info('Calibrating')
cmds=[]
for ms in sorted(glob.glob('*MS')):
    cmds.append('calibrate-stand-alone -f '+ms+'  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_init/bbs-cal_field.parset '+skymodel+' > '+ms+'_cal.log 2>&1')
#thread_cmd(cmds, max_threads)
logging.warning('Bad runs:')
os.system('grep -L success *cal.log')

##############################################
# Clock check and flagging
check_rm('globaldb')
os.system('mkdir globaldb')

logging.info('Running LoSoTo...')
for i, ms in enumerate(sorted(glob.glob('*MS'))):
    num = re.findall(r'\d+', ms)[-1]
    logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')

check_rm('plot')
os.system('mkdir plot')
check_rm('cal.h5')
os.system('H5parm_importer.py -v cal.h5 globaldb')
os.system('losoto.py -v cal.h5 /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_init/losoto.parset')
os.system('H5parm_exporter.py -v cal.h5 globaldb')

logging.info("Done.")
