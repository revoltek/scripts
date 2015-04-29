#!/usr/bin/env python
# make the calib of the calibrator in circular and then run losoto
# to be run in the directory with all the calibrators SBs
# TODO: once circ is in NDPPP the first BBS can be ported to the fast stefcal version

skymodel = '/home/fdg/model/3C295-allfield.skymodel'
skymodel_cal = '/home/fdg/model/3C295_dual.skymodel'

##############################################################

import sys, os, glob, re
from lofar import bdsm
import numpy as np
import lsmtool
import pyrap.tables as pt
from lib_pipeline import *

set_logger()

#################################################
# Clean
logging.info('Cleaning...')
check_rm('*log')

#################################################
# [PARALLEL] Calibrate on the field - SB.MS:DATA -> SB.MS:DATA_SUB_BC (data, field subtracted, beam corrected, linear)
logging.info('Subtract field and correct beam...')
cmds = []
for ms in glob.glob('*MS'):
    cmds.append('calibrate-stand-alone -f '+ms+' /home/fdg/scripts/autocal/VirA_LBA/parsets_init/bbs_sub-beamcorr.parset /home/fdg/model/3C295-allfield.skymodel > '+ms+'-sub-beamcorr.log 2>&1')
thread_cmd(cmds)
logging.warning('Bad runs:')
os.system('grep -L "successfully" *-sub-beamcorr.log')

#################################################
# [PARALLEL] Transform subtracted data to circular - SB.MS:DATA_SUB_BC -> SB.MS:DATA_SUB_CIRC (data, field subtracted, beam corrected, circular)
logging.info('Convert to circular...')
cmds = []
for ms in glob.glob('*MS'):
    cmds.append('/home/fdg/scripts/mslin2circ.py -i '+ms+':DATA_SUB_BC -o '+ms+':DATA_SUB_CIRC > '+ms+'-lin2circ.log 2>&1')
thread_cmd(cmds)

################################################
# [PARALLEL] Calibrate on the calibratora - SB.MS:DATA_SUB_CIRC (no correction)
logging.info('Calibrating the calibrators...')
cmds=[]
for ms in glob.glob('*.MS'):
    check_rm(ms+'/skymodel')
    cmds.append('makesourcedb in='+skymodel_cal+' out='+ms+'/skymodel format="<"  > '+ms+'-makesourcedb.log 2>&1')
thread_cmd(cmds)
cmds=[]
for ms in glob.glob('*.MS'):
    cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_init/NDPPP-sol.parset msin='+ms+' sol.parmdb='+ms+'/instrument sol.sourcedb='+ms+'/skymodel > '+ms+'-sol.log 2>&1')
thread_cmd(cmds)
logging.warning('Bad runs:')
os.system('grep -L "Total NDPPP time" *-cal.log')

##############################################
# Clock check and flagging solutions
check_rm('globaldb')
os.system('mkdir globaldb')

logging.info('Running LoSoTo...')
for i, ms in enumerate(sorted(glob.glob('*MS'))):
    logging.debug('Copy instrument of '+ms)
    num = re.findall(r'\d+', ms)[-1]
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
    print 'cp -r '+ms+'/instrument globaldb/instrument-'+str(num)
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')

check_rm('plot')
os.system('mkdir plot')
check_rm('cal.h5')
os.system('H5parm_importer.py -v cal.h5 globaldb')
os.system('losoto.py -v cal.h5 /home/fdg/scripts/autocal/VirA_LBA/parsets_init/losoto.parset')
os.system('H5parm_exporter.py -v cal.h5 globaldb')

logging.info("Done.")
