#!/usr/bin/python
# initial calibration of the calibrator in circular, sol flag + effects separation
# also correct for beam(centre) + all phase/amp effects and subtract the calibrator

#skymodel = '/home/fdg/scripts/model/3C196-allfield.skymodel' # tooth LBA
#sourcedb = '/home/fdg/scripts/model/3C196-allfield.skydb' # tooth LBA
skymodel = '/home/fdg/scripts/model/3C295-allfield.skymodel' # tooth LBA
sourcedb = '/home/fdg/scripts/model/3C295-allfield.skydb' # tooth LBA
parset_dir = '/home/fdg/scripts/autocal/parset_cal'

###################################################

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from lib_pipeline import *

set_logger()
check_rm('logs')
s = Scheduler(dry=False)
mss = sorted(glob.glob('*MS'))

check_rm('globaldb*')
os.system('mkdir globaldb')

nchan = find_nchan(mss[0])
logging.debug('Channel in the MS: '+str(nchan)+'.')

##############################################
# Initial processing (2/2013->2/2014)
#logging.warning('Fix beam table...')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)

##############################################
# Beam correction DATA -> CORRECTED_DATA
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
s.run(check=True)

##############################################
# Split channels
logging.info('Splitting channels...')
for chan in xrange(nchan):
    logging.debug('Channel: '+str(chan))
    for ms in mss:
        msout = ms.replace('.MS','-chan'+str(chan)+'.MS')
        check_rm(msout)
        s.add('NDPPP '+parset_dir+'/NDPPP-split.parset msin='+ms+' msin.startchan='+str(chan)+' msout='+msout, log=ms+'-'+str(chan)+'_split.log', cmd_type='NDPPP')
    s.run(check=True)
    
mss = sorted(glob.glob('*-chan*.MS'))

##############################################
# Convert to circular DATA -> CIRC_DATA
logging.info('Converting to circular...')
for ms in mss:
    s.add('mslin2circ.py -i '+ms+':DATA -o '+ms+':CIRC_DATA', log=ms+'_circ2lin.log', cmd_type='python')
s.run(check=True)

##############################################
# Avg data CIRC_DATA -> SMOOTHED_DATA (BL-based smoothing)
# NOTE: the WEIGHTED_COLUMN is now smoothed in this dataset, a backup is in WEIGHTED_COLUMN_ORIG
logging.info('BL-averaging...')
for ms in mss:
    s.add('BLavg.py -r -w -i CIRC_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python')
s.run(check=True)

###############################################
# Initial calibrator
# Solve cal_SB.MS:SMOOTHED_DATA
logging.info('Calibrating with skymodel: '+sourcedb+'...')
for ms in mss:
    check_rm(ms+'/instrument')
    s.add('NDPPP '+parset_dir+'/NDPPP-cal.parset msin='+ms+' cal.parmdb='+ms+'/instrument cal.sourcedb='+sourcedb, log=ms+'_cal.log', cmd_type='NDPPP')
s.run(check=True)

##################################################
# Correct/subtract
# CIRC_DATA -> CORRECTED_DATA
logging.info('Correct and subtract...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-corsub.parset msin='+ms+' cor.parmdb='+ms+'/instrument sub.applycal.parmdb='+ms+'/instrument sub.sourcedb='+sourcedb, log=ms+'_corsub.log', cmd_type='NDPPP')
s.run(check=True)

for i, ms in enumerate(mss):
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')

    num = re.findall(r'\d+', ms)[-2]
    chan = re.findall(r'\d+', ms)[-1]
    logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(num)+'-'+str(chan))
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num)+'-'+str(chan))

####################################################

s.add('NDPPP '+parset_dir+'/NDPPP-concat.parset', log='concat.log', cmd_type='NDPPP')
s.run(check=True)

##############################################
# Clock/TEC check and flagging
logging.info('Running LoSoTo...')
check_rm('plots')
os.makedirs('plots')
check_rm('cal.h5')
s.add('H5parm_importer.py -v cal.h5 globaldb', log='losoto.log', cmd_type='python')
s.run(check=True)
#s.add('losoto -v cal.h5 '+parset_dir+'/losoto-flag.parset', log='losoto-flag.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
s.add('losoto -v cal.h5 '+parset_dir+'/losoto-amp.parset', log='losoto-amp.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('losoto -v cal.h5 '+parset_dir+'/losoto-ph.parset', log='losoto-ph.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)

logging.info("Done.")
