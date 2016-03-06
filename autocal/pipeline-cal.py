#!/usr/bin/python
# initial calibration of the calibrator in circular, sol flag + effects separation

skymodel = '/home/fdg/scripts/model/3C196-allfield.skymodel'
sourcedb = '/home/fdg/scripts/model/3C196-allfield.skydb'
patch = '3C196'
#skymodel = '/home/fdg/scripts/model/3C295-allfield.skymodel'
#sourcedb = '/home/fdg/scripts/model/3C295-allfield.skydb'
#patch = '3C295'

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
os.system('mkdir globaldb-fr')
os.system('mkdir globaldb-clock')

nchan = find_nchan(mss[0])
logging.debug('Channel in the MS: '+str(nchan)+'.')

###############################################
# Initial processing (2/2013->2/2014)
#logging.warning('Fix beam table...')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)

# If more than 4 channels then average in freq to 4 chans
# TODO: avg to 5 sec?
if nchan > 4:
    if nchan % 4 != 0:
        logging.error('Channels should be a multiple of 4.')
        sys.exit(1)
    avg_factor = nchan / 4
    logging.info('Average in freq (factor of %i)...' % avg_factor)
    for ms in mss:
        msout = ms.replace('.MS','-avg.MS')
        s.add('NDPPP '+parset_dir+'/NDPPP-avg.parset msin='+ms+' msout='+msout+' msin.datacolumn=DATA avg.freqstep='+str(avg_factor), log=ms+'_avg.log', cmd_type='NDPPP')
    s.run(check=True)
    nchan = nchan / 4
    mss = sorted(glob.glob('*-avg.MS'))

############################################
# Prepare output parmdb
logging.info('Creating fake parmdb...')
for ms in mss:
    s.add('calibrate-stand-alone -f --parmdb-name instrument-clock '+ms+' '+parset_dir+'/bbs-fakeparmdb-clock.parset '+skymodel, log=ms+'_fakeparmdb-clock.log', cmd_type='BBS')
s.run(check=True)
for ms in mss:
    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
s.run(check=True)
for ms in mss:
    s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,24)"', log=ms+'_taql.log', cmd_type='general')
s.run(check=True)

# 1: find the FR and remve it

###############################################
# Beam correction DATA -> CORRECTED_DATA
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
s.run(check=True)

###############################################
# Convert to circular CORRECTED_DATA -> CORRECTED_DATA
logging.info('Converting to circular...')
for ms in mss:
    s.add('mslin2circ.py -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
s.run(check=True)

################################################
# Avg data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
# NOTE: the WEIGHTED_COLUMN is now smoothed in this dataset, a backup is in WEIGHTED_COLUMN_ORIG
logging.info('BL-averaging...')
for ms in mss:
    s.add('BLavg.py -r -w -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python')
s.run(check=True)

###############################################
# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logging.debug('Calibration - iterating on '+str(nchan)+' channels.')
for chan in xrange(nchan):
    logging.debug('Channel: '+str(chan))
    for ms in mss:
        check_rm(ms+'/instrument-'+str(chan))
        s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' msin.startchan='+str(chan)+' msin.nchan=1 cal.parmdb='+ms+'/instrument-'+str(chan)+' cal.sourcedb='+sourcedb+' cal.sources='+patch, log=ms+'-'+str(chan)+'_sol.log', cmd_type='NDPPP')
    s.run(check=True)

################################################
# Prepare and run losoto
for i, ms in enumerate(mss):
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb-fr/')
    num = re.findall(r'\d+', ms)[-1]
    for chan in xrange(nchan):
        logging.debug('Copy instrument-'+str(chan)+' of '+ms+' into globaldb/instrument-'+str(num)+'-'+str(chan))
        os.system('cp -r '+ms+'/instrument-'+str(chan)+' globaldb/instrument-'+str(num)+'-'+str(chan))
    logging.debug('Copy instrument-fr of '+ms+' into globaldb-fr/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument-fr globaldb/instrument-fr-'+str(num))

logging.info('Running LoSoTo...')
check_rm('plots')
os.makedirs('plots')
check_rm('cal.h5')
s.add('H5parm_importer.py -v cal.h5 globaldb', log='losoto.log', cmd_type='python', processors='max')
s.run(check=True)
s.add('losoto -v cal.h5 '+parset_dir+'/losoto-flag.parset', log='losoto-flag.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('losoto -v cal.h5 '+parset_dir+'/losoto-fr.parset', log='losoto-ph.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('H5parm_exporter.py -v cal.h5 globaldb-fr', log='losoto.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)

######################################################
# Correct FR/beam DATA -> CORRECTED_DATA
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-corFR.parset msin='+ms+' cal.parmdb='+ms+'/instrument-fr', log=ms+'_corFR.log', cmd_type='NDPPP')
s.run(check=True)

# 2: recalibrate without FR

###############################################
# Beam correction CORRECTED_DATA -> CORRECTED_DATA
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA', log=ms+'_beam2.log', cmd_type='NDPPP')
s.run(check=True)

################################################
# Avg data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
# NOTE: the WEIGHTED_COLUMN is now smoothed in this dataset, a backup is in WEIGHTED_COLUMN_ORIG
logging.info('BL-averaging...')
for ms in mss:
    s.add('BLavg.py -r -w -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python')
s.run(check=True)

###############################################
# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logging.debug('Calibration - iterating on '+str(nchan)+' channels.')
for chan in xrange(nchan):
    logging.debug('Channel: '+str(chan))
    for ms in mss:
        check_rm(ms+'/instrument-'+str(chan))
        s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' msin.startchan='+str(chan)+' msin.nchan=1 cal.parmdb='+ms+'/instrument-'+str(chan)+' cal.sourcedb='+sourcedb+' cal.sources='+patch, log=ms+'-'+str(chan)+'_cal.log', cmd_type='NDPPP')
    s.run(check=True)

#############################################################
# Prepare and run losoto
for i, ms in enumerate(mss):
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb-clock/')

    num = re.findall(r'\d+', ms)[-1]
    for chan in xrange(nchan):
        logging.debug('Copy instrument-'+str(chan)+' of '+ms+' into globaldb/instrument-'+str(num)+'-'+str(chan))
        os.system('cp -r '+ms+'/instrument-'+str(chan)+' globaldb/instrument-'+str(num)+'-'+str(chan))
    
    # We export clock, need to create a new parmdb
    logging.debug('Copy instrument-clock of '+ms+' into globaldb-clock/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument-clock globaldb-clock/instrument-'+str(num))

logging.info('Running LoSoTo...')
check_rm('plots')
os.makedirs('plots')
check_rm('cal.h5')
s.add('H5parm_importer.py -v cal.h5 globaldb', log='losoto.log', cmd_type='python', processors='max')
s.run(check=True)
s.add('losoto -v cal.h5 '+parset_dir+'/losoto-flag.parset', log='losoto-flag.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
os.system('cp -r cal.h5 cal.h5-flag')
s.add('losoto -v cal.h5 '+parset_dir+'/losoto-amp.parset', log='losoto-amp.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('losoto -v cal.h5 '+parset_dir+'/losoto-ph.parset', log='losoto-ph.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('H5parm_exporter.py -v cal.h5 globaldb-clock', log='losoto.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)

logging.info("Done.")
