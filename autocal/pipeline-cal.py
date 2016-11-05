#!/usr/bin/python
# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation
import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from lib_pipeline import *

skymodel = '/home/fdg/scripts/model/3C196-allfield.skymodel'
sourcedb = '/home/fdg/scripts/model/3C196-allfield.skydb'
patch = '3C196'
#skymodel = '/home/fdg/scripts/model/3C295-allfield.skymodel'
#sourcedb = '/home/fdg/scripts/model/3C295-allfield.skydb'
#patch = '3C295'

parset_dir = '/home/fdg/scripts/autocal/parset_cal'
datadir = '/lofar5/stsf309/LBAsurvey/%s/3c196' % os.getcwd().split('/')[-2] # assumes ~/data/LBAsurvey/c05-o07/3c196
#datadir = '.'

###################################################

set_logger()
check_rm('logs')
s = Scheduler(dry=False)
mss = sorted(glob.glob(datadir+'/*MS'))

###########################################################
# Avg to 4 chan and 4 sec
# Remove internationals
nchan = find_nchan(mss[0])
timeint = find_timeint(mss[0])
if nchan % 4 != 0:
    logging.error('Channels should be a multiple of 4.')
    sys.exit(1)

avg_factor_f = nchan / 4 # to 4 ch/SB
if avg_factor_f < 1: avg_factor_f = 1
avg_factor_t = int(np.round(4/timeint))
if avg_factor_t < 1: avg_factor_t = 1 # to 4 sec

if avg_factor_f != 1 or avg_factor_t != 1:
    logging.info('Average in freq (factor of %i) and time (factor of %i)...' % (avg_factor_f, avg_factor_t))
    for ms in mss:
        msout = ms.replace('.MS','-avg.MS').split('/')[-1]
        if os.path.exists(msout): continue
        s.add('NDPPP '+parset_dir+'/NDPPP-avg.parset msin='+ms+' msout='+msout+' msin.datacolumn=DATA avg.timestep='+str(avg_factor_t)+' avg.freqstep='+str(avg_factor_f), \
                log=msout+'_avg.log', cmd_type='NDPPP')
    s.run(check=True)
    nchan = nchan / avg_factor_f
    timeint = timeint * avg_factor_t
    mss = sorted(glob.glob('*-avg.MS'))
    
###############################################
# Initial processing (2/2013->2/2014)
#logging.warning('Fix beam table...')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)

############################################
# Prepare output parmdb
# TODO: remove as soon as losoto has the proper exporter
logging.info('Creating fake parmdb...')
#for ms in mss:
#    s.add('calibrate-stand-alone -f --parmdb-name instrument-clock '+ms+' '+parset_dir+'/bbs-fakeparmdb-clock.parset '+skymodel, log=ms+'_fakeparmdb-clock.log', cmd_type='BBS')
#s.run(check=True)
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
    s.add('mslin2circ.py -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
s.run(check=True)

#################################################
# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
# NOTE: the WEIGHTED_COLUMN is now smoothed in this dataset, a backup is in WEIGHTED_COLUMN_ORIG
logging.info('BL-smooth...')
for ms in mss:
    s.add('BLavg.py -r -w -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python')
s.run(check=True)

################################################
# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logging.info('Calibrating...')
for ms in mss:
    check_rm(ms+'/instrument')
    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' cal.parmdb='+ms+'/instrument cal.sourcedb='+sourcedb+' cal.sources='+patch, log=ms+'_sol-circ.log', cmd_type='NDPPP')
s.run(check=True)

################################################
# Prepare and run losoto
check_rm('globaldb')
check_rm('globaldb-fr')
os.system('mkdir globaldb')
os.system('mkdir globaldb-fr')
for i, ms in enumerate(mss):
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb-fr/')
    num = re.findall(r'\d+', ms)[-1]
    logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
    logging.debug('Copy instrument-fr of '+ms+' into globaldb-fr/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument-fr globaldb-fr/instrument-fr-'+str(num))

logging.info('Running LoSoTo...')
check_rm('plots')
check_rm('cal1.h5')
s.add('H5parm_importer.py -v cal1.h5 globaldb', log='losoto1.log', cmd_type='python', processors='max')
s.run(check=True)
s.add('losoto -v cal1.h5 '+parset_dir+'/losoto-flag.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
os.system('cp -r cal1.h5 cal1.h5-flag')
s.add('losoto -v cal1.h5 '+parset_dir+'/losoto-fr.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('H5parm_exporter.py -v -t rotationmeasure000 cal1.h5 globaldb-fr', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
check_rm('plots-fr')
os.system('mv plots plots-fr')

for i, ms in enumerate(mss):
    num = re.findall(r'\d+', ms)[-1]
    check_rm(ms+'/instrument-fr')
    logging.debug('Copy globaldb-fr/sol000_instrument-fr-'+str(num)+' into '+ms+'/instrument-fr')
    os.system('cp -r globaldb-fr/sol000_instrument-fr-'+str(num)+' '+ms+'/instrument-fr')

# 2: recalibrate without FR

###############################################
# Beam correction DATA -> CORRECTED_DATA
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam2.log', cmd_type='NDPPP')
s.run(check=True)

######################################################
# Correct FR CORRECTED_DATA -> CORRECTED_DATA
logging.info('Faraday rotation correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-corFR.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr', log=ms+'_corFR.log', cmd_type='NDPPP')
s.run(check=True)

###############################################
# Convert to circular CORRECTED_DATA -> CORRECTED_DATA
# TESTTESTTEST
#logging.warning('Converting to circular...')
#for ms in mss:
#    s.add('mslin2circ.py -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
#s.run(check=True)

################################################
# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logging.info('BL-smoothing...')
for ms in mss:
    s.add('BLavg.py -r -w -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python')
s.run(check=True)

################################################
# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logging.info('Calibrating...')
for ms in mss:
    check_rm(ms+'/instrument')
    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' cal.parmdb='+ms+'/instrument cal.sourcedb='+sourcedb+' cal.sources='+patch, log=ms+'_sol-lin.log', cmd_type='NDPPP')
s.run(check=True)

#############################################################
# Prepare and run losoto
check_rm('globaldb') # remove it as it was used for the fr
#check_rm('globaldb-clock')
os.system('mkdir globaldb')
#os.system('mkdir globaldb-clock')
for i, ms in enumerate(mss):
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
#    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb-clock/')

    num = re.findall(r'\d+', ms)[-1]
    logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
   
#    # We export clock, need to create a new parmdb
#    logging.debug('Copy instrument-clock of '+ms+' into globaldb-clock/instrument-'+str(num))
#    os.system('cp -r '+ms+'/instrument-clock globaldb-clock/instrument-'+str(num))

logging.info('Running LoSoTo...')
check_rm('plots')
os.makedirs('plots')
check_rm('cal2.h5')

s.add('H5parm_importer.py -v cal2.h5 globaldb', log='losoto2.log', cmd_type='python', processors='max')
s.run(check=True)

# TESTTESTTEST
#os.system('cp -r cal2.h5 cal2.h5-bkp')
#s.add('losoto -v cal2.h5 '+parset_dir+'/losoto-fr.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)

s.add('losoto -v cal2.h5 '+parset_dir+'/losoto-flag.parset', log='losoto2.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
os.system('cp -r cal2.h5 cal2.h5-flag')

s.add('losoto -v cal2.h5 '+parset_dir+'/losoto-amp.parset', log='losoto2.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)

s.add('losoto -v cal2.h5 '+parset_dir+'/losoto-ph.parset', log='losoto2.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)

# copy clock+BP
#s.add('H5parm_exporter.py -v -c --soltab amplitudeSmooth000,phase000,clock000 cal2.h5 globaldb-clock', log='losoto2.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)

# copy ph+BP
s.add('H5parm_exporter.py -v -c --soltab amplitudeSmooth000,phaseOrig000 cal2.h5 globaldb', log='losoto2.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)

logging.info("Done.")
