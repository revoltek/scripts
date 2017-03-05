#!/usr/bin/python
# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation
import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from astropy.time import Time
from autocal.lib_pipeline import *

parset_dir = '/home/fdg/scripts/autocal/parset_cal'
skymodel = '/home/fdg/scripts/model/calib-simple.skymodel'

if 'tooth' in os.getcwd():
    calname = '3c196'
    datadir = '../cals-bkp/'
    bl2flag = 'CS031LBA\;RS409LBA'
elif 'bootes' in os.getcwd(): # bootes 2013
    calname = os.getcwd().split('/')[-1]
    datadir = '../cals-bkp/'
    bl2flag = 'CS013LBA\;CS031LBA'
elif 'daycomm' in os.getcwd(): # daytest
    calname = os.getcwd().split('/')[-1]
    datadir = '/data/scratch/COMMISSIONING2017/c07-o01/%s' % calname
    bl2flag = 'CS031LBA'
else:
    obs = os.getcwd().split('/')[-2] # assumes .../c05-o07/3c196
    calname = os.getcwd().split('/')[-1] # assumes .../c05-o07/3c196
    datadir = '/lofar5/stsf309/LBAsurvey/%s/%s' % (obs, calname)
    #bl2flag = 'CS031LBA\;RS409LBA\;RS310LBA\;RS210LBA\;RS407LBA'
    bl2flag = 'CS031LBA'


if calname == '3c196':
    sourcedb = '/home/fdg/scripts/model/3C196-allfield.skydb'
    patch = '3C196'
elif calname == '3c380':
    sourcedb = '/home/fdg/scripts/model/calib-simple.skydb'
    patch = '3C380'
elif calname == 'CygA':
    sourcedb = '/home/fdg/scripts/model/A-team_4_CC.skydb'
    patch = 'CygA'

###################################################

set_logger('pipeline-cal.logging')
check_rm('logs')
s = Scheduler(dry=False)
mss = sorted(glob.glob(datadir+'/*MS'))

############################################################
# Avg to 4 chan and 4 sec
# Remove internationals

nchan = find_nchan(mss[0])
timeint = find_timeint(mss[0])
if nchan % 4 != 0 and nchan != 1:
    logging.error('Channels should be a multiple of 4.')
    sys.exit(1)

avg_factor_f = nchan / 4 # to 4 ch/SB
if avg_factor_f < 1: avg_factor_f = 1
avg_factor_t = int(np.round(4/timeint)) # to 4 sec
if avg_factor_t < 1: avg_factor_t = 1

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
else:
    logging.info('Copy data - no averaging...')
    for ms in mss:
        msout = ms.replace('.MS','-avg.MS').split('/')[-1]
        if os.path.exists(msout): continue
        os.system('cp -r '+ms+' '+msout)

mss = sorted(glob.glob('*.MS'))

###########################################################
## flag below elev 20 and bad stations, flags will propagate
#logging.info('Flagging...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms+' msout=. flag1.baseline='+bl2flag+' msin.datacolumn=DATA', \
#            log=ms+'_flag.log', cmd_type='NDPPP')
#s.run(check=True)
#    
## Initial processing (2/2013->2/2014)
#obs = pt.table(mss[0]+'/OBSERVATION', readonly=True, ack=False)
#t = Time(obs.getcell('TIME_RANGE',0)[0]/(24*3600.), format='mjd')
#time = np.int(t.iso.replace('-','')[0:8])
#obs.close()
#if time > 20130200 and time < 20140300:
#    logging.info('Fix beam table...')
#    for ms in mss:
#        s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#    s.run(check=False)
#
## Prepare output parmdb
## TODO: remove as soon as losoto has the proper exporter
#logging.info('Creating fake parmdb...')
##for ms in mss:
##    s.add('calibrate-stand-alone -f --parmdb-name instrument-clock '+ms+' '+parset_dir+'/bbs-fakeparmdb-clock.parset '+skymodel, log=ms+'_fakeparmdb-clock.log', cmd_type='BBS')
##s.run(check=True)
#for ms in mss:
#    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
#s.run(check=True)
#for ms in mss:
#    s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,24)"', log=ms+'_taql.log', cmd_type='general')
#s.run(check=True)
#
## 1: find CROSS DELAY and remve it
#
## Beam correction DATA -> CORRECTED_DATA
#logging.info('Beam correction...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
#s.run(check=True)
#
## Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
#logging.info('BL-smooth...')
#for ms in mss:
#    s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python')
#s.run(check=True)
#
## Solve cal_SB.MS:SMOOTHED_DATA (only solve)
#logging.info('Calibrating...')
#for ms in mss:
#    check_rm(ms+'/instrument')
#    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' sol.sourcedb='+sourcedb+' sol.sources='+patch, log=ms+'_sol1.log', cmd_type='NDPPP')
#s.run(check=True)
#
## Prepare and run losoto
#logging.info('Running LoSoTo...')
#check_rm('globaldb')
#os.system('mkdir globaldb')
#check_rm('plots-cd')
#for i, ms in enumerate(mss):
#    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
#    tnum = re.findall(r'\d+', ms)[-2]
#    sbnum = re.findall(r'\d+', ms)[-1]
#    #logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(tnum)+'-'+str(sbnum))
#    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(tnum)+'-'+str(sbnum))
#
#check_rm('plots')
#check_rm('cal1.h5')
#s.add('H5parm_importer.py -v cal1.h5 globaldb', log='losoto1.log', cmd_type='python', processors=1)
#s.run(check=True)
#s.add('losoto -v cal1.h5 '+parset_dir+'/losoto-flag.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
#os.system('cp -r cal1.h5 cal1.h5-flag')
#s.add('losoto -v cal1.h5 '+parset_dir+'/losoto-cd.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
#
#s.add('H5parm_exporter.py -v -t amplitude000,crossdelay cal1.h5 globaldb', log='losoto1.log', log_append=True, cmd_type='python', processors=1)
#s.run(check=True)
#os.system('mv plots plots-cd')
#
#for i, ms in enumerate(mss):
#    tnum = re.findall(r'\d+', ms)[-2]
#    num = re.findall(r'\d+', ms)[-1]
#    check_rm(ms+'/instrument-cd')
#    #logging.debug('Copy globaldb/sol000_instrument-'+str(tnum)+'-'+str(num)+' into '+ms+'/instrument-cd')
#    os.system('cp -r globaldb/sol000_instrument-'+str(tnum)+'-'+str(num)+' '+ms+'/instrument-cd')
#
##################################################
## 2: find the FR and remve it
#
## Correct DELAY CORRECTED_DATA (beam corrected) -> CORRECTED_DATA
## TODO: before the beam? definitely before going to circular!
#logging.info('Cross delay correction...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-cd cor.correction=gain', log=ms+'_corCD.log', cmd_type='NDPPP')
#s.run(check=True)
#
## Convert to circular CORRECTED_DATA -> CORRECTED_DATA
#logging.info('Converting to circular...')
#for ms in mss:
#    s.add('mslin2circ.py -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
#s.run(check=True)
#
## Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
#logging.info('BL-smooth...')
#for ms in mss:
#    s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python')
#s.run(check=True)
#
## Solve cal_SB.MS:SMOOTHED_DATA (only solve)
#logging.info('Calibrating...')
#for ms in mss:
#    check_rm(ms+'/instrument')
#    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' sol.sourcedb='+sourcedb+' sol.sources='+patch, log=ms+'_sol2.log', cmd_type='NDPPP')
#s.run(check=True)
#
## Prepare and run losoto
#logging.info('Running LoSoTo...')
#check_rm('globaldb')
#check_rm('globaldb-fr')
#os.system('mkdir globaldb')
#os.system('mkdir globaldb-fr')
#check_rm('plots-fr')
#for i, ms in enumerate(mss):
#    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
#    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb-fr/')
#    tnum = re.findall(r'\d+', ms)[-2]
#    sbnum = re.findall(r'\d+', ms)[-1]
#    #logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(tnum)+'-'+str(sbnum))
#    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(tnum)+'-'+str(sbnum))
#    #logging.debug('Copy instrument-fr of '+ms+' into globaldb-fr/instrument-'+str(tnum)+'-'+str(sbnum))
#    os.system('cp -r '+ms+'/instrument-fr globaldb-fr/instrument-fr-'+str(tnum)+'-'+str(sbnum))
#
#check_rm('plots')
#check_rm('cal2.h5')
#s.add('H5parm_importer.py -v cal2.h5 globaldb', log='losoto2.log', cmd_type='python', processors=1)
#s.run(check=True)
#s.add('losoto -v cal2.h5 '+parset_dir+'/losoto-flag.parset', log='losoto2.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
#os.system('cp -r cal2.h5 cal2.h5-flag')
#s.add('losoto -v cal2.h5 '+parset_dir+'/losoto-fr.parset', log='losoto2.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
## TESTTESTTEST
##s.add('losoto -v cal2.h5 '+parset_dir+'/losoto-amp.parset', log='losoto2.log', log_append=True, cmd_type='python', processors='max')
##s.run(check=True)
##s.add('losoto -v cal2.h5 '+parset_dir+'/losoto-ph.parset', log='losoto2.log', log_append=True, cmd_type='python', processors='max')
##s.run(check=True)
#
#s.add('H5parm_exporter.py -v -t rotationmeasure000 cal2.h5 globaldb-fr', log='losoto2.log', log_append=True, cmd_type='python', processors=1)
#s.run(check=True)
#os.system('mv plots plots-fr')
#
#for i, ms in enumerate(mss):
#    tnum = re.findall(r'\d+', ms)[-2]
#    num = re.findall(r'\d+', ms)[-1]
#    check_rm(ms+'/instrument-fr')
#    #logging.debug('Copy globaldb-fr/sol000_instrument-fr-'+str(tnum)+'-'+str(num)+' into '+ms+'/instrument-fr')
#    os.system('cp -r globaldb-fr/sol000_instrument-fr-'+str(tnum)+'-'+str(num)+' '+ms+'/instrument-fr')

##################################################
# 3: recalibrate without FR

# Beam correction DATA -> CORRECTED_DATA
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam2.log', cmd_type='NDPPP')
s.run(check=True)

# Convert to linear CORRECTED_DATA -> CORRECTED_DATA
#logging.info('Converting to linear...')
#for ms in mss:
#    s.add('mslin2circ.py -w -r -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
#s.run(check=True)

# Correct FR CORRECTED_DATA -> CORRECTED_DATA
logging.info('Faraday rotation correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', log=ms+'_corFR.log', cmd_type='NDPPP')
s.run(check=True)
# Correct DELAY CORRECTED_DATA (beam corrected) -> CORRECTED_DATA
# TODO: before the beam? definitely before going to circular!
logging.info('Cross delay correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-cd cor.correction=gain', log=ms+'_corCD.log', cmd_type='NDPPP')
s.run(check=True)

# Convert to circular CORRECTED_DATA -> CORRECTED_DATA
# NOTE: in linear
#logging.info('Converting to circular...')
#for ms in mss:
#    s.add('mslin2circ.py -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
#s.run(check=True)

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logging.info('BL-smoothing...')
for ms in mss:
    s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python')
s.run(check=True)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logging.info('Calibrating...')
for ms in mss:
    check_rm(ms+'/instrument')
    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' sol.sourcedb='+sourcedb+' sol.sources='+patch, log=ms+'_sol3.log', cmd_type='NDPPP')
s.run(check=True)

# Prepare and run losoto
logging.info('Running LoSoTo...')
check_rm('globaldb') # remove it as it was used for the fr
#check_rm('globaldb-clock')
os.system('mkdir globaldb')
#os.system('mkdir globaldb-clock')
check_rm('plots-final')
for i, ms in enumerate(mss):
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
#    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb-clock/')

    tnum = re.findall(r'\d+', ms)[-2]
    sbnum = re.findall(r'\d+', ms)[-1]
    #logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(tnum)+'-'+str(sbnum))
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(tnum)+'-'+str(sbnum))
   
#    # We export clock, need to create a new parmdb
#    logging.debug('Copy instrument-clock of '+ms+' into globaldb-clock/instrument-'+str(num))
#    os.system('cp -r '+ms+'/instrument-clock globaldb-clock/instrument-'+str(num))

check_rm('plots')
os.makedirs('plots')
check_rm('cal3.h5')

s.add('H5parm_importer.py -v cal3.h5 globaldb', log='losoto3.log', cmd_type='python', processors=1)
s.run(check=True)

#s.add('losoto -v cal3.h5 '+parset_dir+'/losoto-flag.parset', log='losoto3.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
#os.system('cp -r cal3.h5 cal3.h5-flag')

s.add('losoto -v cal3.h5 '+parset_dir+'/losoto-amp.parset', log='losoto3.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)

s.add('losoto -v cal3.h5 '+parset_dir+'/losoto-ph.parset', log='losoto3.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)

# copy clock+BP
#s.add('H5parm_exporter.py -v -c --soltab amplitudeSmooth000,phase000,clock000 cal3.h5 globaldb-clock', log='losoto3.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)

# copy ph+BP
s.add('H5parm_exporter.py -v -c --soltab amplitudeSmooth000,phaseOrig000 cal3.h5 globaldb', log='losoto3.log', log_append=True, cmd_type='python', processors=1)
s.run(check=True)

os.system('mv plots plots-final')

logging.info("Done.")
