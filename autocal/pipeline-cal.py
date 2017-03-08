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
elif calname == '3c295':
    sourcedb = '/home/fdg/scripts/model/calib-simple.skydb'
    patch = '3C295'
elif calname == 'CygA':
    sourcedb = '/home/fdg/scripts/model/A-team_4_CC.skydb'
    patch = 'CygA'

def run_losoto(c, mss, parsets, outtab='', inglobaldb='globaldb', outglobaldb='globaldb', ininstrument='instrument', outinstrument='instrument', putback=False):
    """
    c : cycle name, e.g. "final"
    mss : lists of MS files
    parsets : lists of parsets to execute
    outtab : strings with soltab to output e.g. 'amplitudeSmooth000,phaseOrig000'
    putback : put back in MS the instrument tables
    """

    logging.info('Running LoSoTo...')

    # prepare globaldbs
    check_rm('plots-'+c)
    check_rm(inglobaldb)
    os.system('mkdir '+inglobaldb)
    if inglobaldb != outglobaldb: 
        check_rm(outglobaldb)
        os.system('mkdir '+outglobaldb)

    for i, ms in enumerate(mss):
        if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky '+inglobaldb)
        if inglobaldb != outglobaldb:
            if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky '+outglobaldb)
    
        tnum = re.findall(r'\d+', ms)[-2]
        sbnum = re.findall(r'\d+', ms)[-1]
        os.system('cp -r '+ms+'/'+ininstrument+' '+inglobaldb+'/instrument-'+str(tnum)+'-'+str(sbnum))
       
        if inglobaldb != outglobaldb:
            os.system('cp -r '+ms+'/'+outinstrument+' '+outglobaldb+'/instrument-'+str(tnum)+'-'+str(sbnum))
    
    check_rm('plots')
    os.makedirs('plots')
    check_rm('cal-'+c+'.h5')
    
    s.add('H5parm_importer.py -v cal-'+c+'.h5 globaldb', log='losoto-'+c+'.log', cmd_type='python')
    s.run(check=True)
    
    for parset in parsets:
        logging.debug('-- executing '+parset+'...')
        s.add('losoto -v cal-'+c+'.h5 '+parset_dir+'/'+parset, log='losoto-'+c+'.log', log_append=True, cmd_type='python', processors='max')
        s.run(check=True)

    os.system('mv plots plots-'+c)
    
    if outtab != '':
        s.add('H5parm_exporter.py -v -c --soltab '+outtab+' cal-'+c+'.h5 '+outglobaldb, log='losoto-'+c+'.log', log_append=True, cmd_type='python')
        s.run(check=True)

    if putback:
        for i, ms in enumerate(mss):
            tnum = re.findall(r'\d+', ms)[-2]
            num = re.findall(r'\d+', ms)[-1]
            check_rm(ms+'/'+outinstrument)
            os.system('cp -r '+outglobaldb+'/sol000_instrument-'+str(tnum)+'-'+str(num)+' '+ms+'/'+outinstrument)


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
#     if os.path.exists(ms+'/instrument-clock'): continue
##    s.add('calibrate-stand-alone -f --parmdb-name instrument-clock '+ms+' '+parset_dir+'/bbs-fakeparmdb-clock.parset '+skymodel, log=ms+'_fakeparmdb-clock.log', cmd_type='BBS')
##s.run(check=True)
#for ms in mss:
#    if os.path.exists(ms+'/instrument-fr'): continue
#    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
#s.run(check=True)
#for ms in mss:
#    s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,24)"', log=ms+'_taql.log', cmd_type='general')
#s.run(check=True)

#################################################
# 1: find the FR and remve it

# Beam correction DATA -> CORRECTED_DATA
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
s.run(check=True)

# TESTTESTTEST
# Correct DELAY CORRECTED_DATA (beam corrected) -> CORRECTED_DATA
#logging.info('Cross delay correction...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-cd cor.correction=gain', log=ms+'_corCD.log', cmd_type='NDPPP')
#s.run(check=True)

# Convert to circular CORRECTED_DATA -> CORRECTED_DATA
logging.info('Converting to circular...')
for ms in mss:
    s.add('mslin2circ.py -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
s.run(check=True)

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logging.info('BL-smooth...')
for ms in mss:
    s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth1.log', cmd_type='python')
s.run(check=True)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logging.info('Calibrating...')
for ms in mss:
    check_rm(ms+'/instrument')
    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' sol.sourcedb='+sourcedb+' sol.sources='+patch, log=ms+'_sol1.log', cmd_type='NDPPP')
s.run(check=True)

run_losoto('fr', mss, ['losoto-fr.parset'], outtab='rotationmeasure000', \
    inglobaldb='globaldb', outglobaldb='globaldb-fr', ininstrument='instrument', outinstrument='instrument-fr', putback=True)

#####################################################
# 2: find CROSS DELAY and remve it

# Beam correction DATA -> CORRECTED_DATA
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
s.run(check=True)

# Correct FR CORRECTED_DATA -> CORRECTED_DATA
logging.info('Faraday rotation correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', log=ms+'_corFR.log', cmd_type='NDPPP')
s.run(check=True)

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logging.info('BL-smooth...')
for ms in mss:
    s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth2.log', cmd_type='python')
s.run(check=True)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logging.info('Calibrating...')
for ms in mss:
    check_rm(ms+'/instrument')
    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' sol.sourcedb='+sourcedb+' sol.sources='+patch, log=ms+'_sol2.log', cmd_type='NDPPP')
s.run(check=True)

run_losoto('cd', mss, ['losoto-flag.parset','losoto-cd.parset'], outtab='amplitude000,crossdelay', \
    inglobaldb='globaldb', outglobaldb='globaldb', ininstrument='instrument', outinstrument='instrument-cd', putback=True)

##################################################
# 3: recalibrate without FR

# Beam correction DATA -> CORRECTED_DATA
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam2.log', cmd_type='NDPPP')
s.run(check=True)

# Correct DELAY CORRECTED_DATA (beam corrected) -> CORRECTED_DATA
logging.info('Cross delay correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-cd cor.correction=gain', log=ms+'_corCD.log', cmd_type='NDPPP')
s.run(check=True)

# Correct FR CORRECTED_DATA -> CORRECTED_DATA
logging.info('Faraday rotation correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', log=ms+'_corFR.log', cmd_type='NDPPP')
s.run(check=True)

# Convert to circular CORRECTED_DATA -> CORRECTED_DATA
#logging.info('Converting to circular...')
#for ms in mss:
#    s.add('mslin2circ.py -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
#s.run(check=True)

# Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
logging.info('BL-smooth...')
for ms in mss:
    s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth3.log', cmd_type='python')
s.run(check=True)

# Solve cal_SB.MS:SMOOTHED_DATA (only solve)
logging.info('Calibrating...')
for ms in mss:
    check_rm(ms+'/instrument')
    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' sol.sourcedb='+sourcedb+' sol.sources='+patch, log=ms+'_sol3.log', cmd_type='NDPPP')
s.run(check=True)

run_losoto('final', mss, ['losoto-flag.parset','losoto-amp.parset','losoto-ph.parset'], outtab='amplitudeSmooth000,phaseOrig000', \
    inglobaldb='globaldb', outglobaldb='globaldb', ininstrument='instrument', outinstrument='instrument', putback=False)
#run_losoto('final', mss, ['losoto-flag.parset','losoto-amp.parset','losoto-ph.parset'], outtab='amplitudeSmooth000,phase000,clock000', \
#    inglobaldb='globaldb', outglobaldb='globaldb-clock', ininstrument='instrument-clock', outinstrument='instrument', putback=False)

# TESTTESTTEST
#run_losoto('final', mss, ['losoto-amp.parset'], outtab='amplitudeSmooth000,phaseOrig000', \
#    inglobaldb='globaldb', outglobaldb='globaldb', ininstrument='instrument', outinstrument='instrument', putback=False)

logging.info("Done.")
