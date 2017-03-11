#!/usr/bin/python
# data preparation for selfcal, apply cal solutions and split SB in time and concatenate in freq
# Input:
# Virgin target MSs and a globaldb of the calibrator
# Output:
# set of group*_TC*.MS file with DATA = calibrator corrected data, beam corrected, flagged

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from astropy.time import Time
from autocal.lib_pipeline import *

parset_dir = '/home/fdg/scripts/autocal/parset_timesplit'
initc = 0 # initial tc num (useful for multiple observation of same target) - tooth10==12

if 'tooth' in os.getcwd():
    # tooth
    ngroups = 2 # number of groups (totalSB/SBperFREQgroup)
    datadir = '.'
    globaldb = 'globaldb'
elif 'bootes' in os.getcwd(): # bootes 2013
    ngroups = 2
    datadir = '../tgts-bkp/' 
    globaldb = '../3c295/globaldb-clock/'
else:
    # survey
    ngroups = 1 # number of groups (totalSB/SBperFREQgroup)
    datadir = '/lofar5/stsf309/LBAsurvey/%s/%s' % (os.getcwd().split('/')[-2], os.getcwd().split('/')[-1]) # assumes e.g. ~/data/LBAsurvey/c05-o07/P155+52
    globaldb = '../3c196/globaldb' #TODO: copy form repository

##################################################################################################

set_logger('pipeline-timesplit.logging')
check_rm('logs')
s = Scheduler(dry=False)
assert os.path.isdir(globaldb)

#################################################
# Clear
logging.info('Cleaning...')

check_rm('*group*')
mss = sorted(glob.glob(datadir+'/*MS'))

##############################################
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

################################################
# Initial processing
# TODO: remove, moved to download pipeline
obs = pt.table(mss[0]+'/OBSERVATION', readonly=True, ack=False)
t = Time(obs.getcell('TIME_RANGE',0)[0]/(24*3600.), format='mjd')
time = np.int(t.iso.replace('-','')[0:8])
obs.close()
if time > 20130200 and time < 20140300:
    logging.warning('Fix beam table')
    for ms in mss:
        s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
    s.run(check=False)

# flag below elev 30
logging.info('Flagging elevation...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-flag-init.parset msin='+ms, log=ms+'_flag-init.log', cmd_type='NDPPP')
s.run(check=True)

####################################################
# Beam correction DATA -> CORRECTED_DATA (beam corrected+reweight)
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
s.run(check=True)

# To circular - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (beam corrected, circular)
#logging.info('Convert to circular...')
#for ms in mss:
#    s.add('/home/fdg/scripts/mslin2circ.py -s -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
#s.run(check=True)

# Copy instrument tables
for ms in mss:
    tnum = re.findall(r't\d+', ms)[0][1:]
    sbnum = re.findall(r'SB\d+', ms)[0][2:]
    check_rm(ms+'/instrument')
    logging.debug('cp -r '+globaldb+'/sol000_instrument-'+str(tnum)+'-'+str(sbnum)+' '+ms+'/instrument')
    os.system('cp -r '+globaldb+'/sol000_instrument-'+str(tnum)+'-'+str(sbnum)+' '+ms+'/instrument')

# Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (calibrator corrected+reweight, beam corrected, circular)
logging.info('Apply solutions...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor1.parmdb='+ms+'/instrument'+' cor2.parmdb='+ms+'/instrument', log=ms+'_cor.log', cmd_type='NDPPP')
s.run(check=True)

# Re-set weights of flagged data to 0 - this is necessary if we want to use dysco
# due to NDPPP leaving flagged weight at super-high values compared to unflagged ones
logging.info('Set weight of flagged data to 0...')
for ms in mss:
    s.add('flag_weight_to_zero.py '+ms, log=ms+'_resetweight.log', cmd_type='python')
s.run(check=True)

###################################################################################################
# Create groups
groupnames = []
logging.info('Concatenating in frequency...')
timechunks = set([re.findall(r't\d+', ms)[0][1:] for ms in mss ])
for timechunk in timechunks:
    for i, msg in enumerate(np.array_split(glob.glob('*_t'+timechunk+'_*MS'), ngroups)):
        if ngroups == 1:
            groupname = 'mss_t%s' % timechunk
        else:
            groupname = 'mss_t%s-%02i' % (timechunk, i)
        groupnames.append(groupname)
        check_rm(groupname)
        os.system('mkdir '+groupname)

        # add missing SB with a fake name not to leave frequency holes
        num_init = int(re.findall(r'\d+', msg[0])[-1])
        num_fin = int(re.findall(r'\d+', msg[-1])[-1])
        ms_name_init = msg[0]
        msg = []
        for j in range(num_init, num_fin+1):
            msg.append(ms_name_init.replace('SB%03i' % num_init, 'SB%03i' % j))

        # NOTE: dysco is added here! After fixing high weight of flagged data and remove of autocorr
        # prepare concatenated time chunks (TC) - SB.MS:CORRECTED_DATA -> group#.MS:DATA (cal corr data, beam corrected, circular)
        s.add('NDPPP '+parset_dir+'/NDPPP-concat.parset msin="['+','.join(msg)+']"  msout='+groupname+'/'+groupname+'.MS', \
                    log=groupname+'_NDPPP_concat.log', cmd_type='NDPPP')
    s.run(check=True)

# Flagging on concatenated dataset
logging.info('Flagging...')
for groupname in groupnames:
    s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+groupname+'/'+groupname+'.MS', \
                log=groupname+'_NDPPP_flag.log', cmd_type='NDPPP')
s.run(check=True)

# Create time-chunks
tc = initc
for groupname in groupnames:
    ms = groupname+'/'+groupname+'.MS'
    if not os.path.exists(ms): continue
    t = pt.table(ms, ack=False)
    starttime = t[0]['TIME']
    endtime   = t[t.nrows()-1]['TIME']
    hours = (endtime-starttime)/3600.
    logging.debug(ms+' has length of '+str(hours)+' h.')

    # split this ms into many TCs (#hours, i.e. chunks of 60 min)
    # to re-concat:
    #   t = table(['T0','T1',...])
    #   t.sort('TIME').copy('output.MS', deep = True)
    for timerange in np.array_split(t.getcol('TIME'), round(hours)):
        logging.debug('%02i - Splitting timerange %f %f' % (tc, timerange[0], timerange[-1]))
        t1 = t.query('TIME >= ' + str(timerange[0]) + ' && TIME <= ' + str(timerange[-1]), sortlist='TIME,ANTENNA1,ANTENNA2')
        splitms = groupname+'/TC%02i.MS' % tc
        check_rm(splitms)
        t1.copy(splitms, True)
        t1.close()
        tc += 1
    t.close()

    check_rm(ms) # remove not-timesplitted file

# put everything together
for group in xrange(ngroups):
    if ngroups == 1:
        groupname = 'mss'
    else:
        groupname = 'mss-%02i' % group
 
    check_rm(groupname)
    os.makedirs(groupname)
    os.system('cp -r mss_t*_%02i/*MS %s' % (group, groupname))

logging.info("Done.")
