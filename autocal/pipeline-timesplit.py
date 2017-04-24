#!/usr/bin/python
# data preparation for selfcal, apply cal solutions and split SB in time and concatenate in freq
# Input:
# Virgin target MSs and a globaldb of the calibrator
# Output:
# set of group*_TC*.MS file with DATA = calibrator corrected data, beam corrected, flagged

import sys, os, glob, re
import numpy as np
from autocal.lib_pipeline import *
from astropy.time import Time
import casacore.tables as pt

parset_dir = '/home/fdg/scripts/autocal/parset_timesplit'
initc = 0 # initial tc num (useful for multiple observation of same target) - tooth10==12
clock = True

if 'LBAsurvey' in os.getcwd():
    ngroups = 1 # number of groups (totalSB/SBperFREQgroup)
    datadir = '/lofar5/stsf309/LBAsurvey/%s/%s' % (os.getcwd().split('/')[-2], os.getcwd().split('/')[-1]) # assumes e.g. ~/data/LBAsurvey/c05-o07/P155+52
    globaldb = 'globaldb_'+os.getcwd().split('/')[-2]
    logging.info('Copy: dsk:/disks/paradata/fdg/LBAsurvey/%s -> .' % globaldb)
    os.system('scp dsk:/disks/paradata/fdg/LBAsurvey/%s .' % globaldb) # TODO: move only _tXXX files
else:
    ngroups = 2
    if clock:
        globaldb = '../cals/globaldb-clock'
    else:
        globaldb = '../cals/globaldb'
    datadir = '../tgts-bkp/' 

##################################################################################################
set_logger('pipeline-timesplit.logging')
check_rm('logs')
s = Scheduler(dry=False)
assert os.path.isdir(globaldb)

#################################################
# Clear
logging.info('Cleaning...')

check_rm('mss*')
mss = sorted(glob.glob(datadir+'/*MS'))

##############################################
# Avg to 4 chan and 4 sec
# Remove internationals
# TODO: move to download pipeline
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
    os.system('cp -r '+globaldb+'/sol000_instrument-'+str(tnum)+'-'+str(sbnum)+' '+ms+'/instrument')

# Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (calibrator corrected+reweight, beam corrected, circular)
logging.info('Apply solutions...')
for ms in mss:
    if clock:
        s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' steps=[cor1,cor2] cor1.parmdb='+ms+'/instrument'+' cor2.parmdb='+ms+'/instrument', log=ms+'_cor.log', cmd_type='NDPPP')
    else:
        s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' steps=[cor1] cor1.parmdb='+ms+'/instrument', log=ms+'_cor.log', cmd_type='NDPPP')
s.run(check=True)

# Re-set weights of flagged data to 0 - this is necessary if we want to use dysco
# due to NDPPP leaving flagged weight at super-high values compared to unflagged ones
#logging.info('Set weight of flagged data to 0...')
#for ms in mss:
#    s.add('flag_weight_to_zero.py '+ms, log=ms+'_resetweight.log', cmd_type='python')
#s.run(check=True)

###################################################################################################
# Create groups
groupnames = []
logging.info('Concatenating in frequency...')
timechunks = set([re.findall(r'_t\d+', ms)[0][2:] for ms in mss ])
for timechunk in timechunks:
    for i, msg in enumerate(np.array_split(sorted(glob.glob('*_t'+timechunk+'_*MS')), ngroups)):
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
if ngroups == 1:
    check_rm('mss')
    os.makedirs('mss')
    os.system('mv mss_t*/*MS mss')
else:
    for group in xrange(ngroups):
        groupname = 'mss-%02i' % group
        check_rm(groupname)
        os.makedirs(groupname)
        os.system('mv mss_t*-%02i/*MS %s' % (group, groupname))
    check_rm('mss_t*')

logging.info("Done.")
