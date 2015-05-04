#!/usr/bin/python
# data preparation for selfcal, apply cal solutions and split SB in time and concatenate in freq
# Input:
# Virgin target MSs and a globaldb of the calibrator
# Output:
# set of group*_TC*.MS file with DATA = calibrator corrected data, beam corrupted

ngroups = 12 # number of groups (=total_time/chunk_time)
initc = 0 # initial tc num (useful for multiple observation of same target) - tooth10==12
fakeskymodel = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/toothbrush.fakemodel.skymodel'
globaldb = '../cals/globaldb'

##################################################################################################

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from lib_pipeline import *

set_logger()
s = Scheduler(qsub=True, max_threads=30, dry=False)

#################################################
# Clear
logging.info('Cleaning...')
check_rm('*log')
check_rm('*group*')

##############################################
# Initial processing
logging.info('Fix beam table')
for ms in sorted(glob.glob('*MS')):
    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
s.run(check=False)

#################################################
# Copy cal solution
for ms in sorted(glob.glob('*MS')):
    num = re.findall(r'\d+', ms)[-1]
    check_rm(ms+'/instrument')
    logging.debug('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')
    os.system('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')

###################################################################################################
# [PARALLEL] Apply cal sol - SB.MS:DATA -> SB.MS:CALCOR_DATA (calibrator corrected data, beam corrected, linear)
logging.info('Apply solutions...')
for ms in glob.glob('*MS'):
    s.add('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_timesplit/bbs_correct.parset '+fakeskymodel, log=ms+'_cor.log', cmd_type='BBS')
s.run(check=True)

###################################################################################################
# [PARALLEL] To circular - SB.MS:CALCOR_DATA -> SB.MS:CALCOR_DATA_CIRC (calibrator corrected data, beam corrected, circular)
logging.info('Convert to circular...')
for ms in glob.glob('*MS'):
    s.add('/home/fdg/scripts/mslin2circ.py -i '+ms+':CALCOR_DATA -o '+ms+':CALCOR_DATA_CIRC', log=ms+'_circ2lin.log', cmd_type='python')
s.run(check=True)

# TODO: combine all SB and run aoflagger

###################################################################################################
# split each MS in timechunks of 1 h and create groups TODO: add possibility to create uneaven groups (e.g. more SB at low-freq)
logging.info('Splitting MSs...')
for i, msg in enumerate(np.array_split(sorted(glob.glob('*.MS')), ngroups)):
    logging.info('Working on group: '+str(i))
    check_rm('group'+str(i))
    os.system('mkdir group'+str(i))

    # iterate ms
    for ms in msg:
        t = pt.table(ms)
        starttime = t[0]['TIME']
        endtime   = t[t.nrows()-1]['TIME']
        hours = (endtime-starttime)/3600.
        logging.debug(ms+' has length of '+str(hours)+' h.')

        # split this ms into many TCs (2*#hours)
        for tc, timerange in enumerate(np.array_split(t.getcol('TIME'), 2*round(hours))):
            tc += initc
            logging.debug('Splitting timerange '+str(timerange[0])+' - '+str(timerange[-1]))
            t1 = t.query('TIME >= ' + str(timerange[0]) + ' && \
                      TIME <= ' + str(timerange[-1]), sortlist='TIME,ANTENNA1,ANTENNA2')
            splitms = ms.replace('.MS','_group'+str(i)+'_TC'+str(tc)+'.MS')
            check_rm(splitms)
            t1.copy(splitms, True)
            t1.close()
            os.system('mv '+splitms+' group'+str(i))

        t.close()

    tcnums = []
    for ms in sorted(glob.glob('group'+str(i)+'/*.MS')):
        tcnums.append(re.findall(r'\d+', ms)[-1])
    tcnums = set(tcnums)

    # prepare concatenated time chunks (TC) - SB_group#_TC#.MS:CALCOR_DATA_CIRC -> group#_TC#.MS:DATA (cal corr data, beam corrected, circular)
    logging.info('Concatenating timechunks...')
    for tcnum in tcnums:
        group_tc = sorted(glob.glob('group'+str(i)+'/*_TC'+tcnum+'.MS'))
        s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_timesplit/NDPPP-concat.parset msin="['+','.join(group_tc)+']"  msout=group'+str(i)+'/group'+str(i)+'_TC'+tcnum+'.MS' \
                log='group'+str(i)+'/NDPPP_concat_TC'+tcnum+'.log', run_type='NDPPP')
    s.run(check=True)

    check_rm('group'+str(i)+'/*_group*_TC*.MS') # remove splitted files

logging.info("Done.")
