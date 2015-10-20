#!/usr/bin/python
# data preparation for selfcal, apply cal solutions and split SB in time and concatenate in freq
# Input:
# Virgin target MSs and a globaldb of the calibrator
# Output:
# set of group*_TC*.MS file with DATA = calibrator corrected data, beam corrupted

parset_dir = '/home/stsf309/scripts/autocal/1RXSJ0603_LBA/parset_timesplit'
ngroups = 12 # number of groups (totalSB/SBperFREQgroup) ~ 20 SB/group
initc = 0 # initial tc num (useful for multiple observation of same target) - tooth10==12
fakeskymodel = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/toothbrush.fakemodel.skymodel'
globaldb = 'globaldb-blavg'

##################################################################################################

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from lib_pipeline import *

set_logger()
s = Scheduler(qsub=True, max_threads=24, dry=False)

##################################################
## Clear
#logging.info('Cleaning...')
#check_rm('*log')
#check_rm('*group*')
#
mss = sorted(glob.glob('*MS'))
#
###############################################
## Initial processing
#logging.info('Fix beam table')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)
#
##################################################
## Copy cal solution
#for ms in mss:
#    num = re.findall(r'\d+', ms)[-1]
#    check_rm(ms+'/instrument')
#    logging.debug('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')
#    os.system('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')
#
## Apply cal sol - SB.MS:DATA -> SB.MS:CALCOR_DATA (calibrator corrected data, beam corrected, linear)
#logging.info('Apply solutions...')
#for ms in mss:
#    s.add('calibrate-stand-alone --replace-sourcedb '+ms+' '+parset_dir+'/bbs_correct.parset '+fakeskymodel, log=ms+'_cor.log', cmd_type='BBS')
#s.run(check=True)
#
####################################################################################################
## To circular - SB.MS:CALCOR_DATA -> SB.MS:CALCOR_DATA_CIRC (calibrator corrected data, beam corrected, circular)
#logging.info('Convert to circular...')
#for ms in mss:
#    s.add('/home/fdg/scripts/mslin2circ.py -s -i '+ms+':CALCOR_DATA -o '+ms+':CALCOR_DATA_CIRC', log=ms+'_circ2lin.log', cmd_type='python')
#s.run(check=True)

###################################################################################################
# split each MS in timechunks of 1 h and create groups
# TODO: add possibility to create uneaven groups (e.g. more SB at low-freq)
groupnames = []
logging.info('Concatenating in frequency...')
for i, msg in enumerate(np.array_split(mss, ngroups)):
    groupname = 'group%02i' %i
    groupnames.append(groupname)
#    check_rm(groupname)
#    os.system('mkdir '+groupname)

    # prepare concatenated time chunks (TC) - SB.MS:CALCOR_DATA_CIRC -> group#.MS:DATA (cal corr data, beam corrected, circular)
#    s.add('NDPPP '+parset_dir+'/NDPPP-concat.parset msin="['+','.join(msg)+']"  msout='+groupname+'/'+groupname+'.MS', \
#                log=groupname+'/NDPPP_concat.log', cmd_type='NDPPP')
#s.run(check=True)

# Flagging on concatenated dataset
#logging.info('Flagging...')
#for groupname in groupnames:
#    s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+groupname+'/'+groupname+'.MS', \
#                log=groupname+'/NDPPP_flag.log', cmd_type='NDPPP')
#s.run(check=True)

for groupname in groupnames:
    ms = groupname+'/'+groupname+'.MS'
    if not os.path.exists(ms): continue
    t = pt.table(ms, ack=False)
    starttime = t[0]['TIME']
    endtime   = t[t.nrows()-1]['TIME']
    hours = (endtime-starttime)/3600.
    logging.debug(ms+' has length of '+str(hours)+' h.')

    # split this ms into many TCs (2*#hours)
    # to re-concat:
    #   t = table(['T0','T1',...])
    #   t.sort('TIME').copy('output.MS', deep = True)
    for tc, timerange in enumerate(np.array_split(t.getcol('TIME'), 2*round(hours))):
        tc += initc
        logging.debug('Splitting timerange '+str(timerange[0])+' - '+str(timerange[-1]))
        t1 = t.query('TIME >= ' + str(timerange[0]) + ' && TIME <= ' + str(timerange[-1]), sortlist='TIME,ANTENNA1,ANTENNA2')
        splitms = ms.replace('.MS','_TC%02i.MS' % tc)
        check_rm(splitms)
        t1.copy(splitms, True)
        t1.close()
    t.close()

    check_rm(ms) # remove not-timesplitted file

logging.info("Done.")
