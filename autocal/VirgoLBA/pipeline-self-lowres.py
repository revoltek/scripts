#!/usr/bin/env python
# apply solution from the calibrator and then run selfcal, an initial model
# and initial solutions from a calibrator must be provided
# local dir must contain all the MSs, in touched in linear pol

# initial self-cal model
model = '/home/fdg/scripts/autocal/VirgoA/150702_LBA-VirgoA.model'
# globaldb produced by pipeline-init
globaldb = '../cals/globaldb'
# fake skymodel with pointing direction
fakeskymodel = '/home/fdg/scripts/autocal/VirgoA/virgo.fakemodel.skymodel'

##############################################################

import sys, os, glob, re
from lofar import bdsm
import numpy as np
import lsmtool
import pyrap.tables as pt
from lib_pipeline import *
from make_mask import make_mask

set_logger()
s = Scheduler(qsub=False, max_threads=25, dry=False)

#################################################
# Clear
logging.info('Cleaning...')
check_rm('*log')
check_rm('*last')
check_rm('*h5')
check_rm('concat*')
check_rm('img')
os.makedirs('img')
check_rm('plot*')
check_rm('*avg.MS')

# all MS
mss = sorted(glob.glob('*.MS'))

##############################################
# Initial processing
#logging.info('Fix beam table...')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)

#################################################
# Copy cal solution
#logging.info('Copy solutions...')
#for ms in mss:
#    num = re.findall(r'\d+', ms)[-1]
#    logging.debug(globaldb+'/sol000_instrument-'+str(num)+' -> '+ms+'/instrument')
#    check_rm(ms+'/instrument')
#    os.system('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')

#########################################################################################
# [PARALLEL] apply solutions and beam correction - SB.MS:DATA -> SB.MS:CALCOR_DATA (calibrator corrected data, beam applied, linear)
#logging.info('Correcting target MSs...')
#for ms in mss:
#    s.add('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/VirgoA/parset_self-lowres/bbs-corbeam.parset '+fakeskymodel, \
#          log=ms+'-init_corbeam.log', cmd_type='BBS')
#s.run(check=True)

########################################################################################
# [PARALLEL] average and remove remote stations SB.MS:CALCOR_DATA -> SB-avg.MS:DATA
logging.info('Avg and remove remotes...')
for ms in mss:
    s.add('NDPPP /home/fdg/scripts/autocal/VirgoA/parset_self-lowres/NDPPP-concatavg1.parset msin='+ms+'\
            msin.datacolumn=CALCOR_DATA msout='+ms.replace('.MS','-avg.MS'), \
       log=ms+'-init_avg.log', cmd_type='NDPPP')
s.run(check=True)

# all MS
mss = sorted(glob.glob('*avg.MS'))

#########################################################################################
# [PARALLEL] Transform to circular pol - SB.MS:DATA -> SB-circ.MS:CIRC_DATA (data, beam applied, circular)
logging.info('Convert to circular...')
for ms in mss:
    s.add('/home/fdg/scripts/mslin2circ.py -i '+ms+':DATA -o '+ms+':CIRC_DATA', log=ms+'-init_circ2lin.log', cmd_type='python')
s.run(check=True)

#####################################################################################
# ft model, model is unpolarized CIRC == LIN - SB.MS:MODEL_DATA (best m87 model)
logging.info('Add models...')
for ms in mss:
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':ms, 'model':model}, log=ms+'_ft.log')
    s.run(check=True) # not parallel

#####################################################################################
# [PARALLEL] calibrate - SB.MS:CIRC_DATA (no correction)
logging.info('Calibrate...')
for ms in mss:
    s.add('NDPPP /home/fdg/scripts/autocal/VirgoA/parset_self-lowres/NDPPP-selfcal_modeldata.parset msin='+ms+' cal.parmdb='+ms+'/instrument', \
          log=ms+'_selfcal.log', cmd_type='NDPPP')
s.run(check=True)

#######################################################################################
# Solution rescaling
logging.info('Running LoSoTo to normalize solutions...')
os.makedirs('plot')
check_rm('globaldb')
os.makedirs('globaldb')
for num, ms in enumerate(mss):
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
    if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
h5parm = 'global.h5'

s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto.log', cmd_type='python')
s.run(check=False)
s.add('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/VirgoA/parset_self-lowres/losoto.parset', log='losoto.log', log_append=True, cmd_type='python')
s.run(check=False)
s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto.log', log_append=True, cmd_type='python')
s.run(check=True)

for num, ms in enumerate(mss):
    check_rm(ms+'/instrument')
    os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')

########################################################################################
# [PARALLEL] correct - SB.MS:CIRC_DATA -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
logging.info('Correct...')
for ms in mss:
    s.add('NDPPP /home/fdg/scripts/autocal/VirgoA/parset_self-lowres/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument', \
          log=ms+'_selfcor.log', cmd_type='NDPPP')
s.run(check=True)

###########################################################################################################################
# avg time only - SB.MS:CORRECTED_DATA -> concat-avg.MS:DATA (selfcal corrected data, beam applied, circular)
logging.info('Average...')
check_rm('concat-avg.MS*')
s.add('NDPPP /home/fdg/scripts/autocal/VirgoA/parset_self-lowres/NDPPP-concatavg2.parset msin="['+','.join(mss)+']" msout=concat-avg.MS avg.freqstep=1 avg.timestep=4', \
        log='concatavg.log', cmd_type='NDPPP')
s.run(check=True)

#########################################################################################################
# large FoV image
logging.info('Low-res wide field image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean-lowres.py', \
        params={'msfile':'concat-avg.MS', 'imagename':'img/clean-wide', 'imtype':'wide'}, log='clean-wide.log')
s.run(check=True)
make_mask(image_name = 'img/clean-wide.image.tt0', mask_name = 'img/wide.newmask')
logging.info('Low-res wide field image (with mask)...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean-lowres.py', \
        params={'msfile':'concat-avg.MS', 'imagename':'img/clean-wide-masked', 'imtype':'wide', 'mask':'img/wide.newmask'}, log='clean-wide-masked.log')

#########################################################################################################
# central clean
logging.info('Clean...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean-lowres.py', \
        params={'msfile':'concat-avg.MS', 'imagename':'img/normal'}, log='clean.log')
s.run(check=True)

##########################################################################################################
# uvsub + large FoV image
s.add('taql "update concat-avg.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub
s.run(check=False)

logging.info('Low-res wide field image sub...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean-lowres.py', \
        params={'msfile':'concat-avg.MS', 'imagename':'img/clean-wide-sub', 'imtype':'wide'}, log='clean-wide-sub.log')
s.run(check=True)
make_mask(image_name = 'img/clean-wide-sub.image.tt0', mask_name = 'img/wide-sub.newmask')
logging.info('Low-res wide field image sub (with mask)...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean-lowres.py', \
        params={'msfile':'concat-avg.MS', 'imagename':'img/clean-wide-sub-masked', 'imtype':'wide', 'mask':'img/wide.newmask'}, log='clean-wide-sub-masked.log')

logging.info("Done.")
