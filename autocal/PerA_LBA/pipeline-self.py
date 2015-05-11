#!/usr/bin/env python
# apply solution from the calibrator and then run selfcal, an initial model
# and initial solutions from a calibrator must be provided
# local dir must contain all the MSs, un touched in linear pol
# 1 selfcal
# 2 subtract field
# 3 flag

# initial self-cal model
model = '/home/fdg/scripts/autocal/PerA_LBA/150510_LBA-PerseusA.model'
# globaldb produced by pipeline-init
globaldb = '../cals/globaldb'
# fake skymodel with pointing direction
fakeskymodel = '/home/fdg/scripts/autocal/PerA_LBA/perseus.fakemodel.skymodel'

##############################################################

import sys, os, glob, re
from lofar import bdsm
import numpy as np
import lsmtool
import pyrap.tables as pt
from lib_pipeline import *
from make_mask import make_mask

set_logger()
s = Scheduler(qsub=False, max_threads=23, dry=False)

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

# all MS
mss = sorted(glob.glob('*.MS'))

##############################################
# Initial processing
logging.info('Fix beam table...')
for ms in mss:
    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
s.run(check=False)

#################################################
# Copy cal solution
logging.info('Copy solutions...')
for ms in mss:
    num = re.findall(r'\d+', ms)[-1]
    logging.debug(globaldb+'/sol000_instrument-'+str(num)+' -> '+ms+'/instrument')
    check_rm(ms+'/instrument')
    os.system('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')

#########################################################################################
# [PARALLEL] apply solutions and beam correction - SB.MS:DATA -> SB.MS:CALCOR_DATA (calibrator corrected data, beam applied, linear)
logging.info('Correcting target MSs...')
for ms in mss:
    s.add('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/PerA_LBA/parset_self/bbs-corbeam.parset '+fakeskymodel, \
          log=ms+'-init_corbeam.log', cmd_type='BBS')
s.run(check=True)

# [PARALLEL] apply solutions and beam correction - SB.MS:MODEL_DATA (PerA model, linear==circular, with ARRAY_FACTOR to simulate apparent sky)
#logging.info('Filling MODEL_DATA...')
#for ms in mss:
#    s.add('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/PerA_LBA/parset_self/bbs-pre.parset '+model, \
#          log=ms+'-init_predict.log', cmd_type='BBS')
#s.run(check=True)

#########################################################################################
# [PARALLEL] Transform to circular pol - SB.MS:CALCOR_DATA -> SB-circ.MS:CIRC_DATA (data, beam applied, circular)
logging.info('Convert to circular...')
for ms in mss:
    s.add('/home/fdg/scripts/mslin2circ.py -i '+ms+':CALCOR_DATA -o '+ms+':CIRC_DATA', log=ms+'-init_circ2lin.log', cmd_type='python')
s.run(check=True)
 
# self-cal cycle
for i in xrange(25):
    logging.info('Starting self-cal cycle: '+str(i))

    if i != 0:    
        model = 'img/clean-c'+str(i-1)+'.model'

    #####################################################################################
    # ft model, model is unpolarized CIRC == LIN - SB.MS:MODEL_DATA (best m87 model)
    logging.info('Add models...')
    for ms in mss:
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':ms, 'model':model, 'wproj':512}, log=ms+'_ft-perseus-c'+str(i)+'.log')
        s.run(check=True)

    #####################################################################################
    # [PARALLEL] calibrate - SB.MS:CIRC_DATA (no correction)
    logging.info('Calibrate...')
    for ms in mss:
        s.add('NDPPP /home/fdg/scripts/autocal/PerA_LBA/parset_self/NDPPP-selfcal_modeldata.parset msin='+ms+' cal.parmdb='+ms+'/instrument', \
              log=ms+'_selfcal-c'+str(i)+'.log', cmd_type='NDPPP')
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
    h5parm = 'global-c'+str(i)+'.h5'

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(i)+'.log', cmd_type='python')
    s.run(check=False)
    s.add('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/PerA_LBA/parset_self/losoto.parset', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    for num, ms in enumerate(mss):
        check_rm(ms+'/instrument')
        os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
    os.system('mv plot plot-c'+str(i))

    ########################################################################################
    # [PARALLEL] correct - SB.MS:CIRC_DATA -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
    logging.info('Correct...')
    for ms in mss:
        s.add('NDPPP /home/fdg/scripts/autocal/PerA_LBA/parset_self/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument', \
              log=ms+'_selfcor-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

#######################
#   QUICK TEST LOOP
#    # avg - SB.MS:CORRECTED_DATA -> concat-avg.MS:DATA
#    check_rm('concat-avg.MS*')
#    logging.info('Average...')
#    s.add('NDPPP /home/fdg/scripts/autocal/PerA_LBA/parset_self/NDPPP-concatavg.parset msin="['+','.join(mss)+']" msout=concat-avg.MS', log='concatavg-c'+str(i)+'.log', cmd_type='NDPPP')
#    s.run(check=True)
#    # clean (make a new model of perseus)
#    logging.info('Clean...')
#    imagename = 'img/clean-c'+str(i)
#    s.add_casa('/home/fdg/scripts/autocal/casa_comm/perseusLBA/casa_clean.py', params={'msfile':'concat-avg.MS', 'imagename':imagename}, log='clean-c'+str(i)+'.log')
#    s.run(check=True)
#    make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask')
#    s.add_casa('/home/fdg/scripts/autocal/casa_comm/perseusLBA/casa_clean.py', params={'msfile':'concat-avg.MS', 'imagename':imagename+'-masked', 'mask':imagename+'.newmask'}, log='clean2-c'+str(i)+'.log')
#    s.run(check=True)
#    continue
######################

    ############################################################################################
    # flag residual
    if i % 5 == 0:

        # add SUBTRACTED_DATA the first run
        if i == 0:
            logging.info('Flag - add SUBTRACTED_DATA column...')
            for ms in mss:
                s.add('addcol2ms.py -i '+ms+' -o SUBTRACTED_DATA', log=ms+'_flag-addcol.log', cmd_type='python')
            s.run(check=True)

        # uvsub, MODEL_DATA is still PerA
        logging.info('Flag - UV-Subtracting Per A...')
        for ms in mss:
            os.system('taql "update '+ms+' set SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub

        # Madflagger
        logging.info('Flag - Flagging residuals...')
        for ms in mss:
            s.add('NDPPP /home/fdg/scripts/autocal/PerA_LBA/parset_self/NDPPP-flag.parset msin='+ms, \
                    log=ms+'_flag-c'+str(i)+'.log', cmd_type='NDPPP')
        s.run(check=True)

    ##################################################################################################
    # imaging
    # avg 1chanSB/20s - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, circular)
    logging.info('Average...')
    check_rm('concat-avg.MS*')
    s.add('NDPPP /home/fdg/scripts/autocal/PerA_LBA/parset_self/NDPPP-concatavg.parset msin="['+','.join(mss)+']" msout=concat-avg.MS', \
            log='concatavg-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # clean (make a new model of perseus)
    logging.info('Clean...')
    imagename = 'img/clean-c'+str(i)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/perseusLBA/casa_clean.py', \
            params={'msfile':'concat-avg.MS', 'imagename':imagename}, log='clean-c'+str(i)+'.log')
    s.run(check=True)
    make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/PerA_LBA/PerA.crtf', 'setTo':1}, log='blank-c'+str(i)+'.log')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/perseusLBA/casa_clean.py', params={'msfile':'concat-avg.MS', 'imagename':imagename+'-masked', 'mask':imagename+'.newmask'}, log='clean2-c'+str(i)+'.log')
    s.run(check=True)

##########################################################################################################
# [PARALLEL] concat+avg - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, circ)
logging.info('Concat...')
check_rm('concat.MS*')
s.add('NDPPP /home/fdg/scripts/autocal/PerA_LBA/parset_self/NDPPP-concatavg.parset msin="['+','.join(mss)+']" msout=concat.MS', \
        log='final_concatavg.log', cmd_type='NDPPP')
s.run(check=True)

#########################################################################################################
# group images of PerA
logging.info('Full BW image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/perseusLBA/casa_clean.py', \
        params={'msfile':'concat.MS', 'imagename':'img/clean-all'}, log='final_clean-all.log')
s.run(check=True)

#########################################################################################################
# low-res image
logging.info('Make low-resolution image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/perseusLBA/casa_clean.py', \
        params={'msfile':'concat.MS', 'imagename':'img/clean-lr', 'imtype':'lr'}, log='final_clean-lr.log')
s.run(check=True)

##########################################################################################################
# uvsub + large FoV image
check_rm('concat.MS*')
pt.msutil.msconcat(mss, 'concat.MS', concatTime=False)
logging.info('Ft+uvsub of Perseus A model...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', \
        params={'msfile':'concat.MS', 'model':'img/clean-c'+str(i)+'.model'}, log='final_ft.log')
s.run(check=True)
os.system('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub

logging.info('Low-res wide field image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/perseusLBA/casa_clean.py', \
        params={'msfile':'concat.MS', 'imagename':'img/clean-wide', 'imtype':'wide'}, log='final_clean-wide.log')
s.run(check=True)

logging.info("Done.")
