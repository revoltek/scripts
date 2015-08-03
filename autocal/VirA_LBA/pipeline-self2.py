#!/usr/bin/env python
# apply solution from the calibrator and then run selfcal, an initial model
# and initial solutions from a calibrator must be provided
# local dir must contain all the MSs, un touched in linear pol
# 1 selfcal
# 2 subtract field
# 3 flag

# initial self-cal model
model = '/home/fdg/scripts/autocal/VirA_LBA/150328_LBA-VirgoA.model.tt0'
# globaldb produced by pipeline-init
globaldb = '../cals/globaldb'
# fake skymodel with pointing direction
fakeskymodel = '/home/fdg/scripts/autocal/VirA_LBA/virgo.fakemodel.skymodel'
# SB to skip for the initial selfcal
n = 10
# number of selfcal cycles
cycles = 10

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

# all MS
mss = sorted(glob.glob('*.MS'))
# subset of MSs for selfcal
mss_c = mss[::n]

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
logging.info('Correcting target MSs...')
for ms in mss:
    s.add('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/VirA_LBA/parset_self2/bbs-corbeam.parset '+fakeskymodel, \
          log=ms+'-init_corbeam.log', cmd_type='BBS')
s.run(check=True)

#########################################################################################
# [PARALLEL] Transform to circular pol - SB.MS:CALCOR_DATA -> SB-circ.MS:CIRC_DATA (data, beam applied, circular)
logging.info('Convert to circular...')
for ms in mss:
    s.add('/home/fdg/scripts/mslin2circ.py -i '+ms+':CALCOR_DATA -o '+ms+':CIRC_DATA', log=ms+'-init_circ2lin.log', cmd_type='python')
s.run(check=True)

########################################################################################
# [PARALLEL] Initialize columns - SB.MS:CIRC_DATA_SUB = CIRC_DATA
logging.info('Make new columns...')
for ms in mss:
    s.add('addcol2ms.py -i '+ms+' -o CIRC_DATA_SUB', log=ms+'-init_addcol.log', cmd_type='python')
s.run(check=True)
logging.info('Reset CIRC_DATA_SUB...')
for ms in mss:
    s.add('taql "update '+ms+' set CIRC_DATA_SUB = CIRC_DATA"')
s.run(check=True)
# After all columns are created
logging.info('Concat...')
pt.msutil.msconcat(mss_c, 'concat.MS', concatTime=False)
 
# self-cal cycle -> 5
for i in xrange(cycles):
    logging.info('Starting self-cal cycle: '+str(i))

    # first cycle use given model
    if i != 0:
        model = 'img/clean-c'+str(i-1)+'.model'

    # last cycle do on all SBs
    if i == cycles-1:
        mss_c = mss
        logging.info('Concat...')
        check_rm('concat.MS*')
        pt.msutil.msconcat(mss_c, 'concat.MS', concatTime=False)

    #####################################################################################
    # ft model, model is unpolarized CIRC == LIN - SB.MS:MODEL_DATA (best m87 model)
    logging.info('Add models...')
    if i == cycles-1:
        for ms in mss_c:
            s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':ms, 'model':model}, log=ms+'_final-ft-virgo.log')
            s.run(check=True) # not parallel!
    else:
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':model}, log='ft-virgo-c'+str(i)+'.log')
        s.run(check=True)

    #####################################################################################
    # [PARALLEL] calibrate - SB.MS:CIRC_DATA_SUB (no correction)
    logging.info('Calibrate...')
    for ms in mss_c:
        s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self2/NDPPP-selfcal_modeldata.parset msin='+ms+' cal.parmdb='+ms+'/instrument', \
              log=ms+'_selfcal-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    #############################################################################################
    # create widefield model
    if i%3 == 0 or i == cycles-1:

        # [PARALLEL] apply NDPPP solutions on complete dataset - SB.MS:CIRC_DATA -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
        logging.info('Make widefield model - Correct...')
        for ms in mss_c:
            s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self2/NDPPP-selfcor.parset msin='+ms+' msin.datacolumn=CIRC_DATA cor.parmdb='+ms+'/instrument', \
                    log=ms+'_flag-selfcor-c'+str(i)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # uvsub, MODEL_DATA is still Virgo
        logging.info('Make widefield model - UV-Subtracting Virgo A...')
        s.add('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub
        s.run(check=False)

        # clean, mask, clean
        if i != cycles-1:
            logging.info('Make widefield model - Widefield imaging...')
            imagename = 'img/clean-wide-c'+str(i)
            s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
                    params={'msfile':'concat.MS', 'imagename':imagename, 'imtype':'wide'}, log='clean-wide1-c'+str(i)+'.log')
            s.run(check=True)
            logging.info('Make widefield model - Make mask...')
            make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask')
            s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/VirA_LBA/m87-blank.crtf'}, log='blank-c'+str(i)+'.log')
            s.run(check=True)
            logging.info('Make widefield model - Widefield imaging2...')
            s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
                    params={'msfile':'concat.MS', 'imagename':imagename.replace('wide','wide-masked'), 'mask':imagename+'.newmask', 'imtype':'widemasked'}, log='clean-wide2-c'+str(i)+'.log')
            s.run(check=True)
            widemodel = imagename.replace('wide','wide-masked')+'.model'

        ###############################################################################################################################
        # ft widefield model
        logging.info('Make widefield model - ft() widefield model...')
        if i == cycles-1:
            for ms in mss_c:
                s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':ms, 'model':widemodel, 'wproj':512}, log=ms+'_flag-ft-c'+str(i)+'.log')
                s.run(check=True) # not parallel!
        else:
            s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', \
                        params={'msfile':'concat.MS', 'model':widemodel, 'wproj':512}, log='flag-ft-c'+str(i)+'.log')
            s.run(check=True)

        # subtract widefield model - concat.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA=CORRECTED_DATA-MODEL_DATA (selfcal corrected data, beam applied, circular, field sources subtracted)
        logging.info('Make widefield model - Subtract widefield model...')
        s.run('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub
        s.run(check=False)

        ########################################################################################################################
        # [PARALLEL] Flagging on CORRECTED_DATA
        logging.info('Make widefield model - Flagging residuals...')
        for ms in mss_c:
            s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self2/NDPPP-flag.parset msin='+ms, \
                    log=ms+'_flag-c'+str(i)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        ########################################################################################################################
        # [PARALLEL] subtract widefield model SB.MS:CIRC_DATA -> SB.MS:CIRC_DATA_SUB (uncal data with subtracted the widefield model)
        # if last cycle skip (useless), but if one but last cycle, do on all SBs
        if i != cycles-1:
            if i >= cycles-4:
                logging.info('Make widefield model - subtract from uncal data on all MSs...')
                this_mss = mss
            else:
                logging.info('Make widefield model - subtract from uncal data...')
                this_mss = mss_c

            for ms in this_mss:
                s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/VirA_LBA/parset_self2/bbs-subcorpt.parset', \
                    log=ms+'_subcorpt-c'+str(i)+'.log', cmd_type='BBS')
            s.run(check=True)

    #######################################################################################
    # Solution rescaling
    logging.info('Running LoSoTo to normalize solutions...')
    os.makedirs('plot')
    check_rm('globaldb')
    os.makedirs('globaldb')
    for num, ms in enumerate(mss_c):
        os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
        if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
    h5parm = 'global-c'+str(i)+'.h5'

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(i)+'.log', cmd_type='python')
    s.run(check=False)
    s.add('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/VirA_LBA/parset_self2/losoto.parset', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    for num, ms in enumerate(mss_c):
        check_rm(ms+'/instrument')
        os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
    os.system('mv plot plot-c'+str(i))

    ########################################################################################
    # [PARALLEL] correct - SB.MS:CIRC_DATA_SUB -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
    logging.info('Correct...')
    for ms in mss_c:
        s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self2/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument', \
              log=ms+'_selfcor-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    ###########################################################################################################################
    # avg 1chanSB/20s - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, circular)
    logging.info('Average...')
    check_rm('concat-avg.MS*')
    s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self2/NDPPP-concatavg.parset msin="['+','.join(mss_c)+']" msout=concat-avg.MS', \
            log='concatavg-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # clean (make a new model of virgo)
    logging.info('Clean...')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
            params={'msfile':'concat-avg.MS', 'imagename':'img/clean-c'+str(i)}, log='clean-c'+str(i)+'.log')
    s.run(check=True)

#########################################################################################################
# low-res image
logging.info('Make low-resolution image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
        params={'msfile':'concat-avg.MS', 'imagename':'img/clean-lr', 'imtype':'lr'}, log='final_clean-lr.log')
s.run(check=True)

##########################################################################################################
# uvsub + large FoV image
logging.info('Ft+uvsub of M87 model...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', \
        params={'msfile':'concat-avg.MS', 'model':'img/clean-c'+str(i)+'.model'}, log='final_ft.log')
s.run(check=True)
s.add('taql "update concat-avg.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub
s.run(check=False)

logging.info('Low-res wide field image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
        params={'msfile':'concat-avg.MS', 'imagename':'img/clean-wide', 'imtype':'wide'}, log='final_clean-wide.log')
s.run(check=True)

logging.info("Done.")
