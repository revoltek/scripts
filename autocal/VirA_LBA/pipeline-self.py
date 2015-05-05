#!/usr/bin/env python
# apply solution from the calibrator and then run selfcal, an initial model
# and initial solutions from a calibrator must be provided
# local dir must contain all the MSs, un touched in linear pol
# 1 selfcal
# 2 subtract field
# 3 flag

# initial self-cal model
model = '/home/fdg/scripts/autocal/VirA_LBA/150328_LBA-VirgoA.model'
# globaldb produced by pipeline-init
globaldb = '../cals/globaldb'
# fake skymodel with pointing direction
fakeskymodel = '/home/fdg/scripts/autocal/VirA_LBA/virgo.fakemodel.skymodel'
# SB per block
n = 15

##############################################################

import sys, os, glob, re
from lofar import bdsm
import numpy as np
import lsmtool
import pyrap.tables as pt
from lib_pipeline import *
from make_mask import make_mask

set_logger()
s = Scheduler(qsub=False, max_threads=20, dry=False)

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
#logging.info('Fix beam table...')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)

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
#logging.info('Correcting target MSs...')
#for ms in mss:
#    s.add('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/VirA_LBA/parset_self/bbs-corbeam.parset '+fakeskymodel, \
#          log=ms+'-init_corbeam.log', cmd_type='BBS')
#s.run(check=True)

#########################################################################################
# [PARALLEL] Transform to circular pol - SB.MS:CALCOR_DATA -> SB-circ.MS:CIRC_DATA (data, beam applied, circular)
#logging.info('Convert to circular...')
#for ms in mss:
#    s.add('/home/fdg/scripts/mslin2circ.py -i '+ms+':CALCOR_DATA -o '+ms+':CIRC_DATA', log=ms+'-init_circ2lin.log', cmd_type='python')
#s.run(check=True)
 
# self-cal cycle -> 5
for i in xrange(5):
    logging.info('Starting self-cal cycle: '+str(i))

    # MS for calibration, use all at cycle 0 and 4
    if i == 0 or i == 4:
        mss_c = mss
        mss_clean = mss[::n]
    else:
        mss_c = mss[::n]
        mss_clean = mss[::n]

    if i != 0:    
        model = 'img/clean-c'+str(i-1)+'.model'

    #####################################################################################
    # ft model, model is unpolarized CIRC == LIN - SB.MS:MODEL_DATA (best m87 model)
    # TODO: test if adding wprojplanes here improves the calibration
    logging.info('Add models...')
    check_rm('concat1.MS*')
    pt.msutil.msconcat(mss_c[:len(mss_c)/2], 'concat1.MS', concatTime=False)
    check_rm('concat2.MS*')
    pt.msutil.msconcat(mss_c[len(mss_c)/2:], 'concat2.MS', concatTime=False)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat1.MS', 'model':model, 'wproj':512}, log='ft-virgo-c'+str(i)+'.log')
    s.run(check=True)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat2.MS', 'model':model, 'wproj':512}, log='ft-virgo-c'+str(i)+'.log')
    s.run(check=True)
    #for ms in mss_c:
    #    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':ms, 'model':model, 'wproj':512}, log=ms+'_ft-virgo-c'+str(i)+'.log')
    #    s.run(check=True)

    #####################################################################################
    # [PARALLEL] calibrate - SB.MS:CIRC_DATA (no correction)
    # TODO: calibrated on wide-field subtracted data?
    logging.info('Calibrate...')
    for ms in mss_c:
        s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self/NDPPP-selfcal_modeldata.parset msin='+ms+' cal.parmdb='+ms+'/instrument', \
              log=ms+'_selfcal-c'+str(i)+'.log', cmd_type='NDPPP')
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
    s.add('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/VirA_LBA/parset_self/losoto.parset', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    for num, ms in enumerate(mss_c):
        check_rm(ms+'/instrument')
        os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
    os.system('mv plot plot-c'+str(i))

    ########################################################################################
    # [PARALLEL] correct - SB.MS:CIRC_DATA -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
    logging.info('Correct...')
    for ms in mss_c:
        s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument', \
              log=ms+'_selfcor-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

#######################
#   QUICK TEST LOOP
#    # avg - SB.MS:CORRECTED_DATA -> concat-avg.MS:DATA
#    check_rm('concat-avg.MS*')
#    logging.info('Average...')
#    s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self/NDPPP-concatavg.parset msin="['+','.join(mss_clean)+']" msout=concat-avg.MS', log='concatavg-c'+str(i)+'.log', cmd_type='NDPPP')
#    s.run(check=True)
#    # clean (make a new model of virgo)
#    logging.info('Clean...')
#    s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', params={'msfile':'concat-avg.MS', 'imagename':'img/clean-c'+str(i)}, log='clean-c'+str(i)+'.log')
#    s.run(check=True)
#    continue
######################

    # create widefield model
    if i == 0 or i == 4:

        # concatenate data - MS.MS -> concat.MS (selfcal corrected data, beam applied, circular)
        # NOTE: split concat into two parts should avoid problems
        logging.info('Make widefield model - Concatenating...')
        check_rm('concat1.MS*')
        pt.msutil.msconcat(mss_c[:len(mss_c)/2], 'concat1.MS', concatTime=False)
        check_rm('concat2.MS*')
        pt.msutil.msconcat(mss_c[len(mss_c)/2:], 'concat2.MS', concatTime=False)

        # uvsub, MODEL_DATA is still Virgo
        logging.info('Make widefield model - UV-Subtracting Virgo A...')
        os.system('taql "update concat1.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub
        os.system('taql "update concat2.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub

        # clean, mask, clean
        logging.info('Make widefield model - Widefield imaging...')
        check_rm('concat.MS*')
        pt.msutil.msconcat(mss_clean, 'concat.MS', concatTime=False)
        imagename = 'img/clean-wide-c'+str(i)
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
                params={'msfile':'concat.MS', 'imagename':imagename, 'imtype':'wide'}, log='clean-wide1-c'+str(i)+'.log')
        s.run(check=True)
        make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask')
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/VirA_LBA/m87.crtf'})
        s.run(check=True)
        logging.info('Make widefield model - Widefield imaging2...')
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
                params={'msfile':'concat.MS', 'imagename':imagename.reaplce('wide','wide-masked'), 'mask':imagename+'.newmask', 'imtype':'wide'}, log='clean-wide2-c'+str(i)+'.log')
        s.run(check=True)

        # Subtract widefield model using ft on a virtual concat - concat.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA-MODEL_DATA (selfcal corrected data, beam applied, circular, field sources subtracted)
        logging.info('Flagging - Subtracting wide-field model...')
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', \
                    params={'msfile':'concat1.MS', 'model':'img/wide-masked-c'+str(i)+'.model', 'wproj':512}, log='ft-flag-c'+str(i)+'.log')
        s.run(check=True)
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', \
                    params={'msfile':'concat2.MS', 'model':'img/wide-masked-c'+str(i)+'.model', 'wproj':512}, log='ft-flag-c'+str(i)+'.log')
        s.run(check=True)
        os.system('taql "update concat1.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub
        os.system('taql "update concat2.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub

        logging.info('Flagging - Flagging residuals...')
        run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_flag.py', params={'msfile':'concat1.MS'}, log='flag-c'+str(i)+'.log')
        s.run(check=True)
        run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_flag.py', params={'msfile':'concat2.MS'}, log='flag-c'+str(i)+'.log')
        s.run(check=True)

        # [PARALLEL] reapply NDPPP solutions - SB.MS:CIRC_DATA -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
        # this because with the virtual concat the CORRECTED_DATA have been uvsubbed
        logging.info('Make widefield model - Re-correct...')
        for ms in mss_c:
            s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument', \
                    log=ms+'_selfcor-c'+str(i)+'.log', cmd_type='NDPPP')
        s.run(check=True)

    # Subtract widefield model using ft on a virtual concat - concat.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA-MODEL_DATA (selfcal corrected data, beam applied, circular, field sources subtracted)
    logging.info('Subtracting wide-field model...')
    check_rm('concat.MS*')
    pt.msutil.msconcat(mss_clean, 'concat.MS', concatTime=False)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', \
            params={'msfile':'concat.MS', 'model':'img/wide-c'+str(i)+'.model', 'wproj':512}, log='ft-wide-c'+str(i)+'.log')
    s.run(check=True)
    os.system('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub

    # avg 1chanSB/20s - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, circular)
    logging.info('Average...')
    check_rm('concat.MS*')
    s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self/NDPPP-concatavg.parset msin="['+','.join(mss_clean)+']" msout=concat.MS', \
            log='concatavg-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # clean (make a new model of virgo)
    logging.info('Clean...')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
            params={'msfile':'concat.MS', 'imagename':'img/clean-c'+str(i)}, log='clean-c'+str(i)+'.log')
    s.run(check=True)

##########################################################################################################
# [PARALLEL] concat+avg - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, circ)
logging.info('Concat...')
check_rm('concat*')
s.add('NDPPP /home/fdg/scripts/autocal/PerA_LBA/parset_self/NDPPP-concatavg.parset msin="['+','.join(mss[len(mss_c)/2,:])+']" msout=concat1.MS', \
                log='final_concatavg1.log', cmd_type='NDPPP')
s.add('NDPPP /home/fdg/scripts/autocal/PerA_LBA/parset_self/NDPPP-concatavg.parset msin="['+','.join(mss[:,len(mss_c)/2])+']" msout=concat2.MS', \
                log='final_concatavg2.log', cmd_type='NDPPP')
s.run(check=True)
s.add('NDPPP /home/fdg/scripts/autocal/PerA_LBA/parset_self/NDPPP-concat.parset msin="concat1.MS,concat2.MS" msout=concat.MS', \
                log='final_concatavg.log', cmd_type='NDPPP')
s.run(check=True)

#########################################################################################################
# group images of VirA
logging.info('Full BW image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
        params={'msfile':'concat.MS', 'imagename':'img/clean-all'}, log='final_clean-all.log')
s.run(check=True)

#########################################################################################################
# low-res image
logging.info('Make low-resolution image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
        params={'msfile':'concat.MS', 'imagename':'img/clean-lr', 'imtype':'lr'}, log='final_clean-lr.log')
s.run(check=True)

##########################################################################################################
# uvsub + large FoV image
check_rm('concat.MS*')
pt.msutil.msconcat(mss, 'concat.MS', concatTime=False)
logging.info('Ft+uvsub of M87 model...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', \
        params={'msfile':'concat.MS', 'model':'img/clean-c'+str(i)+'.model'}, log='final_ft.log')
s.run(check=True)
os.system('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub

logging.info('Low-res wide field image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBA/casa_clean.py', \
        params={'msfile':'concat.MS', 'imagename':'img/clean-wide', 'imtype':'wide'}, log='final_clean-wide.log')
s.run(check=True)

logging.info("Done.")
