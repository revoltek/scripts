#!/usr/bin/env python
# apply solution from the calibrator and then run selfcal, an initial model
# and initial solutions from a calibrator must be provided
# local dir must contain all the MSs, un touched in linear pol
# 1 selfcal
# 2 subtract field
# 3 flag

# initial self-cal model
model = '/home/fdg/scripts/autocal/VirA_LBA/150328_LBA-VirgoA-highres.model'
# globaldb produced by pipeline-init
globaldb = '../cals/globaldb'
# fake skymodel with pointing direction
fakeskymodel = '/home/fdg/scripts/autocal/VirA_LBA/virgo.fakemodel.skymodel'
# SB to skip for the initial selfcal
n = 1
# number of selfcal cycles
cycles = 5

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
mss = sorted(glob.glob('VirAis_*.MS'))

#################################################
# Copy cal solution
#logging.info('Copy solutions...')
#for ms in mss:
#    num = re.findall(r'\d+', ms)[-1]
#    logging.debug(globaldb+'/sol000_instrument-'+str(num)+' -> '+ms+'/instrument')
#    check_rm(ms+'/instrument')
# NOTE: no cal for IS obs    
#    os.tystem('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')

#avg.freqstep = 4
#avg.timestep = 2
logging.info('Average...')
for ms in mss:
    msout = ms.replace('VirAis','VirAis-avg')
    check_rm(msout)
    s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self-is/NDPPP-concatavg.parset msin.datacolumn=DATA msin='+ms+' msout='+msout, log=ms+'-init_avg.log', cmd_type='NDPPP')
s.run(check=True)

mss = sorted(glob.glob('VirAis-avg_*.MS'))
mss_c = mss[::n]
 
#########################################################################################
# [PARALLEL] apply solutions and beam correction - SB.MS:DATA -> SB.MS:CALCOR_DATA (calibrator corrected data, beam applied, linear)
logging.info('Correcting target MSs...')
for ms in mss:
    s.add('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/VirA_LBA/parset_self-is/bbs-beam.parset '+fakeskymodel, \
          log=ms+'-init_beam.log', cmd_type='BBS')
s.run(check=True)

#########################################################################################
# [PARALLEL] Transform to circular pol - SB.MS:CALCOR_DATA -> SB-circ.MS:CIRC_DATA (data, beam applied, circular)
logging.info('Convert to circular...')
for ms in mss:
    s.add('/home/fdg/scripts/mslin2circ.py -i '+ms+':CALCOR_DATA -o '+ms+':CIRC_DATA', log=ms+'-init_circ2lin.log', cmd_type='python')
s.run(check=True)

# self-cal cycle -> 5
for i in xrange(cycles):
    logging.info('Starting self-cal cycle: '+str(i))

    # select uv-cut in meters
    if i == 0:
        blmin = 100000
        blmax = 1300000
    elif i == 1:
        blmin = 50000
        blmax = 1300000
    elif i == 2:
        blmin = 25000
        blmax = 1300000
    elif i == 3:
        blmin = 10000
        blmax = 1300000
    else:
        blmin = 0
        blmax = 1300000
    uvrange = str(blmin)+'~'+str(blmax)+'meters'

    # first cycle use given model
    if i != 0:
        model = 'img/clean-c'+str(i-1)+'.model'

    # last cycle do on all SBs
    if i == cycles-1:
        mss_c = mss

    #####################################################################################
    # ft model, model is unpolarized CIRC == LIN - SB.MS:MODEL_DATA (best m87 model)
    logging.info('Ft() ms per ms.')
    for ms in mss_c:
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':ms, 'model':model}, log=ms+'_ft-virgo-c'+str(i)+'.log')
        s.run(check=True) # not in parallel!

    #####################################################################################
    # [PARALLEL] calibrate - SB.MS:CIRC_DATA (no correction)
    logging.info('Calibrate...')
    for ms in mss_c:
        s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self-is/NDPPP-selfcal_modeldata.parset msin='+ms+' cal.parmdb='+ms+'/instrument'+' msin.blrange="['+str(blmin)+','+str(blmax)+']"', \
              log=ms+'_selfcal-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    #######################################################################################
    # flag
    if i%3 == 0:

        # [PARALLEL] correct - SB.MS:CIRC_DATA -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
        logging.info('Correct...')
        for ms in mss_c:
            s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self-is/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument', \
                  log=ms+'_selfcor-c'+str(i)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # uvsub, MODEL_DATA is still Virgo
        logging.info('Flagging - UV-Subtracting Virgo A...')
        check_rm('concat.MS*')
        pt.msutil.msconcat(mss_c, 'concat.MS', concatTime=False)
        s.add('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub
        s.run(check=False)

        logging.info('Flagging - Flagging residuals...')
        for ms in mss_c:
            s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self-is/NDPPP-flag.parset msin='+ms, \
                    log=ms+'_flag-c'+str(i)+'.log', cmd_type='NDPPP')
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
    s.add('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/VirA_LBA/parset_self-is/losoto.parset', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
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
        s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self-is/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument', \
              log=ms+'_selfcor-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # concat using NDPPP
    logging.info('Concat...')
    check_rm('concat.MS*')
    s.add('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parset_self-is/NDPPP-concat.parset msin="['+','.join(mss_c)+']" msout=concat.MS', \
                log='concat-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # clean (make a new model of virgo)
    logging.info('Clean...')
    if i <= 3:
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBAis/casa_clean.py', \
                params={'msfile':'concat.MS', 'imagename':'img/clean-c'+str(i), 'imtype':'cocoon', 'uvrange':uvrange}, log='clean-c'+str(i)+'.log')
    else:
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBAis/casa_clean.py', \
                params={'msfile':'concat.MS', 'imagename':'img/clean-c'+str(i), 'imtype':'normal', 'uvrange':uvrange}, log='clean-c'+str(i)+'.log')
    s.run(check=True)

#########################################################################################################
# low-res image
logging.info('Make low-resolution image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBAis/casa_clean.py', \
        params={'msfile':'concat.MS', 'imagename':'img/clean-lr', 'imtype':'lr'}, log='final_clean-lr.log')
s.run(check=True)

##########################################################################################################
# uvsub + large FoV image
logging.info('Ft+uvsub of M87 model...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', \
        params={'msfile':'concat.MS', 'model':'img/clean-c'+str(i)+'.model'}, log='final_ft.log')
s.run(check=True)
s.add('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"') # uvsub
s.run(check=False)

logging.info('Low-res wide field image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/virgoLBAis/casa_clean.py', \
        params={'msfile':'concat.MS', 'imagename':'img/clean-wide', 'imtype':'wide'}, log='final_clean-wide.log')
s.run(check=True)

logging.info("Done.")
