#!/usr/bin/env python
# apply solution from the calibrator and then run selfcal, an initial model
# and initial solutions from a calibrator must be provided
# local dir must contain all the MSs, un touched in linear pol
# 1 selfcal
# 2 subtract field
# 3 flag

# initial self-cal model
model = '/home/fdg/scripts/autocal/VirgoLBA/150328_LBA-VirgoA.model'
# fake skymodel with pointing direction
fakeskymodel = '/home/fdg/scripts/autocal/VirgoLBA/virgo.fakemodel.skymodel'
# number of selfcal cycles
cycles = 10
# parset directory
parset_dir = '/home/fdg/scripts/autocal/VirgoLBA/parset_self/'
# globaldb produced by pipeline-init
#globaldb = '../cals/globaldb'

##############################################################

import sys, os, glob, re
from lofar import bdsm
import numpy as np
import lsmtool
import pyrap.tables as pt
from lib_pipeline import *
from make_mask import make_mask

set_logger()
check_rm('logs')
s = Scheduler(dry=False)

#################################################
# Clear
logging.info('Cleaning...')
check_rm('*last *pickle *.log')
check_rm('*h5 globaldb')
check_rm('concat*')
check_rm('plots*')
check_rm('img')
os.makedirs('img')

# all MS
mss = sorted(glob.glob('*.MS'))

###############################################
# Initial processing
#logging.info('Fix beam table...')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)

################################################
# Copy cal solution
#logging.info('Copy solutions...')
#for ms in mss:
#    num = re.findall(r'\d+', ms)[-1]
#    logging.debug(globaldb+'/sol000_instrument-'+str(num)+' -> '+ms+'/instrument')
#    check_rm(ms+'/instrument')
#    os.system('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')

#logging.info('Averaging in time/freq...')
#for ms in mss:
#    msout = ms.replace('.MS','-avg.MS')
#    s.add('NDPPP '+parset_dir+'/NDPPP-concatavg.parset msin='+ms+' msout='+msout+' msin.datacolumn=DATA avg.freqstep=4 avg.timestep=2', log=ms+'-init_avg.log', cmd_type='NDPPP')
#s.run(check=True)

#mss = sorted(glob.glob('*-avg.MS'))

#########################################################################################
# apply solutions and beam correction - SB.MS:DATA -> SB.MS:CORRECTED_DATA (calibrator corrected data, beam applied, linear)
# TODO: convert to NDPPP - problem: does not handle DD solution - are losoto flags taken into account? - problem: if we transfer only clock is not ok
# NOTE: only beam correction, no transfer of solutions
#logging.info('Correcting target MSs...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'-init_corbeam.log', cmd_type='NDPPP')
#s.run(check=True)

#########################################################################################
# Transform to circular pol - SB.MS:CORRECTED_DATA -> SB-circ.MS:CIRC_DATA (data, beam applied, circular)
#logging.info('Convert to circular...')
#for ms in mss:
#    s.add('mslin2circ.py -i '+ms+':CORRECTED_DATA -o '+ms+':CIRC_DATA', log=ms+'-init_circ2lin.log', cmd_type='python')
#s.run(check=True)

###########################################################################################
# TODO: add a run of aoflagger only XY YX on combined MSs

########################################################################################
# Initialize columns - SB.MS:CIRC_DATA_SUB = CIRC_DATA
logging.info('Make new columns...')
for ms in mss:
    s.add('addcol2ms.py -m '+ms+' -c CORRECTED_DATA,MODEL_DATA', log=ms+'-init_addcol.log', cmd_type='python')
    #s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA,CIRC_DATA_SUB', log=ms+'-init_addcol.log', cmd_type='python')
s.run(check=True)
#logging.info('Set CIRC_DATA_SUB == CIRC_DATA...')
#for ms in mss:
#    s.add('taql "update '+ms+' set CIRC_DATA_SUB = CIRC_DATA"', log=ms+'-init_taql.log', cmd_type='general')
#s.run(check=True)

# self-cal cycle
for c in xrange(cycles):
    logging.info('Starting self-cal cycle: '+str(c))

    ###########################################################################################
    # BL avg 
    # does not give good results...
    #logging.info('BL-based averaging...')
    #for ms in mss:
    #    s.add('BLavg.py -r -w -i CIRC_DATA_SUB -o SMOOTHED_CIRC_DATA_SUB '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python')
    #s.run(check=True)

    if c == 0:
        # After all columns are created
        logging.info('Concat...')
#        pt.msutil.msconcat(mss, 'concat.MS', concatTime=False)
        # Smaller concat for ft
        for i, msg in enumerate(np.array_split(mss,10)):
            pt.msutil.msconcat(msg, 'concat-'+str(i)+'.MS', concatTime=False)
    else:
        # first cycle use given model
        model = 'img/clean-c'+str(c-1)+'.model'

    #####################################################################################
    # ft model, model is unpolarized CIRC == LIN - SB.MS:MODEL_DATA (best m87 model)
    # TODO: test if adding wprojplanes here improves the calibration, but wproj create cross problem in widefield -> must make larger image!
    logging.info('Add models...')

    for j, msg in enumerate(np.array_split(mss,10)):
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat-'+str(j)+'.MS', 'model':model}, log='ft-virgo-c'+str(c)+'-g'+str(j)+'.log')
        s.run(check=True) # not parallel!

    #####################################################################################
    # calibrate - SB.MS:DATA (no correction)
    logging.info('Calibrate...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-selfcal_modeldata.parset msin='+ms+' cal.parmdb='+ms+'/instrument', \
              log=ms+'_selfcal-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    #############################################################################################
    # create widefield model
    if c%3 == 0 and c != cycles-1 and False: # TODO skip

        logging.info('Entering wide field section:')

        # apply NDPPP solutions on complete dataset - SB.MS:CIRC_DATA -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
        # must be done before the rescaling or not all the flux is subtracted
        logging.info('Make widefield model - Correct...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-selfcor.parset msin='+ms+' msin.datacolumn=CIRC_DATA cor.parmdb='+ms+'/instrument', \
                    log=ms+'_widefield-selfcor-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # uvsub, MODEL_DATA is still Virgo
        logging.info('Make widefield model - UV-Subtracting Virgo A...')
        s.add('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql-uvsub-c'+str(c)+'.log', cmd_type='general') # uvsub
        s.run(check=False)

        ###########################################################################################################################
        # avg 15s - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, circular)
#        logging.info('Make widefield model - Average...')
#        check_rm('concat-avg.MS*')
#        s.add('NDPPP '+parset_dir+'/NDPPP-concatavg.parset msin="['+','.join(mss)+']" msout=concat-avg.MS avg.freqstep=1 avg.timestep=3', \
#                log='widefield_concatavg-c'+str(c)+'.log', cmd_type='NDPPP')
#        s.run(check=True)

        # clean, mask, clean
        logging.info('Make widefield model - Widefield imaging...')
        imagename = 'img/clean-wide-c'+str(c)
#        s.add_casa('/home/fdg/scripts/autocal/casa_comm/VirgoLBA/casa_clean.py', \
#                    params={'msfile':'concat-avg.MS', 'imagename':imagename, 'imtype':'wide'}, log='clean-wide1-c'+str(c)+'.log')
        s.add('wsclean -reorder -name ' + imagename + ' -size 4096 4096 -mem 90 -j '+str(s.max_processors)+' \
                -scale 10arcsec -weight briggs 0.0 -niter 10000 -no-update-model-required -maxuv-l 5000 -mgain 0.85 -joinchannels -channelsout 20 '+' '.join(mss), \
                log='wscleanA-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
        logging.info('Make widefield model - Make mask...')
        make_mask(image_name = imagename+'-MFS-image.fits', mask_name = imagename+'.newmask')
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/VirgoLBA/m87-blank.crtf'}, log='blank-c'+str(c)+'.log')
        s.run(check=True)
        logging.info('Make widefield model - Widefield imaging2...')
#        s.add_casa('/home/fdg/scripts/autocal/casa_comm/VirgoLBA/casa_clean.py', \
#                    params={'msfile':'concat-avg.MS', 'imagename':imagename.replace('wide','wide-masked'), 'mask':imagename+'.newmask', 'imtype':'widemasked'}, log='clean-wide2-c'+str(c)+'.log')
        s.add('wsclean -reorder -name ' + imagename.replace('wide','wide-masked') + ' -size 4096 4096 -mem 90 -j '+str(s.max_processors)+' \
                -scale 10arcsec -weight briggs 0.0 -niter 5000 -update-model-required -maxuv-l 5000 -mgain 0.85 -joinchannels -channelsout 20 -casamask '+imagename+'.newmask '+' '.join(mss), \
                log='wscleanB-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
#        widemodel = imagename.replace('wide','wide-masked')+'.model'

        ###############################################################################################################################
        # ft widefield model with wsclean
#        logging.info('Make widefield model - ft() widefield model...')
#        for ms in mss:
#            s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat-'+str(j)+'.MS', 'model':widemodel, 'wproj':512}, log='flag-ft-virgo-c'+str(c)+'-g'+str(j)+'.log')
#            s.run(check=True) # not parallel!
#            s.add('wsclean -reorder -predict ' + imagename.replace('wide','wide-masked') + ' -size 2500 2500 -mem 90 -j '+str(s.max_processors)+' -scale 10arcsec '+ms, \
#                    log=ms+'_ft-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
#        s.run(check=True)

        # subtract widefield model - concat.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA=CORRECTED_DATA-MODEL_DATA (selfcal corrected data, beam applied, circular, field sources subtracted)
        logging.info('Make widefield model - Subtract widefield model...')
        s.add('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql-uvsub2-c'+str(c)+'.log', cmd_type='general') # uvsub
        s.run(check=False)

        ########################################################################################################################
        # Flagging on CORRECTED_DATA
        logging.info('Make widefield model - Flagging residuals...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms, \
                    log=ms+'_flag-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        ########################################################################################################################
        # subtract widefield model SB.MS:CIRC_DATA -> SB.MS:CIRC_DATA_SUB (uncal data with subtracted the widefield model)
        # if last cycle skip (useless), but if one but last cycle, do on all SBs
        for ms in mss:
            s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' '+parset_dir+'/bbs-subcorpt.parset', \
                  log=ms+'_subcorpt-c'+str(c)+'.log', cmd_type='BBS')
        s.run(check=True)

    #######################################################################################
    # Solution plotting
    logging.info('Running LoSoTo to normalize solutions...')
    os.makedirs('plots')
    check_rm('globaldb')
    os.makedirs('globaldb')
    for num, ms in enumerate(mss):
        os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
        if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
    h5parm = 'global-c'+str(c)+'.h5'

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', cmd_type='python')
    s.run(check=False)
    s.add('losoto -v '+h5parm+' '+parset_dir+'/losoto.parset', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    for num, ms in enumerate(mss):
        check_rm(ms+'/instrument')
        os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
    os.system('mv plots plots-c'+str(c))

    ########################################################################################
    # correct - SB.MS:CIRC_DATA_SUB -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
#    logging.info('Restoring WEIGHT_SPECTRUM')
#    s.add('taql "update concat.MS set WEIGHT_SPECTRUM = WEIGHT_SPECTRUM_ORIG"', log='taql-restweights-c'+str(c)+'.log', cmd_type='general')
#    s.run(check=True)

    logging.info('Correct...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument', \
              log=ms+'_selfcor-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    ###########################################################################################################################
    # avg 1chanSB/30s - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, circular)
    logging.info('Average...')
    check_rm('concat-avg.MS*')
    s.add('NDPPP '+parset_dir+'/NDPPP-concatavg.parset msin="['+','.join(mss)+']" msout=concat-avg.MS avg.timestep=6', \
            log='concatavg-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # clean (make a new model of virgo)
    logging.info('Clean (cycle: '+str(c)+')...')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/VirgoLBA/casa_clean.py', \
            params={'msfile':'concat-avg.MS', 'imagename':'img/clean-c'+str(c)}, log='clean-c'+str(c)+'.log')
    s.run(check=True)

#########################################################################################################
# low-res image
logging.info('Make low-resolution image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/VirgoLBA/casa_clean.py', \
        params={'msfile':'concat-avg.MS', 'imagename':'img/clean-lr', 'imtype':'lr'}, log='final_clean-lr.log')
s.run(check=True)

##########################################################################################################
# uvsub + large FoV image
logging.info('Ft+uvsub of M87 model...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', \
        params={'msfile':'concat-avg.MS', 'model':'img/clean-c'+str(c)+'.model'}, log='final_ft.log')
s.run(check=True)
s.add('taql "update concat-avg.MS set DATA = DATA - MODEL_DATA"') # uvsub
s.run(check=False)

logging.info('Low-res wide field image...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/VirgoLBA/casa_clean.py', \
        params={'msfile':'concat-avg.MS', 'imagename':'img/clean-wide', 'imtype':'wide'}, log='final_clean-wide.log')
s.run(check=True)

logging.info("Done.")
