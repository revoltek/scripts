#!/usr/bin/env python

# number of selfcal cycles
cycles = 10
# parset directory
parset_dir = '/home/fdg/scripts/autocal/AteamLBA/parset_self/'

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

##############################################
# Find right conf
localdir = os.getcwd().split('/')[-2]
if 'Cas' in localdir:
    logging.info('Observation: CasA')
    model = '/home/fdg/scripts/autocal/AteamLBA/160416_LBA-CasA.model'
    datadir = '../tgts?-bkp'
    casa_clean_parset = '/home/fdg/scripts/autocal/casa_comm/AteamLBA/casa_clean-cas.py'
    losoto_parset = 'losoto-cas.parset'

elif 'Cyg' in localdir:
    logging.info('Observation: CygA')
    model = '/home/fdg/scripts/autocal/AteamLBA/HBA-CygA.model'
    datadir = '../tgts?-bkp'
    casa_clean_parset = '/home/fdg/scripts/autocal/casa_comm/AteamLBA/casa_clean-cyg.py'
    losoto_parset = 'losoto-cyg.parset'

elif 'Tau' in localdir:
    logging.info('Observation: TauA')
    model = '/home/fdg/scripts/autocal/AteamLBA/VLA4-TauA.skydb'
    datadir = '../tgts-bkp'
    casa_clean_parset = '/home/fdg/scripts/autocal/casa_comm/AteamLBA/casa_clean-tau.py'
    losoto_parset = 'losoto-tau.parset'

elif 'Vir' in localdir:
    model = '/home/fdg/scripts/autocal/AteamLBA/150328_LBA-VirA.model'
    datadir = '../tgts-bkp'
    if 'is' in localdir:
        logging.info('Observation: VirA')
        casa_clean_parset = '/home/fdg/scripts/autocal/casa_comm/AteamLBA/casa_clean-viris.py'
        losoto_parset = 'losoto-viris.parset'
    else:
        logging.info('Observation: VirA (old)')
        casa_clean_parset = '/home/fdg/scripts/autocal/casa_comm/AteamLBA/casa_clean-vir.py'
        losoto_parset = 'losoto-vir.parset'

#################################################
# Clear
logging.info('Cleaning...')
check_rm('*last *pickle *.log')
check_rm('*h5 globaldb')
check_rm('*concat*')
check_rm('plots*')
check_rm('img')
os.makedirs('img')

###############################################
# Avg to 4 chan and 2 sec
# Remove internationals
mss = sorted(glob.glob(datadir+'/*MS'))
nchan = find_nchan(mss[0])
timeint = find_timeint(mss[0])
if nchan % 4 != 0:
    logging.error('Channels should be a multiple of 4.')
    sys.exit(1)
avg_factor_f = nchan / 4
if avg_factor_f < 1: avg_factor_f = 1
avg_factor_t = int(np.round(2/timeint))
if avg_factor_t < 1: avg_factor_t = 1
logging.info('Average in freq (factor of %i) and time (factor of %i)...' % (avg_factor_f, avg_factor_t))
for ms in mss:
    msout = ms.replace('.MS','-avg.MS').split('/')[-1]
    if os.path.exists(msout): continue
    s.add('NDPPP '+parset_dir+'/NDPPP-avg.parset msin='+ms+' msout='+msout+' msin.datacolumn=DATA avg.timestep='+str(avg_factor_t)+' avg.freqstep='+str(avg_factor_f), \
                log=msout+'_avg.log', cmd_type='NDPPP')
s.run(check=True)
nchan = nchan / avg_factor_f
timeint = timeint * avg_factor_t
mss = sorted(glob.glob('*-avg.MS'))

###############################################
# Initial processing (2/2013->2/2014)
#logging.info('Fix beam table...')
#for ms in mss:
#    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
#s.run(check=False)

##########################################################################################
# beam correction - SB.MS:DATA -> SB.MS:DATA_BEAM (beam applied, linear)
logging.info('Correcting beam...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'-init_corbeam.log', cmd_type='NDPPP')
s.run(check=True)

##########################################################################################
# Transform to circular pol - SB.MS:DATA_BEAM -> SB-circ.MS:DATA_BEAM (data, beam applied, circular)
logging.info('Convert to circular...')
for ms in mss:
    s.add('mslin2circ.py -s -i '+ms+':DATA_BEAM -o '+ms+':DATA_BEAM', log=ms+'-init_lin2circ.log', cmd_type='python')
s.run(check=True)

#########################################################################################
# Initialize columns
logging.info('Make new columns...')
for ms in mss:
    s.add('addcol2ms.py -m '+ms+' -c CORRECTED_DATA,MODEL_DATA,SMOOTHED_DATA', log=ms+'-init_addcol.log', cmd_type='python')
s.run(check=True)

# self-cal cycle
for c in xrange(cycles):
    logging.info('Starting self-cal cycle: '+str(c))

    ###########################################################################################
    # BL avg 
    # does not give good results...
    #logging.info('BL-based averaging...')
    #for ms in mss:
    #    s.add('BLavg.py -r -w -i DATA_BEAM -o SMOOTHED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python')
    #s.run(check=True)

    if c == 0:
        # After all columns are created
        logging.info('Concat...')
        # Smaller concat for ft
        for i, msg in enumerate(np.array_split(mss,10)):
            pt.msutil.msconcat(msg, 'concat-'+str(i)+'.MS', concatTime=False)
        ###########################################################################################
        # TODO: add a run of aoflagger only XY YX on combined MSs
    else:
        # first cycle use given model
        model = 'img/clean-c'+str(c-1)+'.model'

    #####################################################################################
    # ft model, model is unpolarized CIRC == LIN - SB.MS:MODEL_DATA (best m87 model)
    logging.info('Add models...')
    # if only BBS skymodel available
    #if j == 0:
    #   for j, msg in enumerate(np.array_split(mss,10)):
    #        s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin=concat-'+str(j)+'.MS pre.sourcedb='+model, log='ft-c'+str(c)+'-g'+str(j)+'.log', cmd_type='NDPPP')
    #   s.run(check=True) # not parallel!
    #else:

    for j, msg in enumerate(np.array_split(mss,10)):
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat-'+str(j)+'.MS', 'model':model}, log='ft-c'+str(c)+'-g'+str(j)+'.log')
        s.run(check=True) # NOT parallel

    #####################################################################################
    # calibrate - SB.MS:DATA_BEAM (no correction)
    logging.info('Calibrate...')
    for ms in mss:
        check_rm(ms+'/instrument')
        s.add('NDPPP '+parset_dir+'/NDPPP-selfcal_modeldata.parset msin='+ms+' msin.datacolumn=DATA_BEAM cal.parmdb='+ms+'/instrument', \
              log=ms+'_selfcal-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

#    #############################################################################################
#    # create widefield model
#    if c%3 == 0 and c != cycles-1:
#
#        logging.info('Entering wide field section:')
#
#        # apply NDPPP solutions on complete dataset - SB.MS:CIRC_DATA -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
#        # must be done before the rescaling or not all the flux is subtracted
#        logging.info('Make widefield model - Correct...')
#        for ms in mss:
#            s.add('NDPPP '+parset_dir+'/NDPPP-selfcor.parset msin='+ms+' msin.datacolumn=CIRC_DATA cor.parmdb='+ms+'/instrument', \
#                    log=ms+'_widefield-selfcor-c'+str(c)+'.log', cmd_type='NDPPP')
#        s.run(check=True)
#
#        # uvsub, MODEL_DATA is still Ateam
#        logging.info('Make widefield model - UV-Subtracting Ateam...')
#        s.add('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql-uvsub-c'+str(c)+'.log', cmd_type='general') # uvsub
#        s.run(check=False)
#
#        ###########################################################################################################################
#        # avg 15s - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, circular)
##        logging.info('Make widefield model - Average...')
##        check_rm('concat-avg.MS*')
##        s.add('NDPPP '+parset_dir+'/NDPPP-concatavg.parset msin="['+','.join(mss)+']" msout=concat-avg.MS avg.freqstep=1 avg.timestep=3', \
##                log='widefield_concatavg-c'+str(c)+'.log', cmd_type='NDPPP')
##        s.run(check=True)
#
#        # clean, mask, clean
#        logging.info('Make widefield model - Widefield imaging...')
#        imagename = 'img/clean-wide-c'+str(c)
#        s.add('wsclean -reorder -name ' + imagename + ' -size 4096 4096 -mem 90 -j '+str(s.max_processors)+' \
#                -scale 10arcsec -weight briggs 0.0 -niter 10000 -no-update-model-required -maxuv-l 5000 -mgain 0.85 -joinchannels -channelsout 20 '+' '.join(mss), \
#                log='wscleanA-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
#        s.run(check=True)
#        logging.info('Make widefield model - Make mask...')
#        make_mask(image_name = imagename+'-MFS-image.fits', mask_name = imagename+'.newmask')
#        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/AteamLBA/m87-blank.crtf'}, log='blank-c'+str(c)+'.log')
#        s.run(check=True)
#        logging.info('Make widefield model - Widefield imaging2...')
#        s.add('wsclean -reorder -name ' + imagename.replace('wide','wide-masked') + ' -size 4096 4096 -mem 90 -j '+str(s.max_processors)+' \
#                -scale 10arcsec -weight briggs 0.0 -niter 5000 -update-model-required -maxuv-l 5000 -mgain 0.85 -joinchannels -channelsout 20 -casamask '+imagename+'.newmask '+' '.join(mss), \
#                log='wscleanB-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
#        s.run(check=True)
##        widemodel = imagename.replace('wide','wide-masked')+'.model'
#
#        ###############################################################################################################################
#        # ft widefield model with wsclean
##        logging.info('Make widefield model - ft() widefield model...')
##        for ms in mss:
##            s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat-'+str(j)+'.MS', 'model':widemodel, 'wproj':512}, log='flag-ft-virgo-c'+str(c)+'-g'+str(j)+'.log')
##            s.run(check=True) # not parallel!
##            s.add('wsclean -reorder -predict ' + imagename.replace('wide','wide-masked') + ' -size 2500 2500 -mem 90 -j '+str(s.max_processors)+' -scale 10arcsec '+ms, \
##                    log=ms+'_ft-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
##        s.run(check=True)
#
#        # subtract widefield model - concat.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA=CORRECTED_DATA-MODEL_DATA (selfcal corrected data, beam applied, circular, field sources subtracted)
#        logging.info('Make widefield model - Subtract widefield model...')
#        s.add('taql "update concat.MS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql-uvsub2-c'+str(c)+'.log', cmd_type='general') # uvsub
#        s.run(check=False)
#
#        ########################################################################################################################
#        # Flagging on CORRECTED_DATA
#        logging.info('Make widefield model - Flagging residuals...')
#        for ms in mss:
#            s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms, \
#                    log=ms+'_flag-c'+str(c)+'.log', cmd_type='NDPPP')
#        s.run(check=True)
#
#        ########################################################################################################################
#        # subtract widefield model SB.MS:CIRC_DATA -> SB.MS:CIRC_DATA_SUB (uncal data with subtracted the widefield model)
#        # if last cycle skip (useless), but if one but last cycle, do on all SBs
#        for ms in mss:
#            s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' '+parset_dir+'/bbs-subcorpt.parset', \
#                  log=ms+'_subcorpt-c'+str(c)+'.log', cmd_type='BBS')
#        s.run(check=True)

    #######################################################################################
    # Solution plotting
    logging.info('Running LoSoTo to normalize solutions...')
    os.makedirs('plots')
    check_rm('globaldb')
    os.makedirs('globaldb')
    for num, ms in enumerate(mss):
        os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
        if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD globaldb/')
    h5parm = 'global-c'+str(c)+'.h5'

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', cmd_type='python')
    s.run(check=False)
    s.add('losoto -v '+h5parm+' '+parset_dir+'/'+losoto_parset, log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    for num, ms in enumerate(mss):
        check_rm(ms+'/instrument')
        os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
    os.system('mv plots plots-c'+str(c))

    ########################################################################################
    # correct - SB.MS:CIRC_DATA_SUB -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
    #logging.info('Restoring WEIGHT_SPECTRUM')
    #for ms in mss:
    #    s.add('taql "update '+ms+' set WEIGHT_SPECTRUM = WEIGHT_SPECTRUM_ORIG"', log='taql-restweights-c'+str(c)+'.log', cmd_type='general')
    #s.run(check=True)

    logging.info('Correct...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument', \
              log=ms+'_selfcor-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    ###########################################################################################################################
    # avg 1chanSB/30s - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, circular)
    # avg each timechunk separately to avoid NDPPP bug
    p = re.compile('(.*)_SB')
    timechunks = set( [ p.findall(ms)[0] for ms in mss ] )
    logging.info('Observation is divided into %i time chunks.' % len(timechunks))
    for timechunk in timechunks:
        logging.info('Average %s...' % timechunk)
        mss_tc = [ms for ms in mss if timechunk in ms]
        check_rm(timechunk+'_concat-avg.MS*')
        s.add('NDPPP '+parset_dir+'/NDPPP-concatavg.parset msin="['+','.join(mss_tc)+']" msout='+timechunk+'_concat-avg.MS avg.timestep=6 avg.freqstep=4', \
            log='concatavg-tc'+timechunk+'-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # clean (make a new model of Ateam)
    logging.info('Clean (cycle: '+str(c)+')...')
    #uvrange = '0~'+str(7+1.*c)+'klambda'
    s.add_casa(casa_clean_parset, params={'msfile':[timechunk+'_concat-avg.MS' for timechunk in timechunks], 'imagename':'img/clean-c'+str(c)}, log='clean-c'+str(c)+'.log')
    s.run(check=True)

#########################################################################################################
# low-res image
logging.info('Make low-resolution image...')
s.add_casa(casa_clean_parset, params={'msfile':'concat-avg.MS', 'imagename':'img/clean-lr', 'imtype':'lr'}, log='final_clean-lr.log')
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
s.add_casa(casa_clean_parset, params={'msfile':'concat-avg.MS', 'imagename':'img/clean-wide', 'imtype':'wide'}, log='final_clean-wide.log')
s.run(check=True)

logging.info("Done.")
