#!/usr/bin/python
# perform self-calibration on a group of SBs concatenated in TCs. Script must be run in dir with MS.
# number/chan in MS are flexible but the must be concatenable (same chans/freq!)
# Input:
# TCs are blocks of SBs should have calibrator corrected (a+p) data in DATA (beam not applied).
# file format of TCs is: group#_TC###.MS.
# Output:
# TCs with selfcal NOT corrected source subtracted data in SUBTRACTED_DATA
# TCs with selfcal corrected source subtracted data in CORRECTED_DATA
# instrument tables contain gain (slow) + fast (scalarphase+TEC) solutions
# last high/low resolution models are copied in the "self/models" dir
# last high/low resolution images + masks + empty images (CORRECTED_DATA) are copied in the "self/images" dir
# h5parm solutions and plots are copied in the "self/solutions" dir

import sys, os, glob, re
import numpy as np
from lofar import bdsm
import pyrap.tables as pt
from lib_pipeline import *
from make_mask import make_mask

parset_dir = '/home/fdg/scripts/autocal/LBAsurvey/parset_self/'

# Tooth
skymodel = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.LBA.skymodel' # for this model remove beam in parset_self/NDPPP-predict.parset
sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.LBA.skydb'

# Survey
#skymodel = '/home/stsf309/scripts/autocal/LBAsurvey/skymodels/%s_%s.skymodel' % (os.getcwd().split('/')[-2], os.getcwd().split('/')[-1])
#sourcedb = '/home/stsf309/scripts/autocal/LBAsurvey/skymodels/%s_%s.skydb' % (os.getcwd().split('/')[-2], os.getcwd().split('/')[-1])

niter = 2

#######################################################################################

set_logger()
check_rm('logs')
s = Scheduler(dry=False)

def losoto(c, mss, parset, instrument_in='instrument', instrument_out='instrument'):
    """
    c = cycle
    mss = list of mss
    parset = losoto parset
    """
    logging.info('Running LoSoTo...')
    check_rm('plots')
    os.makedirs('plots')
    check_rm('globaldb')
    os.makedirs('globaldb')

    for num, ms in enumerate(mss):
        if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
        os.system('cp -r '+ms+'/'+instrument_in+' globaldb/'+instrument_in+'-'+str(num))

    h5parm = 'global-c'+str(c)+'.h5'
    check_rm(h5parm)

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('losoto -v '+h5parm+' '+parset, log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    for num, ms in enumerate(mss):
        check_rm(ms+'/'+instrument_out)
        os.system('mv globaldb/sol000_'+instrument_in+'-'+str(num)+' '+ms+'/'+instrument_out)
    os.system('mv plots self/solutions/plots-c'+str(c))
    os.system('mv '+h5parm+' self/solutions/')

##################################################
## Clear
#logging.info('Cleaning...')
#
check_rm('*bak *last *pickle')
check_rm('img')
os.makedirs('img')
os.makedirs('logs/mss')

# here images, models, solutions for each group will be saved
check_rm('self')
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/models'): os.makedirs('self/models')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')

mss = sorted(glob.glob('mss/TC*[0-9].MS'))
nchan = find_nchan(mss[0])
concat_ms = 'mss/concat.MS'

###################################################################################################
# Add model to MODEL_DATA
logging.info('Add model to MODEL_DATA...')
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for ms in mss:
    check_rm(ms+'/'+sourcedb_basename)
    os.system('cp -r '+sourcedb+' '+ms)
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.sourcedb='+ms+'/'+sourcedb_basename, log=ms+'_pre.log', cmd_type='NDPPP', processors=3)
s.run(check=True)

## 1. find and remove FR

################################################################################################
# Smooth DATA -> SMOOTHED_DATA (circular, smooth)
logging.info('BL-based smoothing...')
for ms in mss:
    s.add('BLavg.py -r -w -i DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python', processors='max')
s.run(check=True, max_threads=2)

#################################################################################################
# solve+correct TEC - group*_TC.MS:SMOOTHED_DATA -> group*_TC.MS:CORRECTED_DATA (circular, smooth, TEC-calibrated)
# TODO: merge in a single step with new NDPPP
# TODO: calibrate also fast CSA?
logging.info('Calibrating TEC...')
for ms in mss:
    check_rm(ms+'/instrument-tecinit')
    s.add('NDPPP '+parset_dir+'/NDPPP-solTEC.parset msin='+ms+' cal.parmdb='+ms+'/instrument-tecinit', log=ms+'_sol-tecinit.log', cmd_type='NDPPP')
s.run(check=True)
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' msin.datacolumn=SMOOTHED_DATA \
            cor1.parmdb='+ms+'/instrument-tecinit cor2.parmdb='+ms+'/instrument-tecinit', log=ms+'_cor-tecinit.log', cmd_type='NDPPP')
s.run(check=True)

##############################################################################################
# Solve G SB.MS:CORRECTED_DATA (only solve)
# NOTE: test with solint=10 and nchan=8
logging.info('Calibrating G...')
for ms in mss:
    check_rm(ms+'/instrument-ginit')
    s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA cal.parmdb='+ms+'/instrument-ginit cal.solint=10 cal.nchan=8', log=ms+'_sol-g.log', cmd_type='NDPPP')
s.run(check=True)

##################################################################################
# Preapre fake FR parmdb
logging.info('Prepare fake FR parmdb...')
for ms in mss:
    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
s.run(check=True)
for ms in mss:
    s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,24)"', log=ms+'_taql.log', cmd_type='general')
s.run(check=True)

# merge parmdbs
logging.info('Merging instrument tables...')
for ms in mss:
    merge_parmdb(ms+'/instrument-tecinit', ms+'/instrument-ginit', ms+'/instrument', clobber=True)

#################################################
# Prepare and run losoto
check_rm('globaldb')
check_rm('globaldb-fr')
os.system('mkdir globaldb')
os.system('mkdir globaldb-fr')
for i, ms in enumerate(mss):
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb-fr/')
    num = re.findall(r'\d+', ms)[-1]
    logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
    logging.debug('Copy instrument-fr of '+ms+' into globaldb-fr/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument-fr globaldb-fr/instrument-fr-'+str(num))

logging.info('Running LoSoTo...')
check_rm('plots')
check_rm('global-fr.h5')
s.add('H5parm_importer.py -v global-fr.h5 globaldb', log='losoto1.log', cmd_type='python', processors=1)
s.run(check=True)
s.add('losoto -v global-fr.h5 '+parset_dir+'/losoto-plot.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('losoto -v global-fr.h5 '+parset_dir+'/losoto-fr.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('H5parm_exporter.py -v -t rotationmeasure000 global-fr.h5 globaldb-fr', log='losoto1.log', log_append=True, cmd_type='python', processors=1)
s.run(check=True)
os.system('mv plots self/solutions/plots-fr')
os.system('mv global-fr.h5 self/solutions')

for i, ms in enumerate(mss):
    num = re.findall(r'\d+', ms)[-1]
    check_rm(ms+'/instrument-fr')
    logging.debug('Copy globaldb-fr/sol000_instrument-fr-'+str(num)+' into '+ms+'/instrument-fr')
    os.system('cp -r globaldb-fr/sol000_instrument-fr-'+str(num)+' '+ms+'/instrument-fr')

###################################################################################################
# To linear - SB.MS:DATA -> SB.MS:CORRECTED_DATA (linear)
# TODO: better stay in circular?
logging.info('Convert to linear...')
for ms in mss:
    s.add('/home/fdg/scripts/mslin2circ.py -s -w -r -i '+ms+':DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
s.run(check=True, max_threads=1)

####################################################
# Correct FR SB.MS:CORRECTED_DATA->DATA_INIT
logging.info('Faraday rotation correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-corFR.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr', log=ms+'_corFR.log', cmd_type='NDPPP')
s.run(check=True)

###################################################################################################
# To circular - SB.MS:DATA_INIT -> SB.MS:INIT_DATA (circular)
logging.info('Convert to circular...')
for ms in mss:
    s.add('/home/fdg/scripts/mslin2circ.py -s -w -i '+ms+':DATA_INIT -o '+ms+':DATA_INIT', log=ms+'_circ2lin.log', cmd_type='python')
s.run(check=True, max_threads=1)

# 2: recalibrate without FR

###############################################################################################
# Create columns
logging.info('Creating MODEL_DATA_HIGHRES, SUBTRACTED_DATA...')
for ms in mss:
    s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA_HIGHRES,SUBTRACTED_DATA', log=ms+'_addcol.log', cmd_type='python')
s.run(check=True, max_threads=2)

####################################################################################################
## Self-cal cycle
for c in xrange(niter):
    logging.info('Start selfcal cycle: '+str(c))

    #################################################################################################
    # Smooth DATA_INIT -> SMOOTHED_DATA
    logging.info('BL-based smoothing...')
    for ms in mss:
        s.add('BLavg.py -r -w -i DATA_INIT -o SMOOTHED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python', processors='max')
    s.run(check=True, max_threads=2)

    if c == 0:
        # on first cycle concat (need to be done after smoothing)
        logging.info('Concatenating TCs...')
        check_rm(concat_ms+'*')
        pt.msutil.msconcat(mss, concat_ms, concatTime=False)

    # solve TEC - group*_TC.MS:SMOOTHED_DATA
    logging.info('Solving TEC...')
    for ms in mss:
        check_rm(ms+'/instrument-tec')
        s.add('NDPPP '+parset_dir+'/NDPPP-solTEC.parset msin='+ms+' cal.parmdb='+ms+'/instrument-tec', log=ms+'_sol-tec-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # LoSoTo plot
    # TODO: add flagging
    losoto(c, mss, parset_dir+'/losoto-plot.parset', instrument_in='instrument-tec', instrument_out='instrument-tec')

    # correct TEC - group*_TC.MS:DATA_INIT -> group*_TC.MS:CORRECTED_DATA
    logging.info('Correcting TEC...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' msin.datacolumn=DATA_INIT cor1.parmdb='+ms+'/instrument-tec cor2.parmdb='+ms+'/instrument-tec', \
                log=ms+'_cor-tec-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    ################################################################################################
    ## Solve SB.MS:CORRECTED_DATA (only solve)
    ## TESTTESTTEST
    #logging.info('Calibrating G...')
    #for ms in mss:
    #    check_rm(ms+'/instrument-g')
    #    s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA cal.parmdb='+ms+'/instrument-g cal.solint=30 cal.nchan=4', log=ms+'_sol-g-c'+str(c)+'.log', cmd_type='NDPPP')
    #s.run(check=True)
    #losoto(c, mss, parset_dir+'/losoto-fr.parset', instrument_in='instrument-g', instrument_out='instrument-g')

    logging.info('Restoring WEIGHT_SPECTRUM before imaging...')
    s.add('taql "update '+concat_ms+' set WEIGHT_SPECTRUM = WEIGHT_SPECTRUM_ORIG"', log='taql-resetweights-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA (FR/TEC corrected, beam corrected)

    # clean mask clean (cut at 8k lambda) - MODEL_DATA updated
    # -use-differential-lofar-beam -baseline-averaging
    # TEST: go to 5k from 8k and to 10arcsec from 5 arcsec and from 5000 to 2500 in size
    logging.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 3000 3000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 10arcsec -weight briggs 0.0 -auto-threshold 5 -niter 8000 -no-update-model-required -maxuv-l 5000 -mgain 0.75 \
            -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 '+concat_ms, \
            log='wscleanA-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    maskname = imagename+'.newmask'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
               params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/LBAsurvey/tooth_mask.crtf', 'setTo':1}, log='casablank-c'+str(c)+'.log')
    s.run(check=True)
    # TODO: automasking with wsclean 2.1
    logging.info('Cleaning with mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 3000 3000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 10arcsec -weight briggs 0.0 -auto-threshold 5 -niter 5000 -no-update-model-required -maxuv-l 5000 -mgain 0.75 \
            -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -casamask '+maskname+' '+concat_ms, \
            log='wscleanB-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    
    logging.info('Predict...')
    s.add('wsclean -predict -name ' + imagename + ' -size 2500 2500 -mem 90 -j '+str(s.max_processors)+' \
            -scale 10arcsec -channelsout 10 '+concat_ms, \
            log='wscleanPRE-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    ####################################################################
    # FAST VERSION (no low-res)
    #continue
    ####################################################################

    logging.info('Moving MODEL_DATA to MODEL_DATA_HIGHRES...')
    s.add('taql "update '+concat_ms+' set MODEL_DATA_HIGHRES = MODEL_DATA"', log='taql1-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    ############################################################################################################
    # Subtract model from all TCs - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
    logging.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
    s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql3-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    # reclean low-resolution
    logging.info('Cleaning low resolution (cycle: '+str(c)+')...')
    imagename = 'img/wide-lr-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 20arcsec -weight briggs 0.0 -auto-threshold 5 -niter 5000 -no-update-model-required -maxuv-l 2000 -mgain 0.75 \
            -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 '+concat_ms, \
            log='wscleanA-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    # TODO: remove re-imaging and just keep CC into mask
    maskname = imagename+'.newmask'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshpix=6) # a bit higher treshold
    logging.info('Cleaning low resolution with mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-lr-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 20arcsec -weight briggs 0.0 -auto-threshold 5 -niter 3000 -no-update-model-required -maxuv-l 2000 -mgain 0.75 \
            -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -casamask '+maskname+' '+concat_ms, \
            log='wscleanB-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    # resample at high res to avoid FFT problem on long baselines
    logging.info('Predict...')
    for model in glob.glob(imagename+'*model.fits'):
        model_out = model.replace(imagename,imagename+'-resamp')
        s.add('~/opt/src/nnradd/build/nnradd 10asec '+model_out+' '+model, log='resamp-lr-'+str(c)+'.log', log_append=True, cmd_type='general')
    s.run(check=True)
    s.add('wsclean -predict -name ' + imagename + '-resamp -size 8000 8000 -mem 90 -j '+str(s.max_processors)+' \
            -scale 10arcsec -channelsout 10 '+concat_ms, \
            log='wscleanPRE-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    ###############################################################################################################
    # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (empty)
    logging.info('Subtracting low-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
    s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql4-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    # Flag on residuals
    # TODO: add Sarod code
    logging.info('Flagging residuals...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms, \
                log=ms+'_flag-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # Concat models
    logging.info('Adding model data columns (MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA)...')
    s.add('taql "update '+concat_ms+' set MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA"', log='taql5-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

# Perform a final clean to create an inspection image of SUBTRACTED_DATA which should be very empty
logging.info('Empty cleaning...')
s.add('taql "update '+concat_ms+' set SUBTRACTED_DATA = CORRECTED_DATA"', log='taql_sub.log', cmd_type='general')
s.run(check=True)
imagename = 'img/empty'
s.add('wsclean -reorder -name ' + imagename + ' -size 5000 5000 -mem 90 -j '+str(s.max_processors)+' \
        -scale 5arcsec -weight briggs 0.0 -niter 1 -no-update-model-required -maxuv-l 8000 -mgain 0.6 \
        -pol I -cleanborder 0 -datacolumn SUBTRACTED_DATA '+concat_ms, \
        log='wsclean-empty.log', cmd_type='wsclean', processors='max')
s.run(check=True)

# Copy last *model
logging.info('Coadd+copy models...')
# resample at high res to avoid FFT problem on long baselines
for model in glob.glob('img/wideM-'+str(c)+'*-model.fits'):
    if "MFS" in model: continue
    model_lr = model.replace('wideM-'+str(c),'wideM-lr-'+str(c)+'-resamp')
    model_out = model.replace('img/wideM-'+str(c),'self/models/coadd')
    s.add('~/opt/src/nnradd/build/nnradd 10asec '+model_out+' '+model+' '+model_lr, log='final_resamp.log', log_append=True, cmd_type='general')
s.run(check=True) 

# Copy images
[ os.system('mv img/wide-'+str(c)+'.newmask self/images') for c in xrange(niter) ]
[ os.system('mv img/wide-lr-'+str(c)+'.newmask self/images') for c in xrange(niter) ]
[ os.system('mv img/wide-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
[ os.system('mv img/wide-lr-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
[ os.system('mv img/wideM-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
[ os.system('mv img/wideM-lr-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
os.system('mv img/empty-image.fits self/images')
os.system('mv logs self')

logging.info("Done.")
