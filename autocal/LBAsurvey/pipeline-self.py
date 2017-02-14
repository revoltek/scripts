#!/usr/bin/python
# perform self-calibration on a group of SBs concatenated in TCs. Script must be run in dir with MS.
# number/chan in MS are flexible but the must be concatenable (same chans/freq!)
# Input:
# TCs are blocks of SBs should have calibrator corrected (a+p) data in DATA (beam not applied).
# file format of TCs is: group#_TC###.MS.
# Output:
# TCs with selfcal corrected source subtracted data in CORRECTED_DATA
# instrument tables contain gain (slow) + fast (scalarphase+TEC) solutions
# last high/low resolution models are copied in the "self/models" dir
# last high/low resolution images + masks + empty images (CORRECTED_DATA) are copied in the "self/images" dir
# h5parm solutions and plots are copied in the "self/solutions" dir

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
from autocal.lib_pipeline import *
from make_mask import make_mask

parset_dir = '/home/fdg/scripts/autocal/LBAsurvey/parset_self/'
niter = 2

if 'tooth' in os.getcwd():
    # Tooth
    skymodel = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.LBA.skymodel' # for this model remove beam in parset_self/NDPPP-predict.parset
    sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.LBA.skydb'
else:
    # Survey
    skymodel = '/home/fdg/scripts/autocal/LBAsurvey/skymodels/%s_%s.skymodel' % (os.getcwd().split('/')[-2], os.getcwd().split('/')[-1])
    sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/skymodels/%s_%s.skydb' % (os.getcwd().split('/')[-2], os.getcwd().split('/')[-1])

#######################################################################################

set_logger('pipeline-self.logging')
check_rm('logs')
s = Scheduler(dry=False)

def losoto(c, mss, parset, instrument_in='instrument', instrument_out=None):
    """
    c = cycle
    mss = list of mss
    parset = losoto parset
    instrument_in = parmdb name
    instrument_out = if None, do not copy back
    """
    logging.info('Running LoSoTo...')
    check_rm('plots')
    os.makedirs('plots')
    check_rm('globaldb')
    os.makedirs('globaldb')

    for num, ms in enumerate(mss):
        if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
        logging.debug('Copy: '+ms+'/'+instrument_in+' -> globaldb/'+instrument_in+'-'+str(num))
        os.system('cp -r '+ms+'/'+instrument_in+' globaldb/'+instrument_in+'-'+str(num))

    h5parm = 'global-c'+str(c)+'.h5'
    check_rm(h5parm)

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('losoto -v '+h5parm+' '+parset, log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    if instrument_out is not None:
        for num, ms in enumerate(mss):
            check_rm(ms+'/'+instrument_out)
            logging.debug('Move: globaldb/sol000_'+instrument_in+'-'+str(num)+' -> '+ms+'/'+instrument_out)
            os.system('mv globaldb/sol000_'+instrument_in+'-'+str(num)+' '+ms+'/'+instrument_out)
    os.system('mv plots self/solutions/plots-c'+str(c))
    os.system('mv '+h5parm+' self/solutions/')

##################################################
# Clear
logging.info('Cleaning...')

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
concat_ms = 'mss/concat.MS'

###############################################################################################
# Create columns (non compressed)
# TODO: remove when moving to NDPPP DFT
logging.info('Creating MODEL_DATA_HIGHRES...')
for ms in mss:
    s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA_HIGHRES', log=ms+'_addcol.log', cmd_type='python')
s.run(check=True)

####################################################################################################
# Add model to MODEL_DATA
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for ms in mss:
    check_rm(ms+'/'+sourcedb_basename)
    logging.debug('Copy: '+sourcedb+' -> '+ms)
    os.system('cp -r '+sourcedb+' '+ms)
logging.info('Add model to MODEL_DATA...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.sourcedb='+ms+'/'+sourcedb_basename, log=ms+'_pre.log', cmd_type='NDPPP', processors=3)
s.run(check=True)

###################################################################################
# Preapre fake FR parmdb
logging.info('Prepare fake FR parmdb...')
for ms in mss:
    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
s.run(check=True)
for ms in mss:
    s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,24)"', log=ms+'_taql.log', cmd_type='general')
s.run(check=True)

#####################################################################################################
# Self-cal cycle
for c in xrange(niter):
    logging.info('Start selfcal cycle: '+str(c))

    # Smooth DATA -> SMOOTHED_DATA
    logging.info('BL-based smoothing...')
    for ms in mss:
        s.add('BLsmooth.py -r -i DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python', processors='max') # TEST
        #s.add('BLsmooth.py -r -w -i DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python', processors='max')
    s.run(check=True, max_threads=4)

    if c == 0:
        # on first cycle concat (need to be done after smoothing)
        logging.info('Concatenating TCs...')
        check_rm(concat_ms+'*')
        pt.msutil.msconcat(mss, concat_ms, concatTime=False)

    # solve TEC - group*_TC.MS:SMOOTHED_DATA
    logging.info('Solving TEC...')
    for ms in mss:
        check_rm(ms+'/instrument-tec')
        s.add('NDPPP '+parset_dir+'/NDPPP-solTEC.parset msin='+ms+' sol.parmdb='+ms+'/instrument-tec', \
                log=ms+'_sol-tec-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # LoSoTo plot
    # TODO: add flagging?
    losoto(c, mss, parset_dir+'/losoto-plot.parset', instrument_in='instrument-tec')

    # correct TEC - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA
    # separate step, need to apply solutions to DATA, not SMOOTHED_DATA
    logging.info('Correcting TEC...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' msin.datacolumn=DATA cor1.parmdb='+ms+'/instrument-tec cor2.parmdb='+ms+'/instrument-tec', \
                log=ms+'_cor-tec-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    #####################################################################################################
    # Faraday rotation correction
    if c >= 1:

        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logging.info('BL-based smoothing...')
        for ms in mss:
            s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth2-c'+str(c)+'.log', cmd_type='python', processors='max') # TEST
            #s.add('BLsmooth.py -r -w -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth2-c'+str(c)+'.log', cmd_type='python', processors='max')
        s.run(check=True, max_threads=4)

        # Solve G SB.MS:SMOOTHED_DATA (only solve)
        logging.info('Solving G...')
        for ms in mss:
            check_rm(ms+'/instrument-g')
            s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' sol.parmdb='+ms+'/instrument-g sol.solint=10 sol.nchan=8', \
                    log=ms+'_sol-g-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # losoto
        check_rm('globaldb')
        check_rm('globaldb-fr')
        os.system('mkdir globaldb')
        os.system('mkdir globaldb-fr')
        for i, ms in enumerate(mss):
            if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
            if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb-fr/')
            num = re.findall(r'\d+', ms)[-1]
            logging.debug('Copy: '+ms+'/instrument-g -> globaldb/instrument-'+str(num))
            os.system('cp -r '+ms+'/instrument-g globaldb/instrument-'+str(num))
            logging.debug('Copy: '+ms+'/instrument-fr -> globaldb-fr/instrument-'+str(num))
            os.system('cp -r '+ms+'/instrument-fr globaldb-fr/instrument-'+str(num))
        
        logging.info('Running LoSoTo...')
        check_rm('plots')
        check_rm('global-fr.h5')
        s.add('H5parm_importer.py -v global-fr.h5 globaldb', log='losoto-fr-c'+str(c)+'.log', cmd_type='python', processors=1)
        s.run(check=True)
        s.add('losoto -v global-fr.h5 '+parset_dir+'/losoto-fr.parset', log='losoto-fr-c'+str(c)+'.log', log_append=True, cmd_type='python', processors='max')
        s.run(check=True)
        s.add('H5parm_exporter.py -v -t rotationmeasure000 global-fr.h5 globaldb-fr', log='losoto-fr-c'+str(c)+'.log', log_append=True, cmd_type='python', processors=1)
        s.run(check=True)
        os.system('mv plots self/solutions/plots-fr')
        os.system('mv global-fr.h5 self/solutions')
        
        for i, ms in enumerate(mss):
            num = re.findall(r'\d+', ms)[-1]
            check_rm(ms+'/instrument-fr')
            logging.debug('Copy globaldb-fr/sol000_instrument-'+str(num)+' -> '+ms+'/instrument-fr')
            os.system('cp -r globaldb-fr/sol000_instrument-'+str(num)+' '+ms+'/instrument-fr')

        # TODO: better stay in linear?
        # To linear - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (linear)
        logging.info('Convert to linear...')
        for ms in mss:
            s.add('/home/fdg/scripts/mslin2circ.py -s -w -r -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True, max_threads=4)
        
        # Correct FR SB.MS:CORRECTED_DATA->CORRECTED_DATA
        logging.info('Faraday rotation correction...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-corFR.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr', log=ms+'_corFR-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # To circular - SB.MS:DATA_INIT -> SB.MS:INIT_DATA (circular)
        logging.info('Convert to circular...')
        for ms in mss:
            s.add('/home/fdg/scripts/mslin2circ.py -s -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True, max_threads=4)

    ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA (FR/TEC corrected, beam corrected)

    # TODO: fix dysco
    #logging.info('Restoring WEIGHT_SPECTRUM before imaging...')
    #s.add('taql "update '+concat_ms+' set WEIGHT_SPECTRUM = WEIGHT_SPECTRUM_ORIG"', log='taql-resetweights-c'+str(c)+'.log', cmd_type='general')
    #s.run(check=True)

    # do beam-corrected+deeper image at last cycle
    if c == niter-1:
        # beam corrected: -use-differential-lofar-beam' - no baseline avg!
        logging.info('Cleaning beam (cycle: '+str(c)+')...')
        imagename = 'img/wideBeam'
        s.add('wsclean -reorder -name ' + imagename + ' -size 3500 3500 -trim 2500 2500 -mem 90 -j '+str(s.max_processors)+' \
                -scale 12arcsec -weight briggs 0.0 -auto-mask 10 -auto-threshold 1 -niter 100000 -no-update-model-required -mgain 0.8 \
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -apply-primary-beam -use-differential-lofar-beam '+' '.join(mss), \
                log='wscleanBeam-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
        # super low resolution to catch extended emission
        logging.info('Cleaning beam-low (cycle: '+str(c)+')...')
        imagename = 'img/wideBeamLow'
        s.add('wsclean -reorder -name ' + imagename + ' -size 1000 1000 -trim 512 512 -mem 90 -j '+str(s.max_processors)+' \
                -scale 1arcmin -weight natural -auto-mask 5 -auto-threshold 1 -niter 10000 -no-update-model-required -mgain 0.8 -maxuv-l 1000\
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -apply-primary-beam -use-differential-lofar-beam '+' '.join(mss), \
                log='wscleanBeamLow-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)

    # TODO: add multiscale and rms-variable-background
    # clean mask clean (cut at 8k lambda) - MODEL_DATA updated
    # -use-differential-lofar-beam -baseline-averaging
    logging.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 3500 3500 -trim 2500 2500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 12arcsec -weight briggs 0.0 -auto-mask 5 -auto-threshold 1 -niter 100000 -no-update-model-required -maxuv-l 5000 -mgain 0.8 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 '+' '.join(mss), \
            log='wscleanA-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    # make mask
    maskname = imagename+'-mask.fits'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 4)
    # remove CC not in mask
    for modelname in sorted(glob.glob(imagename+'*model.fits')):
        blank_image_fits(modelname, maskname, inverse=True)
    
    # TODO: move to DFT with NDPPP
    # update-model cannot be done in wsclean because of baseline-averaging
    logging.info('Predict...')
    s.add('wsclean -predict -name ' + imagename + ' -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+concat_ms, \
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
    s.add('wsclean -reorder -name ' + imagename + ' -size 5000 5000 -trim 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 20arcsec -weight briggs 0.0 -auto-mask 5 -auto-threshold 1 -niter 100000 -no-update-model-required -maxuv-l 2000 -mgain 0.8 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 '+' '.join(mss), \
            log='wscleanA-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    # make mask
    maskname = imagename+'-mask.fits'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 4)
    # remove CC not in mask
    for modelname in glob.glob(imagename+'*model.fits'):
        blank_image_fits(modelname, maskname, inverse=True)

    # TODO: move to DFT with NDPPP
    # resample at high res to avoid FFT problem on long baselines and predict
    logging.info('Predict...')
    for model in sorted(glob.glob(imagename+'*model.fits')):
        model_out = model.replace(imagename,imagename+'-resamp')
        s.add('~/opt/src/nnradd/build/nnradd 10asec '+model_out+' '+model, log='resamp-lr-'+str(c)+'.log', log_append=True, cmd_type='general')
    s.run(check=True)
    s.add('wsclean -predict -name ' + imagename + '-resamp -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+concat_ms, \
            log='wscleanPRE-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    ###############################################################################################################
    # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (empty)
    logging.info('Subtracting low-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
    s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql4-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    # Flag on residuals (CORRECTED_DATA)
    logging.info('Flagging residuals...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms, log=ms+'_flag-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # Concat models
    logging.info('Adding model data columns (MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA)...')
    s.add('taql "update '+concat_ms+' set MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA"', log='taql5-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

# Perform a final clean to create an inspection image of CORRECTED_DATA which should be empty
logging.info('Empty cleaning...')
imagename = 'img/empty'
s.add('wsclean -reorder -name ' + imagename + ' -size 5000 5000 -mem 90 -j '+str(s.max_processors)+' \
        -scale 5arcsec -weight briggs 0.0 -niter 1 -no-update-model-required -maxuv-l 8000 -mgain 0.6 \
        -pol I -cleanborder 0 -datacolumn CORRECTED_DATA '+concat_ms, \
        log='wsclean-empty.log', cmd_type='wsclean', processors='max')
s.run(check=True)

# Copy last *model
logging.info('Coadd+copy models...')
# resample at high res to avoid FFT problem on long baselines
for model in glob.glob('img/wide-'+str(c)+'*-model.fits'):
    if "MFS" in model: continue
    model_lr = model.replace('wide-'+str(c),'wide-lr-'+str(c)+'-resamp')
    model_out = model.replace('img/wide-'+str(c),'self/models/coadd')
    s.add('~/opt/src/nnradd/build/nnradd 10asec '+model_out+' '+model+' '+model_lr, log='final_resamp.log', log_append=True, cmd_type='general')
s.run(check=True) 

# Copy images
os.system('mv img/wideBeam-MFS-image.fits self/images')
os.system('mv img/wideBeamLow-MFS-image.fits self/images')
[ os.system('mv img/wide-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
[ os.system('mv img/wide-lr-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
os.system('mv img/empty-image.fits self/images')
os.system('mv logs self')

logging.info("Done.")
