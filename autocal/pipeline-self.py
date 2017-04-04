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
import lofar.bdsm # import it befose casacore, bugfix
import casacore.tables as pt
from autocal.lib_pipeline import *
from make_mask import make_mask

parset_dir = '/home/fdg/scripts/autocal/parset_self/'
skymodel = '/home/fdg/scripts/model/calib-simple.skymodel'
niter = 3

if 'tooth' in os.getcwd():
    sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.LBA.skydb'
    apparent = True # no beam correction
elif 'bootes' in os.getcwd():
    sourcedb = '/home/fdg/scripts/model/Bootes_HBA.corr.skydb'
    apparent = False
else:
    # Survey
    sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/skymodels/%s_%s.skydb' % (os.getcwd().split('/')[-2], os.getcwd().split('/')[-1])
    apparent = False

#######################################################################################

set_logger('pipeline-self.logging')
check_rm('logs')
s = Scheduler(dry=False)

##################################################
# Clear
logging.info('Cleaning...')

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

# make beam
phasecentre = get_phase_centre(mss[0])
make_beam_reg(phasecentre[0], phasecentre[1], 5, 'self/beam.reg')

###############################################################################################
# Create columns (non compressed)
# TODO: remove when moving to NDPPP DFT
logging.info('Creating MODEL_DATA_HIGHRES and SUBTRACTED_DATA...')
for ms in mss:
    s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA_HIGHRES,SUBTRACTED_DATA', log=ms+'_addcol.log', cmd_type='python')
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
    if apparent:
        s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.usebeammodel=false pre.sourcedb='+ms+'/'+sourcedb_basename, log=ms+'_pre.log', cmd_type='NDPPP')
    else:
        s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.usebeammodel=true pre.sourcedb='+ms+'/'+sourcedb_basename, log=ms+'_pre.log', cmd_type='NDPPP')
s.run(check=True)

###################################################################################
# Preapre fake FR parmdb
logging.info('Prepare fake FR parmdb...')
for ms in mss:
    if os.path.exists(ms+'/instrument-fr'): continue
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
    # Re-done in case of new flags
    # TEST: higher ionfactor
    if c == 0:
        incol = 'DATA'
    else:
        incol = 'SUBTRACTED_DATA'

    logging.info('BL-based smoothing...')
    for ms in mss:
        s.add('BLsmooth.py -r -f 0.5 -i '+incol+' -o SMOOTHED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python')
    s.run(check=True)

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
    run_losoto(s, str(c), mss, [parset_dir+'/losoto-plot.parset'], ininstrument='instrument-tec', putback=False)
    os.system('mv plots-'+str(c)+' self/solutions/plots-c'+str(c))
    os.system('mv cal-'+str(c)+'.h5 self/solutions/')

    # correct TEC - group*_TC.MS:(SUBTRACTED_)DATA -> group*_TC.MS:CORRECTED_DATA
    logging.info('Correcting TEC...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' msin.datacolumn='+incol+' cor1.parmdb='+ms+'/instrument-tec cor2.parmdb='+ms+'/instrument-tec', \
                log=ms+'_cor-tec-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    #####################################################################################################
    # Cross-delay + Faraday rotation correction
    if c >= 1:

        # To circular - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (circular)
        logging.info('Convert to circular...')
        for ms in mss:
            s.add('/home/fdg/scripts/mslin2circ.py -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True)
 
        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logging.info('BL-based smoothing...')
        for ms in mss:
            s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth2-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True)

        # Solve G SB.MS:SMOOTHED_DATA (only solve)
        logging.info('Solving G...')
        for ms in mss:
            check_rm(ms+'/instrument-g')
            s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' sol.parmdb='+ms+'/instrument-g sol.solint=30 sol.nchan=8', \
                    log=ms+'_sol-g1-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        run_losoto(s, 'fr', mss, [parset_dir+'/losoto-fr.parset'], ininstrument='instrument-g', inglobaldb='globaldb',
            outinstrument='instrument-fr', outglobaldb='globaldb-fr', outtab='rotationmeasure000', putback=True)
        os.system('mv plots-fr self/solutions/')
        os.system('mv cal-fr.h5 self/solutions/')
       
        # To linear - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (linear)
        logging.info('Convert to linear...')
        for ms in mss:
            s.add('/home/fdg/scripts/mslin2circ.py -w -r -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True)
        
        # Correct FR SB.MS:CORRECTED_DATA->CORRECTED_DATA
        logging.info('Faraday rotation correction...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', log=ms+'_corFR-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # Smooth CORRECTED_DATA -> SMOOTHED_DATA
        logging.info('BL-based smoothing...')
        for ms in mss:
            s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth3-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True)

        # Solve G SB.MS:SMOOTHED_DATA (only solve)
        logging.info('Solving G...')
        for ms in mss:
            check_rm(ms+'/instrument-g')
            s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' sol.parmdb='+ms+'/instrument-g sol.solint=30 sol.nchan=8', \
                    log=ms+'_sol-g2-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        run_losoto(s, 'cd', mss, [parset_dir+'/losoto-cd.parset'], ininstrument='instrument-g', inglobaldb='globaldb',
            outinstrument='instrument-cd', outglobaldb='globaldb', outtab='amplitude000,crossdelay', putback=True)
        os.system('mv plots-cd self/solutions/')
        os.system('mv cal-cd.h5 self/solutions/')

        #run_losoto(s, 'amp', mss, [parset_dir+'/losoto-amp.parset'], ininstrument='instrument-g', inglobaldb='globaldb',
        #    outinstrument='instrument-amp', outglobaldb='globaldb', outtab='amplitude000,phase000', putback=True)
        #os.system('mv plots-amp self/solutions/')
        #os.system('mv cal-amp.h5 self/solutions/')

        # Correct FR SB.MS:CORRECTED_DATA->CORRECTED_DATA
        logging.info('Cross-delay correction...')
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-cd cor.correction=Gain', log=ms+'_corCD-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)
        # Correct slow AMP SB.MS:CORRECTED_DATA->CORRECTED_DATA
        #logging.info('Slow amp correction...')
        #for ms in mss:
        #    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-amp cor.correction=Gain', log=ms+'_corAMP-c'+str(c)+'.log', cmd_type='NDPPP')
        #s.run(check=True)
 
        ## To circular - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (circular)
        #logging.info('Convert to circular...')
        #for ms in mss:
        #    s.add('/home/fdg/scripts/mslin2circ.py -w -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin-c'+str(c)+'.log', cmd_type='python')
        #s.run(check=True)

   ###################################################################################################################
    # clen on concat.MS:CORRECTED_DATA (FR/TEC corrected, beam corrected)

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

    # clean mask clean (cut at 5k lambda)
    # no MODEL_DATA update with -baseline-averaging
    logging.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 3500 3500 -trim 2500 2500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 12arcsec -weight briggs 0.0 -auto-mask 5 -auto-threshold 1 -niter 100000 -no-update-model-required -maxuv-l 5000 -mgain 0.8 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -auto-threshold 20 '+' '.join(mss), \
            log='wscleanA-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    # make mask
    maskname = imagename+'-mask.fits'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 4)
    # remove CC not in mask
    for modelname in sorted(glob.glob(imagename+'*model.fits')):
        blank_image_fits(modelname, maskname, inverse=True)

    logging.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/wideM-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 3500 3500 -trim 2500 2500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 12arcsec -weight briggs 0.0 -auto-mask 5 -auto-threshold 1 -niter 100000 -no-update-model-required -maxuv-l 5000 -mgain 0.8 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -auto-threshold 0.1 -fitsmask '+maskname+' '+' '.join(mss), \
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

    # do low-res first cycle and remove it from the data
    if c == 0:
        # Subtract model from all TCs - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logging.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql1-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
    
        # reclean low-resolution
        logging.info('Cleaning low resolution...')
        imagename_lr = 'img/wide-lr'
        s.add('wsclean -reorder -name ' + imagename_lr + ' -size 5000 5000 -trim 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
                -scale 20arcsec -weight briggs 0.0 -auto-mask 5 -auto-threshold 1 -niter 100000 -no-update-model-required -maxuv-l 2000 -mgain 0.8 \
                -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -auto-threshold 1 '+' '.join(mss), \
                log='wsclean-lr.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
    
        # make mask
        maskname = imagename_lr+'-mask.fits'
        make_mask(image_name = imagename_lr+'-MFS-image.fits', mask_name = maskname, threshisl = 4)
        # remove CC not in mask and do not remove anything in the beam
        for modelname in glob.glob(imagename_lr+'*model.fits'):
            blank_image_fits(modelname, maskname, inverse=True)
            blank_image_reg(modelname, 'self/beam.reg', inverse=False)
    
        # resample at high res to avoid FFT problem on long baselines and predict
        logging.info('Predict...')
        for model in sorted(glob.glob(imagename_lr+'*model.fits')):
            model_out = model.replace(imagename_lr, imagename_lr+'-resamp')
            s.add('~/opt/src/nnradd/build/nnradd 10asec '+model_out+' '+model, log='resamp-lr-'+str(c)+'.log', log_append=True, cmd_type='general')
        s.run(check=True)
        s.add('wsclean -predict -name ' + imagename_lr + '-resamp -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+concat_ms, \
                log='wscleanPRE-lr.log', cmd_type='wsclean', processors='max')
        s.run(check=True)

        # corrupt model with TEC solutions ms:MODEL_DATA -> ms:MODEL_DATA
        for ms in mss:
            s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                cor1.parmdb='+ms+'/instrument-tec cor1.invert=false cor2.parmdb='+ms+'/instrument-tec cor2.invert=false', \
                log=ms+'_corrupt.log', cmd_type='NDPPP')
        s.run(check=True)
    
        # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (empty)
        logging.info('Subtracting low-res model (SUBTRACTED_DATA = DATA - MODEL_DATA)...')
        s.add('taql "update '+concat_ms+' set SUBTRACTED_DATA = DATA - MODEL_DATA"', log='taql2-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)

        # put in MODEL_DATA the best available model
        logging.info('Predict...')
        s.add('wsclean -predict -name ' + imagename + ' -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+concat_ms, \
                log='wscleanPRE-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)

    ###############################################################################################################
    # Flag on residuals (CORRECTED_DATA)
    #logging.info('Flagging residuals...')
    #for ms in mss:
    #    s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms, log=ms+'_flag-c'+str(c)+'.log', cmd_type='NDPPP')
    #s.run(check=True
    
# Subtract model from all TCs - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
logging.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql3-c'+str(c)+'.log', cmd_type='general')
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
    # if this pixel scale is changed, change also the resampling in the peel pipeline
    s.add('~/opt/src/nnradd/build/nnradd 10asec '+model_out+' '+model+' '+model_lr, log='final_resamp.log', log_append=True, cmd_type='general')
s.run(check=True) 

# Copy images
os.system('mv img/wideBeam-MFS-image.fits self/images')
os.system('mv img/wideBeamLow-MFS-image.fits self/images')
[ os.system('mv img/wide-'+str(c)+'-MFS-image.fits self/images') for c in xrange(niter) ]
[ os.system('mv img/wide-lr-MFS-image.fits self/images') for c in xrange(niter) ]
os.system('mv img/empty-image.fits self/images')
os.system('mv logs self')

logging.info("Done.")
