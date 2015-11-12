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

parset_dir = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/'
skymodel = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/toothbrush.GMRT150.skymodel'
niter = 2

#######################################################################################

import sys, os, glob, re
import numpy as np
from lofar import bdsm
import pyrap.tables as pt
from lib_pipeline import *
from make_mask import make_mask

set_logger()
s = Scheduler(dry=False)

def losoto(c, mss, g, parset):
    """
    c = cycle
    mss = list of mss
    g = group number
    parset = losoto parset
    """
    logging.info('Running LoSoTo...')
    check_rm('plots')
    os.makedirs('plots')
    check_rm('globaldb')
    os.makedirs('globaldb')

    for num, ms in enumerate(mss):
        os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
        if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
    h5parm = 'global-c'+str(c)+'.h5'

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('losoto -v '+h5parm+' '+parset, log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    for num, ms in enumerate(mss):
        check_rm(ms+'/instrument')
        os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
    os.system('mv plots self/solutions/g'+g+'/plots-c'+str(c))
    os.system('mv '+h5parm+' self/solutions/g'+g)


# here images, models, solutions for each group will be saved
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/models'): os.makedirs('self/models')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')

for group in sorted(glob.glob('group*'))[::-1]:

    mss_orig = sorted(glob.glob(group+'/group*_TC*[0-9].MS'))
    concat_ms = group+'/concat.MS'
    concat_ms_orig = group+'/concat_orig.MS'
    g = str(re.findall(r'\d+', mss_orig[0])[0])
    logging.info('Working on group: '+g+'...')
    
    ################################################################################################
    # Clear
    logging.info('Cleaning...')
    check_rm(group+'/*log *log *bak *last *pickle')
    check_rm(group+'/plots* plots')
    check_rm(group+'/*h5 *h5 globaldb')
    check_rm('*last')
    check_rm('img')
    os.makedirs('img')
    check_rm('self/images/g'+g)
    os.makedirs('self/images/g'+g)
    check_rm('self/solutions/g'+g)
    os.makedirs('self/solutions/g'+g)

    #################################################################################################
    # Create columns
    logging.info('Creating MODEL_DATA_HIGHRES, SUBTRACTED_DATA, MODEL_DATA and CORRECTED_DATA...')
    for ms in mss_orig:
        s.add('addcol2ms.py -i '+ms+' -o MODEL_DATA,CORRECTED_DATA,MODEL_DATA_HIGHRES,SUBTRACTED_DATA', log=ms+'_addcol.log', cmd_type='python')
    s.run(check=True)

    ###################################################################################################
    # Separate LL and RR
    logging.info('Separating RR and LL...')
    for ms in mss_orig:
        msLL = ms.replace('.MS','-LL.MS')
        if os.path.exists(msLL): os.system('rm -r '+msLL)
        os.system( 'cp -r '+ms+' '+msLL )

        msRR = ms.replace('.MS','-RR.MS')
        if os.path.exists(msRR): os.system('rm -r '+msRR)
        os.system( 'cp -r '+ms+' '+msRR )

    for ms in mss_orig:
        msLL = ms.replace('.MS','-LL.MS')
        msRR = ms.replace('.MS','-RR.MS')
        s.add('taql "update '+msRR+' set DATA[,3]=DATA[,0]"', log=ms+'_init-taql.log', cmd_type='general', log_append=True)
        s.add('taql "update '+msLL+' set DATA[,0]=DATA[,3]"', log=ms+'_init-taql.log', cmd_type='general', log_append=True)
    s.run(check=True)
    for ms in mss_orig:
        msLL = ms.replace('.MS','-LL.MS')
        msRR = ms.replace('.MS','-RR.MS')
        s.add('taql "update '+msRR+' set WEIGHT_SPECTRUM[,3]=WEIGHT_SPECTRUM[,0]"', log=ms+'_init-taql2.log', cmd_type='general', log_append=True)
        s.add('taql "update '+msLL+' set WEIGHT_SPECTRUM[,0]=WEIGHT_SPECTRUM[,3]"', log=ms+'_init-taql2.log', cmd_type='general', log_append=True)
    s.run(check=True)
    for ms in mss_orig:
        msLL = ms.replace('.MS','-LL.MS')
        msRR = ms.replace('.MS','-RR.MS')
        s.add('taql "update '+msRR+' set FLAG[,3]=FLAG[,0]"', log=ms+'_init-taql3.log', cmd_type='general', log_append=True)
        s.add('taql "update '+msLL+' set FLAG[,0]=FLAG[,3]"', log=ms+'_init-taql3.log', cmd_type='general', log_append=True)
    s.run(check=True)

    mss = sorted(glob.glob(group+'/group*_TC*[0-9]-*.MS'))
    mssrr = sorted(glob.glob(group+'/group*_TC*[0-9]-RR.MS'))
    mssll = sorted(glob.glob(group+'/group*_TC*[0-9]-LL.MS'))

    ####################################################################################################
    # Self-cal cycle
    for c in xrange(niter):
        logging.info('Start selfcal cycle: '+str(c))
        
        if c != 0:
            # copy FLAG and MODEL_DATA in RR/LL - group*_TC-*-[RR|LL].MS:MODEL_DATA = group*_TC.MS:MODEL_DATA
            logging.info('Copying FLAG and MODEL_DATA to the averaged datasets...')
            for ms in mss_orig:
                msrr = ms.replace('.MS','-RR.MS')
                s.add('taql "update '+msrr+', '+ms+' as orig set MODEL_DATA[,0]=orig.MODEL_DATA[,0]" && \
                       taql "update '+msrr+', '+ms+' as orig set MODEL_DATA[,3]=orig.MODEL_DATA[,0]" && \
                       taql "update '+msrr+', '+ms+' as orig set FLAG[,0]=orig.FLAG[,0]" && \
                       taql "update '+msrr+', '+ms+' as orig set FLAG[,3]=orig.FLAG[,0]"', \
                       log=ms+'taql_copymodel-c'+str(c)+'.log', cmd_type='general')
            s.run(check=True)
            for ms in mss_orig:
                msll = ms.replace('.MS','-LL.MS')
                s.add('taql "update '+msll+', '+ms+' as orig set MODEL_DATA[,0]=orig.MODEL_DATA[,3]" && \
                       taql "update '+msll+', '+ms+' as orig set MODEL_DATA[,3]=orig.MODEL_DATA[,3]" && \
                       taql "update '+msll+', '+ms+' as orig set FLAG[,0]=orig.FLAG[,3]" && \
                       taql "update '+msll+', '+ms+' as orig set FLAG[,3]=orig.FLAG[,3]"', \
                       log=ms+'taql_copymodel-c'+str(c)+'.log', cmd_type='general')
            s.run(check=True)

        #################################################################################################
        # Smooth
        logging.info('Smoothing...')
        for ms in mss:
            s.add('BLavg.py -r -w -i DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True)

        if c == 0:

            # on first cycle concat (need to be done after smoothing)
            logging.info('Concatenating TCs...')
            check_rm(concat_ms+'*')
            pt.msutil.msconcat(mss, concat_ms, concatTime=False)

            # calibrate phase-only (only solve) - group*_TC.MS:SMOOTHED_DATA (beam: ARRAY_FACTOR)
            logging.info('Calibrating phase...')
            for ms in mss:
                s.add('calibrate-stand-alone -f '+ms+' '+parset_dir+'/bbs-sol_tec.parset '+skymodel, \
                      log=ms+'_soltec-c'+str(c)+'.log', cmd_type='BBS')
            s.run(check=True)

            # plots
            losoto(str(c)+'rr', mssrr, g, parset_dir+'/losoto-plot.parset')
            losoto(str(c)+'ll', mssll, g, parset_dir+'/losoto-plot.parset')

            # correct phase-only - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase corrected, beam corrected)
            logging.info('Correcting phase...')
            for ms in mss:
                s.add('calibrate-stand-alone '+ms+' '+parset_dir+'/bbs-cor_tec.parset '+skymodel, \
                log=ms+'_cortec-c'+str(c)+'.log', cmd_type='BBS')
            s.run(check=True)

        else:

            # calibrate phase-only (only solve) - group*_TC.MS:SMOOTHED_DATA @ MODEL_DATA
            logging.info('Calibrating phase...')
            for ms in mss:
                s.add('calibrate-stand-alone -f --parmdb-name instrument_tec '+ms+' '+parset_dir+'/bbs-sol_tec-preamp.parset '+skymodel, \
                      log=ms+'_solpreamp-c'+str(c)+'.log', cmd_type='BBS')
            s.run(check=True)

            # calibrate phase-only - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase corrected, beam corrected)
            logging.info('Correcting phase...')
            for ms in mss:
                s.add('calibrate-stand-alone --parmdb-name instrument_tec '+ms+' '+parset_dir+'/bbs-cor_tec-preamp.parset '+skymodel, \
                log=ms+'_corpreamp-c'+str(c)+'.log', cmd_type='BBS')
            s.run(check=True)

            # Smooth
            logging.info('Smoothing...')
            for ms in mss:
                s.add('BLavg.py -r -w -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth-preamp-c'+str(c)+'.log', cmd_type='python')
            s.run(check=True)

            # calibrate amplitude (only solve) - group*_TC.MS:SMOOTHED_DATA @ MODEL_DATA
            logging.info('Calibrating amplitude...')
            for ms in mss:
                s.add('calibrate-stand-alone -f --parmdb-name instrument_amp '+ms+' '+parset_dir+'/bbs-sol_amp.parset '+skymodel, \
                      log=ms+'_calamp-c'+str(c)+'.log', cmd_type='BBS')
            s.run(check=True)

            # merge parmdbs
            logging.info('Merging instrument tables...')
            for ms in mss:
                merge_parmdb(ms+'/instrument_tec', ms+'/instrument_amp', ms+'/instrument', clobber=True)
    
            ########################################################
            # LoSoTo Amp rescaling
            losoto(str(c)+'rr', mssrr, g, parset_dir+'/losoto.parset')
            losoto(str(c)+'ll', mssll, g, parset_dir+'/losoto.parset')
        
            # correct - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase+amp corrected, beam corrected)
            logging.info('Correcting...')
            for ms in mss:
                s.add('calibrate-stand-alone '+ms+' '+parset_dir+'/bbs-cor_amptec.parset '+skymodel, \
                      log=ms+'_coramptec-c'+str(c)+'.log', cmd_type='BBS')
            s.run(check=True)
        
        logging.info('Restoring WEIGHT_SPECTRUM before imging...')
        s.add('taql "update '+concat_ms+' set WEIGHT_SPECTRUM = WEIGHT_SPECTRUM_ORIG"', log='taql-restweights-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
    
        # join RR and LL
        logging.info('Reconstructing polarizations...')
        for ms in mss_orig:
            msRR = ms.replace('.MS','-RR.MS')
            s.add('taql "update '+ms+', '+msRR+' as rr set CORRECTED_DATA[,0]=rr.CORRECTED_DATA[,0]"', log=ms+'_taql1-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
        for ms in mss_orig:
            msLL = ms.replace('.MS','-LL.MS')
            s.add('taql "update '+ms+', '+msLL+' as ll set CORRECTED_DATA[,3]=ll.CORRECTED_DATA[,3]"', log=ms+'_taql1-c'+str(c)+'.log', cmd_type='general', log_append=True)
        s.run(check=True)
        # do also FLAG in case losoto made them
        for ms in mss_orig:
            msRR = ms.replace('.MS','-RR.MS')
            s.add('taql "update '+ms+', '+msRR+' as rr set FLAG[,0]=rr.FLAG[,0]"', log=ms+'_taql2-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
        for ms in mss_orig:
            msLL = ms.replace('.MS','-LL.MS')
            s.add('taql "update '+ms+', '+msLL+' as ll set FLAG[,3]=ll.FLAG[,3]"', log=ms+'_taql2-c'+str(c)+'.log', cmd_type='general', log_append=True)
        s.run(check=True)

        # after columns creation
        logging.info('Concatenating TCs...')
        check_rm(concat_ms_orig+'*')
        pt.msutil.msconcat(mss_orig, concat_ms_orig, concatTime=False)
    
        ###################################################################################################################
        # concat all TCs in one MS - group*_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected)
    
        # clean mask clean (cut at 8k lambda) - MODEL_DATA updated
        logging.info('Cleaning (cycle: '+str(c)+')...')
        imagename = 'img/wide-'+str(c)
        s.add('wsclean_1.8 -reorder -name ' + imagename + ' -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
                -scale 5arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 8000 -mgain 0.85 '+concat_ms_orig, \
                log='wscleanA-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
        make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask')
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
                   params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/1RXSJ0603_LBA/tooth_mask.crtf', 'setTo':1}, log='casablank-c'+str(c)+'.log')
        s.run(check=True)
        logging.info('Cleaning low resolution (cycle: '+str(c)+')...')
        s.add('wsclean_1.8 -reorder -name ' + imagename + '-masked -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
                -scale 5arcsec -weight briggs 0.0 -niter 20000 -update-model-required -maxuv-l 8000 -mgain 0.85 -casamask '+imagename+'.newmask '+concat_ms_orig, \
                log='wscleanB-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
       
        logging.info('Moving MODEL_DATA to MODEL_DATA_HIGHRES...')
        s.add('taql "update '+concat_ms_orig+' set MODEL_DATA_HIGHRES = MODEL_DATA"', log='taql1-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
    
        ####################################################################
        # FAST VERSION (no low-res)
        #continue
        ####################################################################
    
        ############################################################################################################
        # Subtract model from all TCs - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logging.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        s.add('taql "update '+concat_ms_orig+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql3-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)

        # reclean low-resolution
        logging.info('Cleaning low resolution (cycle: '+str(c)+')...')
        imagename = 'img/wide-lr-'+str(c)
        s.add('wsclean_1.8 -reorder -name ' + imagename + ' -size 4000 4000 -mem 30 -j '+str(s.max_processors)+'\
                -scale 15arcsec -weight briggs 0.0 -niter 50000 -no-update-model-required -maxuv-l 2500 -mgain 0.85 '+concat_ms_orig, \
                log='wscleanA-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
        make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask', threshpix=6) # a bit higher treshold
        logging.info('Cleaning low resolution with mask (cycle: '+str(c)+')...')
        s.add('wsclean_1.8 -reorder -name ' + imagename + '-masked -size 4000 4000 -mem 30 -j '+str(s.max_processors)+' \
                -scale 15arcsec -weight briggs 0.0 -niter 10000 -update-model-required -maxuv-l 2500 -mgain 0.85 -casamask '+imagename+'.newmask '+concat_ms_orig, \
                log='wscleanB-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)

        ###############################################################################################################
        # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (empty)
        logging.info('Subtracting low-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        s.add('taql "update '+concat_ms_orig+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql4-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)

        # Flag on residuals
        logging.info('Flagging residuals...')
        for ms in mss_orig:
            s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms, \
                    log=ms+'_flag-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)
    
        # Concat models
        logging.info('Adding model data columns (MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA)...')
        s.add('taql "update '+concat_ms_orig+' set MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA"', log='taql5-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
    
    # Perform a final clean to create an inspection image of SUBTRACTED_DATA which should be very empty
    logging.info('Empty cleaning...')
    s.add('taql "update '+concat_ms_orig+' set SUBTRACTED_DATA = CORRECTED_DATA"', log='taql_sub.log', cmd_type='general')
    s.run(check=True)
    imagename = 'img/empty'
    s.add('wsclean_1.8 -reorder -name ' + imagename + ' -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
            -scale 5arcsec -weight briggs 0.0 -niter 1 -no-update-model-required -maxuv-l 8000 -mgain 0.85 -datacolumn SUBTRACTED_DATA '+concat_ms_orig, \
            log='wsclean-empty.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    
    # Copy last *model
    logging.info('Copying models/images...')
    os.system('mv img/wide-'+str(c)+'-masked-model.fits self/models/wide-g'+g+'.model')
    os.system('mv img/wide-lr-'+str(c)+'-masked-model.fits self/models/wide-lr-g'+g+'.model')
    # Copy images
    [ os.system('mv img/wide-'+str(c)+'.newmask self/images/g'+g) for c in xrange(niter) ]
    [ os.system('mv img/wide-lr-'+str(c)+'.newmask self/images/g'+g) for c in xrange(niter) ]
    [ os.system('mv img/wide-'+str(c)+'-image.fits self/images/g'+g) for c in xrange(niter) ]
    [ os.system('mv img/wide-lr-'+str(c)+'-image.fits self/images/g'+g) for c in xrange(niter) ]
    [ os.system('mv img/wide-'+str(c)+'-masked-image.fits self/images/g'+g) for c in xrange(niter) ]
    [ os.system('mv img/wide-lr-'+str(c)+'-masked-image.fits self/images/g'+g) for c in xrange(niter) ]
    os.system('mv img/empty-image.fits self/images/g'+g)
    os.system('mv *log '+group)

# TODO: add final imaging with all SBs

logging.info("Done.")
