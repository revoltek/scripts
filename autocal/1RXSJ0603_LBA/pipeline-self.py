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

def losoto(c, mss, mssavg, g, parset):
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

    for num, ms in enumerate(mssavg):
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


# here an image+model for each group will be saved
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/models'): os.makedirs('self/models')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')

for group in sorted(glob.glob('group*'))[::-1]:

    mss_orig = sorted(glob.glob(group+'/group*_TC*[0-9].MS'))
    concat_ms = group+'/concat.MS'
    g = str(re.findall(r'\d+', mss_orig[0])[0])
    logging.info('Working on group: '+g+'...')
    
    ################################################################################################
    # Clear
    logging.info('Cleaning...')
    check_rm(group+'/*-BLavg.MS')
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

        s.add('taql "update '+msRR+' set DATA[,3]=DATA[,0]"', log=ms+'_init-taql.log', cmd_type='general')
        s.add('taql "update '+msLL+' set DATA[,0]=DATA[,3]"', log=ms+'_init-taql.log', cmd_type='general', log_append=True)
    s.run(check=True)

    mss = sorted(glob.glob(group+'/group*_TC*[0-9]-*.MS'))
    mssrr = sorted(glob.glob(group+'/group*_TC*[0-9]-RR.MS'))
    mssll = sorted(glob.glob(group+'/group*_TC*[0-9]-LL.MS'))

    #################################################################################################
    # Smooth
    logging.info('Smoothing...')
    for ms in mss:
        s.add('BLavg.py -m '+ms, log=ms+'_smooth.log', cmd_type='python')
    s.run(check=True)
    mssavg = sorted(glob.glob(group+'/group*_TC*-BLavg.MS'))
    mssavgrr = sorted(glob.glob(group+'/group*_TC*[0-9]-RR-BLavg.MS'))
    mssavgll = sorted(glob.glob(group+'/group*_TC*[0-9]-LL-BLavg.MS'))

    ####################################################################################################
    # Self-cal cycle
    for i in xrange(niter):
        logging.info('Start selfcal cycle: '+str(i))
   
        if i == 0:
            # calibrate phase-only - group*_TC.MS:DATA (beam: ARRAY_FACTOR) -> group*_TC.MS:CORRECTED_DATA (selfcal phase corrected, beam corrected)
            logging.info('Calibrating phase...')
            for ms in mssavg:
                s.add('calibrate-stand-alone -f '+ms+' '+parset_dir+'/bbs-sol_tec.parset '+skymodel, \
                      log=ms+'_soltec-c'+str(i)+'.log', cmd_type='BBS')
            s.run(check=True)
            losoto(str(i)+'rr', mssrr, mssavgrr, g, parset_dir+'/losoto-plot.parset')
            losoto(str(i)+'ll', mssll, mssavgll, g, parset_dir+'/losoto-plot.parset')
            logging.info('Correcting phase...')
            for ms in mss:
                s.add('calibrate-stand-alone '+ms+' '+parset_dir+'/bbs-cor_tec.parset '+skymodel, \
                log=ms+'_cortec-c'+str(i)+'.log', cmd_type='BBS')
            s.run(check=True)
        else:
            # copy FLAG and MODEL_DATA in avgBL - group*_TC-avgBL.MS:MODEL_DATA = group*_TC.MS:MODEL_DATA
            logging.info('Copying FLAG and MODEL_DATA to the averaged datasets...')
            for ms in mss_orig:
                msavg = ms.replace('.MS','-RR-BLavg.MS')
                s.add('taql "update '+msavg+', '+ms+' as orig set MODEL_DATA[,0]=orig.MODEL_DATA[,0]" && \
                       taql "update '+msavg+', '+ms+' as orig set MODEL_DATA[,3]=orig.MODEL_DATA[,0]" && \
                       taql "update '+msavg+', '+ms+' as orig set FLAG[,0]=orig.FLAG[,0]" && \
                       taql "update '+msavg+', '+ms+' as orig set FLAG[,3]=orig.FLAG[,3]"', \
                       log=ms+'taql_copymodel-c'+str(i)+'.log', cmd_type='general')
            s.run(check=True)
            for ms in mss_orig:
                msavg = ms.replace('.MS','-LL-BLavg.MS')
                s.add('taql "update '+msavg+', '+ms+' as orig set MODEL_DATA[,0]=orig.MODEL_DATA[,3]" && \
                       taql "update '+msavg+', '+ms+' as orig set MODEL_DATA[,3]=orig.MODEL_DATA[,3]" && \
                       taql "update '+msavg+', '+ms+' as orig set FLAG[,0]=orig.FLAG[,3]" && \
                       taql "update '+msavg+', '+ms+' as orig set FLAG[,3]=orig.FLAG[,3]"', \
                       log=ms+'taql_copymodel-c'+str(i)+'.log', cmd_type='general')
            s.run(check=True)

            # calibrate phase-only - group*_TC.MS:DATA @ MODEL_DATA -> group*_TC.MS:CORRECTED_DATA_PHASE (selfcal phase corrected, beam corrected)
            logging.info('Calibrating phase...')
            for ms in mssavg:
                s.add('calibrate-stand-alone -f --parmdb-name instrument_tec '+ms+' '+parset_dir+'/bbs-solcor_tec.parset '+skymodel, \
                      log=ms+'_calpreamp-c'+str(i)+'.log', cmd_type='BBS')
            s.run(check=True)

            # TODO: problem here as TEC are applied to smoothed data!!!
    
            # calibrate amplitude (only solve) - group*_TC.MS:CORRECTED_DATA_PHASE @ MODEL_DATA
            logging.info('Calibrating amplitude...')
            for ms in mssavg:
                s.add('calibrate-stand-alone -f --parmdb-name instrument_amp '+ms+' '+parset_dir+'/bbs-sol_amp.parset '+skymodel, \
                      log=ms+'_calamp-c'+str(i)+'.log', cmd_type='BBS')
            s.run(check=True)

            # merge parmdbs
            logging.info('Merging instrument tables...')
            for ms in mssavg:
                merge_parmdb(ms+'/instrument_tec', ms+'/instrument_amp', ms+'/instrument', clobber=True)
    
            ########################################################
            # LoSoTo Amp rescaling
            losoto(str(i)+'rr', mssrr, mssavgrr, g, parset_dir+'/losoto.parset')
            losoto(str(i)+'ll', mssll, mssavgll, g, parset_dir+'/losoto.parset')
        
            # correct - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase+amp corrected, beam corrected)
            logging.info('Correcting...')
            for ms in mss:
                s.add('calibrate-stand-alone '+ms+' '+parset_dir+'/bbs-cor_amptec.parset '+skymodel, \
                      log=ms+'_coramptec-c'+str(i)+'.log', cmd_type='BBS')
            s.run(check=True)
    
        # join RR and LL
        logging.info('Reconstructing polarizations...')
        for ms in mss_orig:
            msRR = ms.replace('.MS','-RR.MS')
            s.add('taql "update '+ms+', '+msRR+' as rr set CORRECTED_DATA[,0]=rr.CORRECTED_DATA[,0]"', log=ms+'_taql-c'+str(i)+'.log', cmd_type='general')
        s.run(check=True)
        for ms in mss_orig:
            msLL = ms.replace('.MS','-LL.MS')
            s.add('taql "update '+ms+', '+msLL+' as ll set CORRECTED_DATA[,3]=ll.CORRECTED_DATA[,3]"', log=ms+'_taql-c'+str(i)+'.log', cmd_type='general', log_append=True)
        s.run(check=True)

        # after columns creation
        logging.info('Concatenating TCs...')
        check_rm(concat_ms+'*')
        pt.msutil.msconcat(mss_orig, concat_ms, concatTime=False)
    
        ###################################################################################################################
        # concat all TCs in one MS - group*_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected)
    
        # clean mask clean (cut at 8k lambda) - MODEL_DATA updated
        logging.info('Cleaning (cycle: '+str(i)+')...')
        imagename = 'img/wide-'+str(i)
        s.add('wsclean_1.8 -reorder -name ' + imagename + ' -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
                -scale 5arcsec -weight briggs 0.0 -niter 100000 -mgain 1 -no-update-model-required -maxuv-l 8000 -mgain 0.85 '+concat_ms, \
                log='wscleanA-c'+str(i)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
        make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask')
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
                   params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/1RXSJ0603_LBA/tooth_mask.crtf', 'setTo':1}, log='casablank-c'+str(i)+'.log')
        s.run(check=True)
        logging.info('Cleaning low resolution (cycle: '+str(i)+')...')
        s.add('wsclean_1.8 -reorder -name ' + imagename + '-masked -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
                -scale 5arcsec -weight briggs 0.0 -niter 20000 -mgain 1 -update-model-required -maxuv-l 8000 -mgain 0.85 -casamask '+imagename+'.newmask '+concat_ms, \
                log='wscleanB-c'+str(i)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
       
        logging.info('Moving MODEL_DATA to MODEL_DATA_HIGHRES...')
        s.add('taql "update '+concat_ms+' set MODEL_DATA_HIGHRES = MODEL_DATA"', log='taql1-c'+str(i)+'.log', cmd_type='general')
        s.run(check=True)
    
        ####################################################################
        # FAST VERSION (no low-res)
        #continue
        ####################################################################
    
        ############################################################################################################
        # Subtract model from all TCs - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        logging.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql2-c'+str(i)+'.log', cmd_type='general')
        s.run(check=True)

        # reclean low-resolution
        logging.info('Cleaning low resolution (cycle: '+str(i)+')...')
        imagename = 'img/wide-lr-'+str(i)
        s.add('wsclean_1.8 -reorder -name ' + imagename + ' -size 4000 4000 -mem 30 -j '+str(s.max_processors)+'\
                -scale 15arcsec -weight briggs 0.0 -niter 50000 -mgain 1 -no-update-model-required -maxuv-l 2500 -mgain 0.85 '+concat_ms, \
                log='wscleanA-lr-c'+str(i)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
        make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask', threshpix=6) # a bit higher treshold
        logging.info('Cleaning low resolution with mask (cycle: '+str(i)+')...')
        s.add('wsclean_1.8 -reorder -name ' + imagename + '-masked -size 4000 4000 -mem 30 -j '+str(s.max_processors)+' \
                -scale 15arcsec -weight briggs 0.0 -niter 10000 -mgain 1 -update-model-required -maxuv-l 2500 -mgain 0.85 -casamask '+imagename+'.newmask '+concat_ms, \
                log='wscleanB-lr-c'+str(i)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)

        ###############################################################################################################
        # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (empty)
        logging.info('Subtracting low-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
        s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql3-c'+str(i)+'.log', cmd_type='general')
        s.run(check=True)

        # Flag on residuals
        logging.info('Flagging residuals...')
        for ms in mss_orig:
            s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms, \
                    log=ms+'_flag-c'+str(i)+'.log', cmd_type='NDPPP')
        s.run(check=True)
    
        # Concat models
        logging.info('Adding model data columns (MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA)...')
        s.add('taql "update '+concat_ms+' set MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA"', log='taql4-c'+str(i)+'.log', cmd_type='general')
        s.run(check=True)
    
    # Perform a final clean to create an inspection image of SUBTRACTED_DATA which should be very empty
    logging.info('Empty cleaning...')
    s.add('taql "update '+concat_ms+' set SUBTRACTED_DATA = CORRECTED_DATA"', log='taql5.log', cmd_type='general')
    s.run(check=True)
    imagename = 'img/empty'
    s.add('wsclean_1.8 -reorder -name ' + imagename + ' -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
            -scale 5arcsec -weight briggs 0.0 -niter 1 -mgain 1 -no-update-model-required -maxuv-l 8000 -mgain 0.85 -datacolumn SUBTRACTED_DATA '+concat_ms, \
            log='wscleanA-c'+str(i)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    
    # Copy last *model
    logging.info('Copying models/images...')
    os.system('mv img/wide-'+str(i)+'-masked-model.fits self/models/wide-g'+g+'.model')
    os.system('mv img/wide-lr-'+str(i)+'-masked-model.fits self/models/wide-lr-g'+g+'.model')
    # Copy images
    [ os.system('mv img/wide-'+str(i)+'.newmask self/images/g'+g) for i in xrange(niter) ]
    [ os.system('mv img/wide-lr-'+str(i)+'.newmask self/images/g'+g) for i in xrange(niter) ]
    [ os.system('mv img/wide-'+str(i)+'-image.fits self/images/g'+g) for i in xrange(niter) ]
    [ os.system('mv img/wide-lr-'+str(i)+'-image.fits self/images/g'+g) for i in xrange(niter) ]
    [ os.system('mv img/wide-'+str(i)+'-masked-image.fits self/images/g'+g) for i in xrange(niter) ]
    [ os.system('mv img/wide-lr-'+str(i)+'-masked-image.fits self/images/g'+g) for i in xrange(niter) ]
    os.system('mv img/empty-image.fits self/images/g'+g)
    os.system('mv *log '+group)

# TODO: add final imaging with all SBs

logging.info("Done.")
