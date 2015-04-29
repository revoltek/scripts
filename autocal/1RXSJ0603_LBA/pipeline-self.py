#!/usr/bin/python
# perform self-calibration on a group of SBs concatenated in TCs. Script must be run in dir with MS.
# number/chan in MS are flexible but the must be concatenable (same chans/freq!)
# Input:
# TCs are blocks of SBs should have calibrator corrected (a+p) data in DATA (beam not applied).
# file format of TCs is: group#_TC###.MS.
# Output:
# TCs with selfcal corrected source subtracted data in SUBTRACTED_DATA
# instrument tables contain gain (slow) + fast (scalarphase+TEC) solutions
# last high/low resolution models are copied in the "self/models" dir
# last high/low resolution images are copied in the "self/images" dir

skymodel = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/toothbrush.GMRT150.skymodel'
max_threads = 40

#######################################################################################

import sys, os, glob, re
import numpy as np
from lofar import bdsm
import pyrap.tables as pt
import lsmtool
from lib_pipeline import *
from make_mask import make_mask

set_logger()

# here an image+model for each group will be saved
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/models'): os.makedirs('self/models')

for group in sorted(glob.glob('group*')):

    mss = sorted(glob.glob(group+'/group*_TC*.MS'))
    g = str(re.findall(r'\d+', mss[0])[0])
    logging.info('Working on group: '+g+'...')
    concat_ms = group+'/concat.MS'
    
    ################################################################################################
    # Clear
    logging.info('Cleaning...')
    check_rm(group+'/*log *log')
    check_rm(group+'/plot* plot')
    check_rm(group+'/*h5 *h5')
    check_rm('*last')
    check_rm('img')
    os.makedirs('img')
    
    #################################################################################################
    # create a fake parmdb to be used later for merging slow-amp and fast-phase parmdbs
    logging.info('Creating fake parmdb...')
    cmds = []
    for ms in mss:
        cmds.append('calibrate-stand-alone -f --parmdb-name instrument_empty '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-fakeparmdb.parset '+skymodel+' > '+ms+'_fakeparmdb.log 2>&1')
    thread_cmd(cmds, max_threads=max_threads)
    
    ####################################################################################################
    # self-cal cycle
    for i in xrange(3):
        logging.info('Start selfcal cycle: '+str(i))
    
        # TEST for circular
        # separate LL and RR
        #msLL = ms.replace('.MS','-LL.MS')
        #if os.path.exists(msLL): os.system('rm -r '+msLL)
        #msRR = ms.replace('.MS','-RR.MS')
        #if os.path.exists(msRR): os.system('rm -r '+msRR)
        #os.system( 'cp -r '+ms+' '+msLL )
        #os.system( 'mv '+ms+' '+msRR )
        #print 'taql \'update '+msRR+' set DATA[,3]=DATA[,0]\''
        #os.system( 'taql \'update '+msRR+' set DATA[,3]=DATA[,0]\'' )
        #print 'taql \'update '+msLL+' set DATA[,0]=DATA[,3]\''
        #os.system( 'taql \'update '+msLL+' set DATA[,0]=DATA[,3]\'' )
    
        if i == 0:
            # calibrate phase-only - group*_TC.MS:DATA (beam: ARRAY_FACTOR) -> group*_TC.MS:CORRECTED_DATA (selfcal phase corrected, beam corrected)
            logging.info('Calibrating phase...')
            cmds = []
            for ms in mss:
                cmds.append('calibrate-stand-alone -f '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-solcor.parset '+skymodel+' > '+ms+'_cal-c'+str(i)+'.log 2>&1')
            thread_cmd(cmds, max_threads=max_threads)
        else:
            # calibrate phase-only - group*_TC.MS:DATA @ MODEL_DATA -> group*_TC.MS:CORRECTED_DATA_PHASE (selfcal phase corrected, beam corrected)
            logging.info('Calibrating phase...')
            cmds = []
            for ms in mss:
                cmds.append('calibrate-stand-alone -f --parmdb-name instrument_csp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-solcor_csp.parset '+skymodel+' > '+ms+'_calpreamp-c'+str(i)+'.log 2>&1')
            thread_cmd(cmds, max_threads=max_threads)
    
            # calibrate amplitude (only solve) - group*_TC.MS:CORRECTED_DATA_PHASE @ MODEL_DATA
            logging.info('Calibrating amplitude...')
            cmds = []
            for ms in mss:
                cmds.append('calibrate-stand-alone -f --parmdb-name instrument_amp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-sol_amp.parset '+skymodel+' > '+ms+'_calamp-c'+str(i)+'.log 2>&1')
            thread_cmd(cmds, max_threads=max_threads)
    
            # merge parmdbs
            logging.info('Merging instrument tables...')
            for ms in mss:
                check_rm(ms+'/instrument')
                merge_parmdb(ms+'/instrument_amp', ms+'/instrument_csp', ms+'/instrument_empty', ms+'/instrument')
    
            ########################################################
            # LoSoTo Amp rescaling
            os.makedirs('plot')
            check_rm('globaldb')
            os.makedirs('globaldb')
            for num, ms in enumerate(mss):
                os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
                if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
            h5parm = 'global-c'+str(i)+'.h5'
            os.system('H5parm_importer.py -v '+h5parm+' globaldb > losoto-c'+str(i)+'.log 2>&1')
            os.system('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/losoto.parset >> losoto-c'+str(i)+'.log 2>&1')
            os.system('H5parm_exporter.py -v -c '+h5parm+' globaldb >> losoto-c'+str(i)+'.log 2>&1')
            for num, ms in enumerate(mss):
                check_rm(ms+'/instrument')
                os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
            os.system('mv plot '+group+'/plot-c'+str(i))
            os.system('mv '+h5parm+' '+group)
    
            # correct amplitude - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase+amp corrected, beam corrected)
            logging.info('Correcting amplitude...')
            cmds = []
            for ms in mss:
                cmds.append('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-cor_ampcsp.parset '+skymodel+' > '+ms+'_corampcsp-c'+str(i)+'.log 2>&1')
            thread_cmd(cmds, max_threads=max_threads)
    
        # TEST for circular
        # join RR and LL
        #if os.path.exists(ms): os.system('rm -r '+ms)
        #os.system('mv '+msRR+' '+ms)
        #print 'taql \'update '+ms+', '+msLL+' as ll set DATA[3,]=ll.DATA[3,]\''
        #os.system('taql \'update '+ms+', '+msLL+' as ll set DATA[3,]=ll.DATA[3,]\'')
        #os.system('rm -r '+msLL)
    
        ###################################################################################################################
        # concat all TCs in one MS - group*_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected)
        # NOTE: observations must have the same channels i.e. same SBs freq!
        # NOTE: virtual concat or casa concat both fails if second observation SB is the first of the concat list (mss are now sorted), very strange bug
        check_rm(concat_ms+'*')
        logging.info('Concatenating TCs...')
        pt.msutil.msconcat(mss, concat_ms, concatTime=False)
    
        # clean mask clean (cut at 8k lambda) - MODEL_DATA updated
        logging.info('Cleaning 1...')
        imagename = 'img/wide-'+str(i)
        os.system('/opt/cep/WSClean/wsclean-1.7/build/wsclean -reorder -name ' + imagename + ' -size 5000 5000 \
                -scale 5arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -no-update-model-required -maxuv-l 8000 '+concat_ms+' > wscleanA-c'+str(i)+'.log 2>&1')
        make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask')
        run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_blank.py', params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/1RXSJ0603_LBA/tooth_mask.crtf', 'setTo':1})
        logging.info('Cleaning 2...')
        os.system('/opt/cep/WSClean/wsclean-1.7/build/wsclean -reorder -name ' + imagename + '-masked -size 5000 5000 \
                -scale 5arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -update-model-required -maxuv-l 8000 -casamask '+imagename+'.newmask '+concat_ms+'  > wscleanB-c'+str(i)+'.log 2>&1')
    
        # copy model column
        if i == 0:
            cmds = []
            logging.info('Creating MODEL_DATA_HIGHRES...')
            for ms in mss:
                cmds.append('addcol2ms.py -i '+ms+' -o MODEL_DATA_HIGHRES')
            thread_cmd(cmds, max_threads=max_threads)
            cmds = []
        logging.info('Moving MODEL_DATA to MODEL_DATA_HIGHRES...')
        for ms in mss:
            cmds.append('taql "update '+ms+' set MODEL_DATA_HIGHRES = MODEL_DATA"')
        thread_cmd(cmds, max_threads=max_threads)
    
        ####################################################################
        # FAST VERSION (no low-res)
        #continue
        ####################################################################
    
        ####################################################################################################################################################
        # Subtract model from all TCs - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        cmds = []
        logging.info('Subtracting high-res model...')
        for ms in mss:
            cmds.append('taql "update '+ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"')
        thread_cmd(cmds, max_threads=max_threads)
    
        # concat all TCs in one MS - group*_TC*.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
        check_rm(concat_ms+'*')
        logging.info('Concatenating TCs...')
        pt.msutil.msconcat(mss, concat_ms, concatTime=False)
    
        # reclean low-resolution
        logging.info('Cleaning low resolution 1...')
        imagename = 'img/wide-lr-'+str(i)
        os.system('/opt/cep/WSClean/wsclean-1.7/build/wsclean -reorder -name ' + imagename + ' -size 4000 4000 \
                -scale 15arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -no-update-model-required -maxuv-l 2500 '+concat_ms+' > wscleanA-lr-c'+str(i)+'.log 2>&1')
        make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask')
        logging.info('Cleaning low resolution 2...')
        os.system('/opt/cep/WSClean/wsclean-1.7/build/wsclean -reorder -name ' + imagename + '-masked -size 4000 4000 \
                -scale 15arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -update-model-required -maxuv-l 2500 -casamask '+imagename+'.newmask '+concat_ms+' > wscleanB-lr-c'+str(i)+'.log 2>&1')
    
        # Concat models
        cmds = []
        logging.info('Adding model data columns (MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA)...')
        for ms in mss:
            cmds.append('taql "update '+ms+' set MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA"')
        thread_cmd(cmds, max_threads=max_threads)
    
    # Final subtract of the best model - group*_TC*.MS:DATA - MODEL_DATA -> group*_TC*.MS:SUBTRACTED_DATA (not corrected data - all source subtracted, beam corrected, circular)
    cmds = []
    logging.info('Final subtraction...')
    for ms in mss:
        cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-subfinal.parset '+skymodel+' > '+ms+'_final-sub.log 2>&1')
    thread_cmd(cmds, max_threads=max_threads)
    
    # Copy last *model and *image
    logging.info('Copying models/images...')
    os.system('mv img/wide-'+str(i)+'-masked-model.fits self/models/wide-g'+g+'.model')
    os.system('mv img/wide-lr-'+str(i)+'-masked-model.fits self/models/wide-lr-g'+g+'.model')
    os.system('mv img/wide-'+str(i)+'-masked-image.fits self/images/wide-g'+g+'.image')
    os.system('mv img/wide-lr-'+str(i)+'-masked-image.fits self/images/wide-lr-g'+g+'.image')
    os.system('mv *log '+group)

logging.info("Done.")
