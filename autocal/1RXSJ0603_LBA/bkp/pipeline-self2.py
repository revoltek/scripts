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
fake_skymodel = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/toothbrush.fakemodel.skymodel'
max_threads = 16

#######################################################################################

import sys, os, glob, re
import numpy as np
from lofar import bdsm
import pyrap.tables as pt
import lsmtool
from lib_pipeline import *

set_logger()

mss = sorted(glob.glob('group*_TC*.MS'))
g = str(re.findall(r'\d+', mss[0])[0])
logging.info('Working on group: '+g+'...')

#############################################################################################
# Clear
logging.info('Cleaning...')
check_rm('concat*')
check_rm('*log')
check_rm('*h5')
check_rm('*last')
check_rm('img')
os.makedirs('img')
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/models'): os.makedirs('self/models')

#################################################################################################
# create a fake parmdb to be used later for merging slow-amp and fast-phase parmdbs
logging.info('Creating fake parmdb...')
cmds = []
for ms in mss:
    cmds.append('calibrate-stand-alone -f --parmdb-name instrument_empty '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/bbs-fakeparmdb.parset '+fake_skymodel+' > '+ms+'_fakeparmdb.log 2>&1')
thread_cmd(cmds, max_threads=max_threads)

# TODO: a run of AOflagger on concatenated data

####################################################################################################
# self-cal cycle
for i in xrange(2):
    logging.info('Start selfcal cycle: '+str(i))

#    if i != 0: skymodel = 'widefield-comb_g'+g+'.skymodel'

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
            cmds.append('calibrate-stand-alone -f '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/bbs-solcor.parset '+skymodel+' > '+ms+'_cal-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads=max_threads)
    else:
        # calibrate phase-only - group*_TC.MS:DATA @ MODEL_DATA -> group*_TC.MS:CORRECTED_DATA_PHASE (selfcal phase corrected, beam corrected)
        logging.info('Calibrating phase...')
        cmds = []
        for ms in mss:
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_csp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/bbs-solcor_csp.parset '+skymodel+' > '+ms+'_calpreamp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads=max_threads)

        # calibrate amplitude (only solve) - group*_TC.MS:CORRECTED_DATA_PHASE @ MODEL_DATA
        logging.info('Calibrating amplitude...')
        cmds = []
        for ms in mss:
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_amp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/bbs-sol_amp.parset '+skymodel+' > '+ms+'_calamp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads=max_threads)

        ########################################################
        # LoSoTo Amp rescaling
        for ms in mss:
            h5parm = ms.replace('.MS','.h5')
            check_rm(h5parm)
            os.system('H5parm_importer.py -v -i instrument_amp '+h5parm+' '+ms+' > '+ms+'_losoto-c'+str(i)+'.log 2>&1')
            os.system('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/losoto.parset >> '+ms+'_losoto-c'+str(i)+'.log 2>&1')
            os.system('H5parm_exporter.py -v -i instrument_amp -c '+h5parm+' '+ms+' >> '+ms+'_losoto-c'+str(i)+'.log 2>&1')
            check_rm(ms+'/instrument_amp')
            os.system('mv '+ms+'/sol000_instrument_amp'+' '+ms+'/instrument_amp')

        # merge parmdbs
        logging.info('Merging instrument tables...')
        for ms in mss:
            check_rm(ms+'/instrument')
            merge_parmdb(ms+'/instrument_amp', ms+'/instrument_csp', ms+'/instrument_empty', ms+'/instrument')

        # correct amplitude - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase+amp corrected, beam corrected)
        logging.info('Correcting amplitude...')
        cmds = []
        for ms in mss:
            cmds.append('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/bbs-cor_ampcsp.parset '+skymodel+' > '+ms+'_corampcsp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads=max_threads)

    # TEST for circular
    # join RR and LL
    #if os.path.exists(ms): os.system('rm -r '+ms)
    #os.system('mv '+msRR+' '+ms)
    #print 'taql \'update '+ms+', '+msLL+' as ll set DATA[3,]=ll.DATA[3,]\''
    #os.system('taql \'update '+ms+', '+msLL+' as ll set DATA[3,]=ll.DATA[3,]\'')
    #os.system('rm -r '+msLL)

    ##############################################################################################################################################
    # concat all TCs in one MS - group*_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected)
    # NOTE: observations must have the same channels i.e. same SBs freq!
    # NOTE: virtual concat or casa concat both fails if second observation SB is the first of the concat list (mss are now sorted), very strange bug
    check_rm('concat.MS*')
    logging.info('Concatenating TCs...')
    pt.msutil.msconcat(mss, 'concat.MS', concatTime=False)
#    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/casa_concat.py '+msregexp+' concat.MS > split1-c'+str(i)+'.log 2>&1')

    # clean mask clean
    # TODO: taper at 7-8klambda
    logging.info('Cleaning...')
    imagename = 'img/wide-'+str(i)
    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/casa_clean.py concat.MS '+imagename+'-casa > casaclean.log 2>&1')
    os.system('/opt/cep/WSClean/wsclean-1.7/build/wsclean -reorder -name ' + imagename + '-WSC -size 5000 5000 \
            -scale 5arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -no-update-model-required concat.MS > waclean.log 2>&1')
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager.sh concat.MS '+imagename+'-awimager > awimager-old.log 2>&1')
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager-new.sh concat.MS '+imagename+'-awimager > awimager-new.log 2>&1')
    sys.exit('clean test')
    os.system('make_mask.py -p 6 -i 4 -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/tooth.mask -m '+imagename+'.newmask '+imagename+'.restored.corr > make_mask-c'+str(i)+'.log  2>&1')
    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/casa_clean.py concat.MS '+imagename+'-masked '+imagename+'.newmask > cleanB-c'+str(i)+'.log 2>&1')
#    os.system('/opt/cep/WSClean/wsclean-1.7/build/wsclean -reorder -name ' + imagename + '-masked -size 5000 5000 \
#            -scale 5arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -no-update-model-required -casamask '+imagename+'.newmask concat.MS')

    # copy model column
    cmds = []
    logging.info('Creating MODEL_DATA_HIGHRES...')
    for ms in mss:
        cmds.append('addcol2ms.py -i '+ms+' -o MODEL_DATA_HIGHRES')
    thread_cmd(cmds, max_threads=max_threads)

    ####################################################################
    # FAST VERSION (no low-res)
    #continue
    ####################################################################

    ####################################################################################################################################################
    # Subtract model from all TCs - concat.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
    cmds = []
    logging.info('Subtracting high-res model...')
    for ms in mss:
            cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/casa_uvsub.py '+ms+' > '+ms+'_uvsub-c'+str(i)+'.log 2>&1')
    thread_cmd(cmds, max_threads=max_threads)

    # concat all TCs in one MS - group*_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
    check_rm('concat.MS*')
    logging.info('Concatenating TCs...')
    pt.msutil.msconcat(mss, 'concat.MS', concatTime=False)
    #os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/casa_concat.py '+msregexp+' concat.MS > split2-c'+str(i)+'.log 2>&1')

    # reclean low-resolution
    logging.info('Cleaning low resolution...')
    imagename = 'img/wide-lr-'+str(i)
    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/casa_clean-lr.py concat.MS '+imagename+' > cleanA-lr-c'+str(i)+'.log 2>&1')
#    os.system('/opt/cep/WSClean/wsclean-1.7/build/wsclean -reorder -name ' + imagename + ' -size 2500 2500 \
#            -scale 15arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -no-update-model-required -maxuv-l 3500 concat.MS')
    os.system('make_mask.py -p 6 -i 4 -m '+imagename+'.newmask '+imagename+'.restored.corr > make_mask-lr-c'+str(i)+'.log  2>&1')
    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/casa_clean-lr.py concat.MS '+imagename+'-masked '+imagename+'.newmask > cleanB-lr-c'+str(i)+'.log 2>&1')
#    os.system('/opt/cep/WSClean/wsclean-1.7/build/wsclean -reorder -name ' + imagename + '-masked -size 2500 2500 \
#            -scale 15arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -no-update-model-required -maxuv-l 3500 -casamask '+imagename+'.newmask concat.MS')

    # Concat models
    cmds = []
    logging.info('Adding model data columns...')
    for ms in mss:
        cmds.append('taql "update '+ms+' set MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA"')
    thread_cmd(cmds, max_threads=max_threads)

# Final subtract of the best model - group*_TC.MS:DATA - MODEL_DATA -> group*_TC.MS:SUBTRACTED_DATA (all source subtracted, beam corrected, circular)
cmds = []
logging.info('Final subtraction...')
for ms in mss:
    cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self2/bbs-subfinal.parset > '+ms+'_final-sub.log 2>&1')
thread_cmd(cmds, max_threads=max_threads)

# Copy .model[tt*] images
[os.system('cp -r '+model+' self/models/wide_g'+g+'.'+model.split('.')[1]) for model in glob.glob('img/wide-'+str(i)+'.model*')]
[os.system('cp -r '+model+' self/models/wide-lr_g'+g+'.'+model.split('.')[1]) for model in glob.glob('img/wide-lr-'+str(i)+'.model*')]
[os.system('cp -r '+image+' self/images/wide_g'+g+'.'+model.split('.')[1]) for image in glob.glob('img/wide-'+str(i)+'.image*')]
[os.system('cp -r '+image+' self/images/wide-lr_g'+g+'.'+model.split('.')[1]) for image in glob.glob('img/wide-lr-'+str(i)+'.image*')]

logging.info("Done.")
