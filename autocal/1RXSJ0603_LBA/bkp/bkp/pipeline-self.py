#!/usr/bin/python
# perform self-calibration on a group of SBs concatenated in TCs. Script must be run in dir with MS.
# number/chan in MS are flexible but the must be concatenable (same chans/freq!)
# Input:
# TCs are blocks of SBs should have calibrator corrected (a+p) data in DATA (beam not applied).
# file format of TCs is: group#_TC###.MS.
# Output:
# TCs with selfcal corrected source subtracted data in SUBTRACTED_DATA
# instrument tables contain gain (slow) + fast (scalarphase+TEC) solutions

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
check_rm('*skymodel')
check_rm('*last')
check_rm('img')
os.system('mkdir img')

#################################################################################################
# create a fake parmdb to be used later for merging slow-amp and fast-phase parmdbs
logging.info('Creating fake parmdb...')
cmds = []
for ms in mss:
    cmds.append('calibrate-stand-alone -f --parmdb-name instrument_empty '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-fakeparmdb.parset '+fake_skymodel+' > '+ms+'_fakeparmdb.log 2>&1')
thread_cmd(cmds, max_threads=max_threads)

# TODO: a run of AOflagger on concatenated data

####################################################################################################
# self-cal cycle
for i in xrange(4):
    logging.info('Start selfcal cycle: '+str(i))

    if i != 0: skymodel = 'widefield-comb_g'+g+'.skymodel'

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

    if i < 2:
        # calibrate phase-only - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase corrected, beam corrected)
        logging.info('Calibrating phase...')
        cmds = []
        for ms in mss:
            cmds.append('calibrate-stand-alone -f '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-solcor.parset '+skymodel+' > '+ms+'_cal-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads=max_threads)
    else:
        # calibrate phase-only - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA_PHASE (selfcal phase corrected, beam corrupted)
        logging.info('Calibrating phase...')
        cmds = []
        for ms in mss:
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_csp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-solcor_csp.parset '+skymodel+' > '+ms+'_calpreamp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads=max_threads)

        # calibrate amplitude (only solve) - group*_TC.MS:CORRECTED_DATA_PHASE
        logging.info('Calibrating amplitude...')
        cmds = []
        for ms in mss:
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_amp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-sol_amp.parset '+skymodel+' > '+ms+'_calamp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads=max_threads)

        ########################################################
        # LoSoTo Amp rescaling
        for ms in mss:
            h5parm = ms.replace('.MS','.h5')
            check_rm(h5parm)
            os.system('H5parm_importer.py -v -i instrument_amp '+h5parm+' '+ms+' > '+ms+'_losoto-c'+str(i)+'.log 2>&1')
            os.system('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/losoto.parset >> '+ms+'_losoto-c'+str(i)+'.log 2>&1')
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
            cmds.append('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-cor_ampcsp.parset '+skymodel+' > '+ms+'_corampcsp-c'+str(i)+'.log 2>&1')
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
#    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/casa_concat.py '+msregexp+' concat.MS > split1-c'+str(i)+'.log 2>&1')

    # clean mask clean
    logging.info('Cleaning...')
    imagename = 'img/wide-'+str(i)
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager.sh concat.MS '+imagename+' > awimagerA-c'+str(i)+'.log 2>&1')
    os.system('make_mask.py -p 6 -i 4 -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/tooth.mask -m '+imagename+'.mask '+imagename+'.restored.corr > make_mask-c'+str(i)+'.log  2>&1')
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager.sh concat.MS '+imagename+'-masked '+imagename+'.mask > awimagerB-c'+str(i)+'.log 2>&1')

    # making new model
    logging.info('Creating model...')
    img = bdsm.process_image(imagename+'-masked.restored.corr', mean_map='zero', adaptive_rms_box=True, thresh_isl=4, thresh_pix=6,\
            rms_box_bright=(20, 7), rms_box=(120, 40), atrous_do=False, ini_method='curvature', advanced_opts=True, blank_limit=1e-5)
    img.export_image(outfile=imagename+'-masked.BDSMresidual', clobber=True, img_type='gaus_resid')
    img.export_image(outfile=imagename+'-masked.BDSMmodel', clobber=True, img_type='gaus_model')
    img.write_catalog(outfile='widefield_g'+g+'.skymodel', catalog_type='gaul', format='bbs', bbs_patches='source', clobber=True)

    ####################################################################
    # FAST VERSION
    #os.system('cp widefield_g'+g+'.skymodel widefield-comb_g'+g+'.skymodel')
    #continue
    ####################################################################

    ####################################################################################################################################################
    # Subtract model from all TCs - concat.MS:DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
    cmds = []
    logging.info('Subtracting high-res model...')
    for ms in mss:
        if i < 2:
            cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-sub_csp.parset widefield_g'+g+'.skymodel > '+ms+'_sub-c'+str(i)+'.log 2>&1')
        else:
            cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/bbs-sub_ampcsp.parset widefield_g'+g+'.skymodel > '+ms+'_sub-c'+str(i)+'.log 2>&1')
    thread_cmd(cmds, max_threads=max_threads)

    # concat all TCs in one MS - group*_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
    check_rm('concat.MS*')
    logging.info('Concatenating TCs...')
    pt.msutil.msconcat(mss, 'concat.MS', concatTime=False)
    #os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_self/casa_concat.py '+msregexp+' concat.MS > split2-c'+str(i)+'.log 2>&1')

    # reclean low-resolution
    logging.info('Cleaning low resolution...')
    imagename = 'img/wide-lr-'+str(i)
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager-lr.sh concat.MS '+imagename+' > awimagerA-lr-c'+str(i)+'.log 2>&1')
    os.system('make_mask.py -p 6 -i 4 -m '+imagename+'.mask '+imagename+'.restored.corr > make_mask-lr-c'+str(i)+'.log  2>&1')
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager-lr.sh concat.MS '+imagename+'-masked '+imagename+'.mask > awimagerB-lr-c'+str(i)+'.log 2>&1')

    # make new low resolution models
    logging.info('Creating low resolution model...')
    img = bdsm.process_image(imagename+'-masked.restored.corr', mean_map='zero', adaptive_rms_box=True, thresh_isl=4, thresh_pix=6,\
            rms_box_bright=(20, 7), rms_box=(120, 40), atrous_do=False, ini_method='curvature', advanced_opts=True, blank_limit=1e-5)
    img.export_image(outfile=imagename+'-masked.BDSMresidual', clobber=True, img_type='gaus_resid')
    img.export_image(outfile=imagename+'-masked.BDSMmodel', clobber=True, img_type='gaus_model')
    img.write_catalog(outfile='widefield-lr_g'+g+'.skymodel', catalog_type='gaul', format='bbs', bbs_patches='source', clobber=True)

    # Concat models
    logging.info('Concatenating skymodels...')
    skymod = lsmtool.load('widefield_g'+g+'.skymodel')
    skymod.concatenate('widefield-lr_g'+g+'.skymodel')
    skymod.write('widefield-comb_g'+g+'.skymodel', clobber=True)

# Final subtract of the best model - group*_TC.MS:DATA -> group*_TC.MS:SUBTRACTED_DATA (beam corrupted - all sources subtracted)
cmds = []
logging.info('Final subtraction...')
for ms in mss:
    cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/bbs-self_subfinal.parset widefield-comb_g'+g+'.skymodel > '+ms+'_final-sub.log 2>&1')
thread_cmd(cmds, max_threads=max_threads)

# TODO: add a very large CASA run and subtract model from CORRECTED_DATA (must be beam-corrected!), then beam corrupt the very empty dataset
# this beacause awimages is not able to make widewide fov images

logging.info("Done.")
