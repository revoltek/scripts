#!/usr/bin/python
# perform self-calibration on a group of SBs concatenated in TCs. Script must be run in dir with MS.
# number/chan in MS are flexible but the must be concatenable (same chans/freq!)
# TCs should have calibrator corrected (a+p) data in DATA (beam not applied).
# file format of TCs is: group#_TC###.MS.

skymodel = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/toothbrush.GMRT150.skymodel'
fake_skymodel = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/toothbrush.fakemodel.skymodel'

#######################################################################################

import sys, os, glob, re
from Queue import Queue
from lofar import bdsm
import pyrap.tables as pt
import numpy as np
import lsmtool
from pipeline_lib import *

if os.path.exists('img'): os.system('rm -r img')
os.system('rm -r concat*MS *log *skymodel *last')
os.system('mkdir img')

set_logger()

# create a fake parmdb to be used later for merging slow-amp and fast-phase parmdbs
logging.info('Creating fake parmdb...')
cmds = []
for ms in glob.glob('group*_TC*.MS'):
    cmds.append('calibrate-stand-alone -f --parmdb-name instrument_empty '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/bbs-fakeparmdb.parset '+fake_skymodel+' > '+ms+'-fakeparmdb.log 2>&1')
thread_cmd(cmds)

# TODO: a run of AOflagger on concatenated data

# self-cal cycle
for i in xrange(4):
    logging.info('Start selfcal cycle: '+str(i))

    if i != 0: skymodel = 'widefield-comb.skymodel'

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
        cmds = []
        # calibrate phase-only - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase corrected, beam corrected)
        logging.info('Calibrating phase...')
        for ms in glob.glob('group*_TC*.MS'):
            cmds.append('calibrate-stand-alone -f '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/bbs-self_solcor.parset '+skymodel+' > '+ms+'-cal-'+str(i)+'.log 2>&1')
        thread_cmd(cmds)
    else:
        cmds = []
        # calibrate phase-only - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA_PHASE (selfcal phase corrected, beam corrupted)
        logging.info('Calibrating phase...')
        for ms in glob.glob('group*_TC*.MS'):
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_csp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/bbs-self_solcor_csp.parset '+skymodel+' > '+ms+'-calpreamp-'+str(i)+'.log 2>&1')
        thread_cmd(cmds)

        cmds = []
        # calibrate amplitude (only solve) - group*_TC.MS:CORRECTED_DATA_PHASE
        logging.info('Calibrating amplitude...')
        for ms in glob.glob('group*_TC*.MS'):
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_amp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/bbs-self_sol_amp.parset '+skymodel+' > '+ms+'-calamp-'+str(i)+'.log 2>&1')
        thread_cmd(cmds)

        # merge parmdbs
        logging.info('Merging instrument tables...')
        for ms in glob.glob('group*_TC*.MS'):
            if os.path.exists(ms+'/instrument'): os.system('rm -r '+ms+'/instrument')
            merge_parmdb(ms+'/instrument_amp', ms+'/instrument_csp', ms+'/instrument_empty', ms+'/instrument')

        # TODO: rescaling
        #for ms in glob.glob('concat-TC*.MS'):
            #if os.path.exists('tgt.h5'): os.system('rm tgt.h5')
            #os.system('H5parm_importer.py -v tgt.h5 '+ms)
            #os.system('losoto.py -v tgt.h5 /home/fdg/scripts/autocal/1RXSJ0603_LBA/losoto-tgt.parset')
            #os.system('H5parm_exporter.py -v tgt.h5 '+ms)
            #   os.system('rm -r '+ms+'/instrument')
            #   os.system('mv '+ms+'/instrument-sol000'+' '+ms+'/instrument')

        cmds = []
        # correct amplitude - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase+amp corrected, beam corrected)
        logging.info('Correcting amplitude...')
        for ms in glob.glob('group*_TC*.MS'):
            cmds.append('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/bbs-self_cor_ampcsp.parset '+skymodel+' > '+ms+'-corampcsp-'+str(i)+'.log 2>&1')
        thread_cmd(cmds)

    # TEST for circular
    # join RR and LL
    #if os.path.exists(ms): os.system('rm -r '+ms)
    #os.system('mv '+msRR+' '+ms)
    #print 'taql \'update '+ms+', '+msLL+' as ll set DATA[3,]=ll.DATA[3,]\''
    #os.system('taql \'update '+ms+', '+msLL+' as ll set DATA[3,]=ll.DATA[3,]\'')
    #os.system('rm -r '+msLL)

    # concat all TCs in one MS - concat_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected)
    # NOTE: observations must have the same channels i.e. same SBs freq!
    logging.info('Concatenating TCs...')
    if os.path.exists('concat.MS'): os.system('rm -r concat.MS')
    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/casa_concat.py group*_TC*.MS concat.MS > split1-'+str(i)+'.log 2>&1')

    # clean mask clean TODO: try casa-clean with nterms=2 when using 30 SBs (anyway we are dominated by DDE) - try tapering to lower resolution
    logging.info('Cleaning...')
    imagename = 'img/wide-'+str(i)
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager.sh concat.MS '+imagename+' > awimager1-'+str(i)+'.log 2>&1')
    os.system('make_mask.py -p 6 -i 4 -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/tooth.mask -m '+imagename+'.mask '+imagename+'.restored.corr > make_mask-'+str(i)+'.log  2>&1')
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager.sh concat.MS '+imagename+'-masked '+imagename+'.mask > awimager2-'+str(i)+'.log 2>&1')

    # making new model
    logging.info('Creating model...')
    img = bdsm.process_image(imagename+'-masked.restored.corr', mean_map='zero', adaptive_rms_box=True, thresh_isl=4, thresh_pix=6,\
            rms_box_bright=(20, 7), rms_box=(120, 40), atrous_do=False, ini_method='curvature', advanced_opts=True, blank_limit=1e-5)
    img.export_image(outfile=imagename+'-masked.BDSMresidual', clobber=True, img_type='gaus_resid')
    img.export_image(outfile=imagename+'-masked.BDSMmodel', clobber=True, img_type='gaus_model')
    img.write_catalog(outfile='widefield.skymodel', catalog_type='gaul', format='bbs', bbs_patches='source', clobber=True)

    ####################################################################
    # DEBUG TEST
    #os.system('cp widefield.skymodel widefield-comb.skymodel')
    #continue
    ####################################################################

    # Subtract model from all TCs - concat.MS:DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
    cmds = []
    logging.info('Subtracting high-res model...')
    for ms in glob.glob('group*_TC*.MS'):
        if i < 2:
            cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/bbs-self_sub_csp.parset widefield.skymodel > '+ms+'-sub-'+str(i)+'.log 2>&1')
        else:
            cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/bbs-self_sub_ampcsp.parset widefield.skymodel > '+ms+'-sub-'+str(i)+'.log 2>&1')
    thread_cmd(cmds)

    # concat all TCs in one MS - concat_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
    logging.info('Concatenating TCs...')
    if os.path.exists('concat.MS'): os.system('rm -r concat.MS')
    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/casa_concat.py group*_TC*.MS concat.MS > split2-'+str(i)+'.log 2>&1')

    # reclean low-resolution
    logging.info('Cleaning low resolution...')
    imagename = 'img/wide-lr-'+str(i)
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager-lr.sh concat.MS '+imagename+' > awimager1-lr-'+str(i)+'.log 2>&1')
    os.system('make_mask.py -p 6 -i 4 -m '+imagename+'.mask '+imagename+'.restored.corr > make_mask-lr-'+str(i)+'.log  2>&1')
    os.system('/home/fdg/scripts/autocal/1RXSJ0603_LBA/awimager-lr.sh concat.MS '+imagename+'-masked '+imagename+'.mask > awimager2-lr-'+str(i)+'.log 2>&1')

    # make new low resolution models
    logging.info('Creating low resolution model...')
    img = bdsm.process_image(imagename+'-masked.restored.corr', mean_map='zero', adaptive_rms_box=True, thresh_isl=4, thresh_pix=6,\
            rms_box_bright=(20, 7), rms_box=(120, 40), atrous_do=False, ini_method='curvature', advanced_opts=True, blank_limit=1e-5)
    img.export_image(outfile=imagename+'-masked.BDSMresidual', clobber=True, img_type='gaus_resid')
    img.export_image(outfile=imagename+'-masked.BDSMmodel', clobber=True, img_type='gaus_model')
    img.write_catalog(outfile='widefield-lr.skymodel', catalog_type='gaul', format='bbs', bbs_patches='source', clobber=True)

    # Concat models
    logging.info('Concatenating skymodels...')
    skymod = lsmtool.load('widefield.skymodel')
    skymod.concatenate('widefield-lr.skymodel')
    skymod.write('widefield-comb.skymodel', clobber=True)

# Final subtract of the best model - concat_TC.MS:DATA -> concat_TC.MS:CORRECTED_DATA (beam corrupted - all sources subtracted)
cmds = []
logging.info('Final subtraction...')
for ms in glob.glob('group*_TC*.MS'):
    cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/bbs-self_subfinal.parset widefield-comb.skymodel > '+ms+'-sub-final.log 2>&1')
thread_cmd(cmds)

# TODO: add a very large (does it make sense in LBA with iono so strong?) CASA run and subtract model from CORRECTED_DATA (must be beam-corrected), then beam corrupt the very empty dataset

logging.info("Done.")
