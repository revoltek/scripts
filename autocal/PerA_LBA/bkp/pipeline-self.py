#!/usr/bin/env python
# apply solution from the calibrator and then run selfcal, an initial model
# and initial solutions from a calibrator must be provided
# local dir must contain all the MSs, un touched in linear pol
# 1 selfcal
# 2 subtract field
# 3 flag

# initial self-cal model
model = '/home/fdg/scripts/autocal/VirA_LBA/150228_LBA-VirgoA.model'
# globaldb produced by pipeline-init
globaldb = '../cals-lin/globaldb'
# fake skymodel with pointing direction
fakeskymodel = '/home/fdg/scripts/autocal/VirA_LBA/virgo.fakemodel.skymodel'
# SB per block
n = 20

##############################################################

import sys, os, glob, re
from lofar import bdsm
import numpy as np
import lsmtool
import pyrap.tables as pt
from lib_pipeline import *

set_logger()

#################################################
# Clear
logging.info('Cleaning...')
check_rm('*log')
check_rm('*last')
check_rm('*h5')
check_rm('concat*MS')
check_rm('block*MS')
check_rm('img')
os.system('mkdir img')

# all MS
mss = sorted(glob.glob('*.MS'))
Nblocks = len(mss)/n
logging.debug("Number of blocks: "+str(Nblocks))
logging.debug("Blocks:")
for j, mss_block in enumerate(np.array_split(mss, Nblocks)):
    logging.debug(str(j)+": "+str(mss_block)+" - len: "+str(len(mss_block)))

#################################################
# Copy cal solution
for ms in sorted(mss):
    num = re.findall(r'\d+', ms)[-1]
    check_rm(ms+'/instrument')
    print 'cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument'
    os.system('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')

#########################################################################################
# [PARALLEL] Correct beam - SB.MS:DATA -> SB.MS:DATA_BC (data, beam applied, linear)
#logging.info('Apply beam...')
#cmds = []
#for ms in mss:
#    append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/VirA_LBA/parsets_self/bbs-beamcorr.parset /home/fdg/scripts/autocal/VirA_LBA/virgo.fakemodel.skymodel > '+ms+'-init_beamcorr.log 2>&1')
#thread_cmd(cmds)
#logging.warning('Bad runs:')
#os.system('grep -L "successfully" *-init_beamcorr.log')
   
#########################################################################################
# [PARALLEL] apply solutions and beam correction - SB.MS:DATA -> SB.MS:CALCOR_DATA (calibrator corrected data, beam applied, circular)
logging.info('Correcting target MSs...')
cmds=[]
for ms in mss:
    cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/VirA_LBA/parsets_self/bbs-corbeam.parset '+fakeskymodel+' > '+ms+'-init_corbeam.log 2>&1')
#thread_cmd(cmds)
logging.warning('Bad runs:')
os.system('grep -L "successfully" *-init_corbeam.log')

#########################################################################################
# [PARALLEL] Transform to circular pol - SB.MS:CALCOR_DATA -> SB-circ.MS:CIRC_DATA (data, beam applied, circular)
logging.info('Convert to circular...')
cmds = []
for ms in mss:
    cmds.append('/home/fdg/scripts/mslin2circ.py -i '+ms+':CALCOR_DATA -o '+ms+':CIRC_DATA > '+ms+'-init_circ2lin.log 2>&1')
#thread_cmd(cmds)
 
# self-cal cycle
for i in xrange(5):
    logging.info('Starting self-cal cycle: '+str(i))

    # MS for calibration, use all at cycle 0 and 4
    if i == 0 or i == 4:
        mss_c = mss
        mss_clean = mss[::n]
    else:
        mss_c = mss[::n]
        mss_clean = mss[::n]

    check_rm('concat*.MS')

    #####################################################################################
    # ft model (not parallel, problem accessing .model multiple times with casa) - model is unpolarized CIRC == LIN
    logging.info('Add models...')
    if i != 0:    
        model = 'img/normal_clean-iter'+str(i-1)+'.model'
    for ms in mss_c:
        os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_ft.py '+ms+' '+model+' > '+ms+'_ft-c'+str(i)+'.log 2>&1')

    #####################################################################################
    # [PARALLEL] calibrate - SB.MS:CIRC_DATA (no correction)
    # TODO: calibrated on wide-field subtracted data (impossible now because field can be subtracted only in LINEAR and calibration to corrupt model is available in CIRCULAR)
    logging.info('Calibrate...')
    cmds=[]
    for ms in mss_c:
        cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-selfcal_modeldata.parset msin='+ms+' cal.parmdb='+ms+'/instrument > '+ms+'_selfcal-c'+str(i)+'.log 2>&1')
    thread_cmd(cmds)
    logging.warning('Bad runs:')
    os.system('grep -L "Total NDPPP time" *_selfcal_c'+str(i)+'.log')

    #######################################################################################
    # Solution rescaling
    logging.info('Running LoSoTo to normalize solutions...')
    for ms in sorted(mss_c):
       h5parm = ms.replace('.MS','.h5')
       check_rm(h5parm)
       os.system('H5parm_importer.py -v '+h5parm+' '+ms+' > '+ms+'_losoto-c'+str(i)+'.log 2>&1')
       os.system('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/VirA_LBA/parsets_self/losoto.parset >> '+ms+'_losoto-c'+str(i)+'.log 2>&1')
       os.system('H5parm_exporter.py -v -c '+h5parm+' '+ms+' >> '+ms+'_losoto-c'+str(i)+'.log 2>&1')
       check_rm(ms+'/instrument')
       os.system('mv '+ms+'/sol000_instrument'+' '+ms+'/instrument')
    
    ########################################################################################
    # [PARALLEL] correct - SB.MS:CIRC_DATA -> SB.MS:CIRC_CORRECTED_DATA (selfcal corrected data, beam applied, circular)
    logging.info('Correct...')
    cmds=[]
    for ms in mss_c:
        cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument > '+ms+'_selfcor-c'+str(i)+'.log 2>&1')
    thread_cmd(cmds)
    logging.warning('Bad runs:')
    os.system('grep -L "Total NDPPP time" *_selfcor-c'+str(i)+'.log')

#######################
#   QUICK TEST LOOP
#    # avg - SB.MS:CIRC_CORRECTED_DATA -> concat-avg.MS:DATA
    logging.info('Average...')
    os.system('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-concatavg.parset msin="['+','.join(mss_clean)+']" msout=concat-avg.MS \
            msin.datacolumn=CIRC_CORRECTED_DATA > concatavg-c'+str(i)+'.log 2>&1')

    # clean (make a new model of virgo)
    logging.info('Clean...')
    os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_clean.py concat-avg.MS iter'+str(i)+' > clean-c'+str(i)+'.log 2>&1')

    continue
#
######################

    # [PARALLEL] Transform back to linear pol - SB.MS:CIRC_CORRECTED_DATA -> SB.MS:CORRECTED_DATA (selfcal corrected data, beam applied, linear)
    logging.info('Convert back to linear...')
    cmds = []
    for ms in mss_c:
        cmds.append('/home/fdg/scripts/mslin2circ.py -r -i '+ms+':CIRC_CORRECTED_DATA -o '+ms+':CORRECTED_DATA > '+ms+'_circ2lin-c'+str(i)+'.log 2>&1')
    thread_cmd(cmds)

    # create widefield model
    if i == 0 or i == 4:
        check_rm('block_*MS')

        # concatenate data - MS.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, linear)
        logging.info('Make widefield model - Concatenating...')
        cmds = []
        for j, mss_block in enumerate(np.array_split(mss, Nblocks)):
            cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-concat.parset msin="['+','.join(mss_block)+']" msout=block_'+str(j)+'.MS > wide_concat-b'+str(j)+'.log 2>&1')
        thread_cmd(cmds)
        logging.warning('Bad runs:')
        os.system('grep -L "Total NDPPP time" wide_concat-b*.log')

        # subtract virgo - concat.MS:DATA -> concat.MS:CORRECTED_DATA (selfcal corrected data virgo subtracted, beam applied, linear)
        logging.info('Make widefield model - create columns...')
        cmds = []
        for j in xrange(Nblocks):
            cmds.append('awimager ms=block_'+str(j)+'.MS image=/nonexistant > wide_fakeawimager-b'+str(j)+'.log 2>&1') # creates MODEL_DATA and CORRECTED_DATA in the way awimages likes
        thread_cmd(cmds)
        # ft, not parallel!
        logging.info('Make widefield model - ft() Virgo A...')
        for j in xrange(Nblocks):
            os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_ft.py block_'+str(j)+'.MS '+model+' > wide_ft-b'+str(j)+'.log 2>&1')
        # uvsub not parallel!
        logging.info('Make widefield model - UV-Subtracting Virgo A...')
        for j in xrange(Nblocks):
            os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_uvsub.py block_'+str(j)+'.MS > wide_uvsub-b'+str(j)+'.log 2>&1')

        # clean, mask, clean
        logging.info('Make widefield model - Widefield imaging...')
        for j in xrange(Nblocks):
            os.system('/home/fdg/scripts/autocal/VirA_LBA/awimager.sh block_'+str(j)+'.MS widefield_b'+str(j)+'-iter'+str(i)+' > wide_awimager1-b'+str(j)+'.log 2>&1')
            os.system('/home/fdg/scripts/autocal/VirA_LBA/make_mask.py -t -m img/widefield_b'+str(j)+'-iter'+str(i)+'.mask img/widefield_b'+str(j)+'-iter'+str(i)+'.restored.corr > wide_makemask-b'+str(j)+'.log 2>&1')
            os.system('/home/fdg/scripts/autocal/VirA_LBA/awimager.sh block_'+str(j)+'.MS widefieldM_b'+str(j)+'-iter'+str(i)+' img/widefield_b'+str(j)+'-iter'+str(i)+'.mask  > wide_awimager2-b'+str(j)+'.log 2>&1')

        logging.info('Make widefield model - Creating model...')
        for j in xrange(Nblocks):
            img = bdsm.process_image('img/widefieldM_b'+str(j)+'-iter'+str(i)+'.restored.corr', mean_map='zero', adaptive_rms_box=True, \
                    rms_box_bright=(20, 7), rms_box=(120, 40), atrous_do=False, ini_method='curvature', advanced_opts=True, blank_limit=1e-5)
            img.export_image(outfile='img/widefieldM_b'+str(j)+'-iter'+str(i)+'.BDSMresidual', clobber=True, img_type='gaus_resid')
            img.export_image(outfile='img/widefieldM_b'+str(j)+'-iter'+str(i)+'.BDSMmodel', clobber=True, img_type='gaus_model')
            img.write_catalog(outfile='widefield_b'+str(j)+'.skymodel', catalog_type='gaul', format='bbs', bbs_patches='source', clobber=True)

            # LSM tool remove central part (m87) from source to subtract in case of false detections
            skymod = lsmtool.load('widefield_b'+str(j)+'.skymodel')
            dist = skymod.getDistance(187.7059304, 12.3911231) # m87 coords
            skymod.select(dist > 0.16) # keep only sources outside central 1/4 deg
            skymod.write('widefield_b'+str(j)+'.skymodel', clobber=True)

    # [PARALLEL] beam corrupt - concat.MS:DATA -> concat.MS:DATA_NOBEAM (selfcal corrected data, beam corrupted, linear)
    logging.info('Corrupt for the beam...')
    cmds = []
    for ms in mss_c:
        cmds.append('/home/fdg/scripts/autocal/VirA_LBA/runTAMMOenv.sh NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-beamcorrupt.parset msin='+ms+' > '+ms+'_beamcorrupt-c'+str(i)+'.log 2>&1')
    thread_cmd(cmds)

    # [PARALLEL] subtract skymodel sb per sb and correct for beam - concat.MS:DATA_NOBEAM -> concat.MS:CORRECTED_DATA (selfcal corrected data, beam applied, linear, source subtracted)
    logging.info('Subtracting skymodel and beam correct...')
    cmds = []
    for j, ms in enumerate(mss_c):
        print 'CKECK: sub widefield_b'+str(j)+'.skymodel from '+ms
        cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/VirA_LBA/parsets_self/bbs-subfield.parset widefield_b'+str(j)+'.skymodel > '+ms+'_subbeamcorr-c'+str(i)+'.log 2>&1')
    thread_cmd(cmds)

    # Flag on the residual
    if i==0 or i == 4:

        # concat all SBs - SB.MS:CORRECTED_DATA -> concat_all.MS:DATA (selfcal corrected data, beam applied, linear)
        logging.info('Concat...')
        os.system('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-concat.parset msin="['+','.join(mss)+']" msout=concat_all.MS > flag_concat.log 2>&1')
        logging.info('Flagging residuals...')

        os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBAiparsets_self/casa_ft.py concat_all.MS '+model+' > flag_ft.log 2>&1')
        os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_uvsub.py concat_all.MS > flag_uvsub.log 2>&1')
        os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_flag.py concat_all.MS > flag_flag.log 2>&1')

        logging.info('Copy flags...')
        tc = pt.table('concat_all.MS')
        flag = tc.getcol('FLAG')

        logging.info('Flagging: '+str(100.*np.count_nonzero(flag)/flag.size)+'%')

        # copy on the SB.MS
        n_chan_per_sb = len(flag[1])/len(mss)
        print "Found "+str(n_chan_per_sb)+" channel per SB."
        for j, ms in enumerate(mss):
            logging.debug('Copy flags into '+ms)
            tc2 = pt.table(ms, readonly=False)
            tc2.putcol('FLAG', flag[:,n_chan_per_sb*j:n_chan_per_sb*j+n_chan_per_sb,:])
            tc2.close()
        tc.close()

    # avg 1chanSB/20s - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, linear)
    logging.info('Average...')
    os.system('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-concatavg.parset msin="['+','.join(mss_clean)+']" msout=concat.MS > concatavg-c'+str(i)+'.log 2>&1')

    # clean (make a new model of virgo)
    logging.info('Clean...')
    os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_clean.py concat.MS iter'+str(i)+' > clean-c'+str(i)+'.log 2>&1')

##########################################################################################################
# [PARALLEL] concat+avg - SB.MS:CORRECTED_DATA -> concat.MS:DATA (selfcal corrected data, beam applied, linear)
logging.info('Concat...')
cmds = []
for j, mss_block in enumerate(np.array_split(mss, Nblocks)):
    cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-concatavg.parset msin="['+','.join(mss_block)+']" msout=concat-avg_b'+str(j)+'.MS > final_concatavg-b'+str(j)+'.log 2>&1')
    # fast version
    #cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-concatavg.parset msin="['+','.join(mss_block)+']" msout=concat-avg_b'+str(j)+'.MS msin.datacolumn=CIRC_CORRECTED_DATA > final_concatavg-b'+str(j)+'.log 2>&1')
thread_cmd(cmds)
# avg using the block, otherwise NDPPP crashes for too many file open
os.system('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-concatavg.parset msin=concat-avg_b*MS msout=concat-avg_all.MS avg.freqstep=1 > final_concatavg-all.log 2>&1')
# fast version
#os.system('NDPPP /home/fdg/scripts/autocal/VirA_LBA/parsets_self/NDPPP-concatavg.parset msin=concat-avg_b*MS msout=concat-avg_all.MS avg.freqstep=1 msin.datacolumn=DATA > final_concatavg-all.log 2>&1')

#########################################################################################################
# group images of VirA
logging.info('Block images...')
os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_clean.py concat-avg_all.MS all > final_clean-all.log 2>&1')
cmds = []
for j in xrange(Nblocks):
    cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_clean.py concat-avg_b'+str(j)+'.MS block'+str(j)+' > final_clean-b'+str(j)+'.log 2>&1')
thread_cmd(cmds, max_threads=3)

#########################################################################################################
# low-res image (to be done at the end -> uvsub)
logging.info('Make low-resolution image...')
os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/parsets_self/casa_clean_lr.py concat-avg_all.MS all > final_lrclean_all.log 2>&1')

logging.info("Done.")
