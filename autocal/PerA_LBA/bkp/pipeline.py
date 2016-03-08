#!/usr/bin/env python

model = 'normal_clean-group0-iter3.model'
ngroups = 1 # TODO: put to 1 for testing

##############################################################

import sys, os, glob, re
from Queue import Queue
from threading import Thread
import subprocess
from lofar import bdsm
import logging
import numpy as np
import lsmtool
import pyrap.tables as pt
from pipeline_lib import *

#################################################
# Clean
logging.info('Cleaning...')
os.system('rm -r img/* *log cals/*log tgts/*log')

set_logger()

################################################
# Integrity check
for cal in sorted(glob.glob('cals/*MS')):
    if not os.path.exists(cal.replace('cals','tgts').replace('3C295','VirgoA')):
        logging.warning("Cannot find " + cal.replace('cals','tgts').replace('3C295','VirgoA')+' remove '+cal)
for tgt in sorted(glob.glob('tgts/*MS')):
    if not os.path.exists(tgt.replace('tgts','cals').replace('VirgoA','3C295')):
        logging.warning("Cannot find " + tgt.replace('tgts','cals').replace('VirgoA','3C295')+' remove '+tgt)

################################################
# [PARALLEL] initial calibration on calibrator
logging.info('Calibrating the calibrators...')
cmds=[]
for i, ms in enumerate(sorted(glob.glob('cals/*MS'))):
    if i < 40:
        cmds.append('calibrate-stand-alone -f '+ms+'  /home/fdg/scripts/autocal/VirA_LBA/bbs-cal_field-low.parset /home/fdg/model/3C295-allfield.skymodel > '+ms+'.log 2>&1')
    else:
        cmds.append('calibrate-stand-alone -f '+ms+'  /home/fdg/scripts/autocal/VirA_LBA/bbs-cal_field.parset /home/fdg/model/3C295-allfield.skymodel > '+ms+'.log 2>&1')
#thread_cmd(cmds)
#logging.warning('Bad runs:')
#os.system('grep -L success cals/*log')

##############################################
# Clock check and flagging and sol transfer
#if not os.path.exists('cals/globaldb'):
#    os.system('mkdir cals/globaldb')
#
#logging.info('Running LoSoTo...')
#for ms in sorted(glob.glob('cals/*MS')):
#    logging.debug('Copy instrument of '+ms)
#    num = re.findall(r'\d+', ms)[-1]
#    os.system('cp -r '+ms+'/instrument cals/globaldb/instrument-'+str(num))
#    print 'cp -r '+ms+'/instrument cals/globaldb/instrument-'+str(num)
#    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky cals/globaldb/')
#
#if os.path.exists('cal.h5'): os.system('rm cal.h5')
#if os.path.exists('plot'): os.system('rm plot')
#os.system('mkdir plot')
#os.system('H5parm_importer.py -v cal.h5 cals/globaldb')
#os.system('losoto.py -v cal.h5 /home/fdg/scripts/autocal/VirA_LBA/losoto.parset')
#os.system('H5parm_exporter.py -v cal.h5 cals/globaldb')
#for ms in sorted(glob.glob('tgts/*MS')):
#    num = re.findall(r'\d+', ms)[-1]
#    if os.path.exists(ms+'/instrument'):
#        os.system('rm -r '+ms+'/instrument')
#    os.system('cp -r cals/globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
#    print 'cp -r cals/globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument'
    
#########################################################################################
# [PARALLEL] NDPPP applys only gain, not directionalgain, use BBS - SB.MS:DATA -> SB.MS:CORRECTED_DATA (calibrator corrected data, beam applied, linear)
logging.info('Correcting target MSs...')
cmds=[]
for i, ms in enumerate(sorted(glob.glob('tgts/*MS'))):
    cmds.append('calibrate-stand-alone --replace-sourcedb '+ms+' /home/fdg/scripts/autocal/VirA_LBA/bbs-cal_correct.parset /home/fdg/scripts/autocal/VirA_LBA/virgo.fakemodel.skymodel > '+ms+'-apply.log 2>&1')
#thread_cmd(cmds)
#logging.warning('Bad runs:')
#os.system('grep -L success tgts/*log')

###########################################################################################
# [PARALLEL] Transform to circular pol - SB.MS:CORRECTED_DATA -> SB-circ.MS:DATA (calibrator corrected data, beam applied, circular)
logging.info('Convert to circular...')
cmds = []
for ms in glob.glob('tgts/*MS'):
    cmds.append('/home/fdg/scripts/mslin2circ.py -i '+ms+':CORRECTED_DATA -o '+ms.replace('.MS','-circ.MS')+':DATA > '+ms+'-circ2lin.log 2>&1')
#thread_cmd(cmds)

# self-cal cycle
logging.info('Starting self-cal cycles...')
for i in xrange(5):

    # ft model (not parallel, problem accessing .model multiple times with casa) - model is unpolarized CIRC == LIN
    logging.info('Add models...')
    for g, group_mss in enumerate(np.array_split(sorted(glob.glob('tgts/*-circ.MS')), ngroups)):
        if i != 0:    
            model = 'img/normal_clean-group'+str(g)+'-iter'+str(i-1)+'.model'
        for ms in group_mss:
            os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_ft.py '+ms+' '+model+' 2 > '+ms+'-ft.log 2>&1')

    # [PARALLEL] calibrate TODO: calibrated on wide-field subtracted data
    logging.info('Calibrate...')
    cmds=[]
    for j, ms in enumerate(sorted(glob.glob('tgts/*-circ.MS'))):
        cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/NDPPP-selfcal_modeldata.parset msin='+ms+' cal.parmdb='+ms+'/instrument > '+ms+'-cal.log 2>&1')
    thread_cmd(cmds)
    logging.warning('Bad runs:')
    os.system('grep -L "Total NDPPP time" tgts/*-cal.log')

    # Solution flagging TODO: debug losoto
    #logging.info('Running LoSoTo to normalize solutions...')
    #for ms in sorted(glob.glob('tgts/*-circ.MS')):
    #   h5parm = ms.replace('.MS','.h5')
    #   if os.path.exists(h5parm): os.system('rm '+h5parm)
    #   os.system('H5parm_importer.py -v '+h5parm+' '+ms)
    #   os.system('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/VirA_LBA/losoto-tgt.parset')
    #   os.system('H5parm_exporter.py -v '+h5parm+' '+ms)
    #   os.system('rm -r '+ms+'/instrument')
    #   os.system('mv '+ms+'/instrument-sol000'+' '+ms+'/instrument')
    
    # [PARALLEL] correct - SB-circ.MS:DATA -> SB-circ.MS:CORRECTED_DATA (selfcal corrected data, beam applied, circular)
    logging.info('Correct...')
    cmds=[]
    for j, ms in enumerate(sorted(glob.glob('tgts/*-circ.MS'))):
        cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/NDPPP-selfcor.parset msin='+ms+' cor.parmdb='+ms+'/instrument > '+ms+'-corr.log 2>&1')
    thread_cmd(cmds)
    logging.warning('Bad runs:')
    os.system('grep -L "Total NDPPP time" tgts/*-corr.log')

#######################
#   QUICK TEST LOOP
#    # avg - concat*-circ.MS:CORRECTED_DATA -> concat-avg.MS:DATA
#    os.system('rm -r tgts/concat*-avg.MS')
#    logging.info('Average...')
#    cmds=[]
#    for g, group_mss in enumerate(np.array_split(sorted(glob.glob('tgts/*-circ.MS')), ngroups)):
#        cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/NDPPP-avg.parset msin="['+','.join(group_mss)+']" msout=tgts/concat-'+str(g)+'-avg.MS > tgts/concat-'+str(g)+'.MS-concat.log 2>&1')
#    thread_cmd(cmds)

    # clean (make a new model of virgo)
#    logging.info('Clean...')
#    cmds=[]
#    for g, concat_ms in enumerate(glob.glob('tgts/concat*-avg.MS')):
#        cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_clean.py '+concat_ms+' group'+str(g)+'-iter'+str(i)+' > '+concat_ms+'-clean.log 2>&1')
#    thread_cmd(cmds, max_threads=3)

#    continue
######################

    # [PARALLEL] Transform back to linear pol - SB-circ.MS:CORRECTED_DATA -> SB-lin.MS:DATA (selfcal corrected data, beam applied, linear)
    logging.info('Convert back to linear...')
    os.system('rm -r tgts/*-lin.MS')
    cmds = []
    for ms in glob.glob('tgts/*-circ.MS'):
        cmds.append('/home/fdg/scripts/mslin2circ.py -r -i '+ms+':CORRECTED_DATA -o '+ms.replace('-circ.MS','-lin.MS')+':DATA > '+ms+'-circ2lin.log 2>&1')
    thread_cmd(cmds)

    # [PARALLEL] concat - SB-lin.MS:DATA -> concat.MS:DATA (selfcal corrected data, beam applied, linear)
    logging.info('Concat...')
    os.system('rm -r tgts/concat*MS')
    cmds=[]
    for g, group_mss in enumerate(np.array_split(sorted(glob.glob('tgts/*-lin.MS')), ngroups)):
        cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/NDPPP-concat.parset msin="['+','.join(group_mss)+']" msout=tgts/concat-'+str(g)+'.MS > tgts/concat-'+str(g)+'.MS-concat.log 2>&1')
    thread_cmd(cmds)

##########################
#    # DEBUG 2 clean (make a new model of virgo)
#    logging.info('Clean...')
#    cmds=[]
#    for g, concat_ms in enumerate(glob.glob('tgts/concat*.MS')):
#        cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_clean.py '+concat_ms+' debug2-group'+str(g)+'-iter'+str(i)+' > '+concat_ms+'-clean.log 2>&1')
#    thread_cmd(cmds, max_threads=3)
#########################

    # subtraction of other sources - concat.MS:DATA -> concat.MS:CORRECTED_DATA (selfcal corrected data, beam applied, linear, virgo subtracted)
    if i == 0 or i == 4:
        for g, concat_ms in enumerate(glob.glob('tgts/concat*.MS')):
            logging.info('Subtracting Virgo A...')
            os.system('awimager ms='+concat_ms+' image=/nonexistant > '+concat_ms+'-fake_awimager.log 2>&1') # creates MODEL_DATA and CORRECTED_DATA in the way awimages likes
            os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_ft.py '+concat_ms+' '+model+' 2 > '+concat_ms+'-ft.log 2>&1')
            os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_uvsub.py '+concat_ms+' > '+concat_ms+'-uvsub.log 2>&1')

            # clean, make mask, reclean again
            logging.info('Widefield imaging...')
            os.system('/home/fdg/scripts/autocal/VirA_LBA/awimager.sh '+concat_ms+' widefield-group'+str(g)+'-iter'+str(i)+' > '+concat_ms+'-awimager1.log 2>&1')
            os.system('/home/fdg/scripts/autocal/VirA_LBA/make_mask.py -t -m img/widefield-group'+str(g)+'-iter'+str(i)+'.mask img/widefield-group'+str(g)+'-iter'+str(i)+'.restored.corr > '+concat_ms+'-makemask.log 2>&1')
            os.system('/home/fdg/scripts/autocal/VirA_LBA/awimager.sh '+concat_ms+' widefieldM-group'+str(g)+'-iter'+str(i)+' img/widefield-group'+str(g)+'-iter'+str(i)+'.mask  > '+concat_ms+'-awimager2.log 2>&1')

            logging.info('Creating model...')
            img = bdsm.process_image('img/widefieldM-group'+str(g)+'-iter'+str(i)+'.restored.corr', mean_map='zero', adaptive_rms_box=True, \
                    rms_box_bright=(20, 7), rms_box=(120, 40), atrous_do=False, ini_method='curvature', advanced_opts=True, blank_limit=1e-5)
            img.export_image(outfile='img/widefieldM-group'+str(g)+'-iter'+str(i)+'.BDSMresidual', clobber=True, img_type='gaus_resid')
            img.export_image(outfile='img/widefieldM-group'+str(g)+'-iter'+str(i)+'.BDSMmodel', clobber=True, img_type='gaus_model')
            img.write_catalog(outfile='widefield-'+str(g)+'.skymodel', catalog_type='gaul', format='bbs', bbs_patches='source', clobber=True)

            logging.info('Removing possible surious sources...')
            # LSM tool remove central part (m87) from source to subtract in case of false detections
            skymod = lsmtool.load('widefield-'+str(g)+'.skymodel')
            dist = skymod.getDistance(187.7059304, 12.3911231) # m87 coords
            skymod.select(dist > 0.16) # keep only sources outside central 1/4 deg
            skymod.write('widefield-'+str(g)+'.skymodel', clobber=True)

    # [PARALLEL] beam corrupt - concat.MS:DATA -> concat.MS:DATA_NOBEAM (selfcal corrected data, beam corrupted, linear)
    logging.info('Corrupt for the beam...')
    cmds=[]
    for g, concat_ms in enumerate(glob.glob('tgts/concat*.MS')):
        cmds.append('/home/fdg/scripts/autocal/VirA_LBA/runTAMMOenv.sh NDPPP /home/fdg/scripts/autocal/VirA_LBA/NDPPP-beamcorrupt.parset msin='+concat_ms+' > '+concat_ms+'-beamcorrupt.log 2>&1')
    thread_cmd(cmds)

    # [PARALLEL] subtract skymodel and correct for beam - concat.MS:DATA_NOBEAM -> concat.MS:CORRECTED_DATA (selfcal corrected data, beam applied, linear, source subtracted)
    logging.info('Subtracting skymodel and beam correct...')
    cmds=[]
    for g, concat_ms in enumerate(glob.glob('tgts/concat*.MS')):
        cmds.append('calibrate-stand-alone --replace-sourcedb '+concat_ms+' /home/fdg/scripts/autocal/VirA_LBA/bbs-tgt_subfield.parset widefield-'+str(g)+'.skymodel > '+concat_ms+'-subeamcorr.log 2>&1')
    thread_cmd(cmds)

##################
#    # DEBUG 3 clean
#    logging.info('Clean...')
#    cmds=[]
#    for g, concat_ms in enumerate(glob.glob('tgts/concat*.MS')):
#        cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_clean.py '+concat_ms+' debug3-group'+str(g)+'-iter'+str(i)+' > '+concat_ms+'-clean.log 2>&1')
#    thread_cmd(cmds, max_threads=3)
##################

    # Flag on the residual
    if i==0 or i == 4:
        all_mss = np.array_split(sorted(glob.glob('tgts/*-circ.MS')), ngroups)
        for g, concat_ms in enumerate(glob.glob('tgts/concat*.MS')):
            logging.info('Flagging residuals...')
            concat_ms_flag = 'tgts/concat-'+str(g)+'-flag.MS' # need to create a new MS because we do a uvsub
            os.system('NDPPP /home/fdg/scripts/autocal/VirA_LBA/NDPPP-split.parset msin='+concat_ms+' msout='+concat_ms_flag+' > '+concat_ms+'-flagsplit.log 2>&1')
            os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_ft.py '+concat_ms_flag+' '+model+' 2 > '+concat_ms_flag+'-flagft.log 2>&1')
            os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_uvsub.py '+concat_ms_flag+' > '+concat_ms_flag+'-flaguvsub.log 2>&1')
            os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_flag.py '+concat_ms_flag+' > '+concat_ms_flag+'-flag.log 2>&1')

            logging.info('Copy flags...')
            tc = pt.table(concat_ms_flag)
            flag = tc.getcol('FLAG')

            # copy on the concat#.MS
            logging.debug('Copy flags into '+concat_ms)
            tc2 = pt.table(concat_ms, readonly=False)
            flag2 = tc2.getcol('FLAG')
            tc2.putcol('FLAG', flag)
            tc2.close()
            logging.info('Flagging: '+str(100.*np.count_nonzero(flag2)/flag2.size)+' -> '+str(100.*np.count_nonzero(flag)/flag.size))

            # copy on the -circ.MS
            n_chan_per_sb = len(flag[1])/len(all_mss[g])
            #print "Found "+str(n_chan_per_sb)+" channel per SB."
            for j, ms in enumerate(all_mss[g]): # these are the mss related with this group
                logging.debug('Copy flags into '+ms)
                tc2 = pt.table(ms, readonly=False)
                tc2.putcol('FLAG', flag[:,n_chan_per_sb*j:n_chan_per_sb*j+n_chan_per_sb,:])
                tc2.close()
            tc.close()
            os.system('rm -r '+concat_ms_flag)

    # [PARALLEL] avg - concat.MS:CORRECTED_DATA -> concat-avg.MS:DATA
    logging.info('Average...')
    cmds=[]
    for g, concat_ms in enumerate(glob.glob('tgts/concat*.MS')):
        cmds.append('NDPPP /home/fdg/scripts/autocal/VirA_LBA/NDPPP-avg.parset msin='+concat_ms+' msout=tgts/concat-'+str(g)+'-avg.MS > '+concat_ms+'-avg.log 2>&1')
    thread_cmd(cmds)

    # [PARALLEL] clean (make a new model of virgo)
    logging.info('Clean...')
    cmds=[]
    for g, concat_ms in enumerate(glob.glob('tgts/concat*-avg.MS')):
        cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_clean.py '+concat_ms+' group'+str(g)+'-iter'+str(i)+' > '+concat_ms+'-clean.log 2>&1')
    thread_cmd(cmds, max_threads=3)

# [PARALLEL] low-res image
logging.info('Make low-resolution image...')
cmds=[]
for g, concat_ms in enumerate(glob.glob('tgts/concat*-avg.MS')):
    cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/VirA_LBA/casa_clean_lr.py '+concat_ms+' group'+str(g)+'-lr > '+concat_ms+'-lrclean.log 2>&1')
thread_cmd(cmds, max_threads=3)

logging.info("Done.")
