#!/usr/bin/python
# perform peeling on a group of TCs. Script should be run in the directory with MSs inside.
# widefield-group#.skymodel is required for each group
# Input:
# group#_TC###.MS must have instrument table with amp+phase+TEC and beam uncorrected "empty" data in SUBTRACTED_DATA
# it should be possible to collect MSs after the run of pipeline-self.py
# .model images are in the "self/models" dir made by pipeline-self
# Output:
# DD-calibrated and subtracted data in SUBTRACTED_DATA
# Image of the facet

# coordinate of the region to find the DD
#coord = [90.833333,42.233333] # toorhbrush
#coord = [91.733333,41.680000] # strong pts
dd = {'name': 'src1', 'coord':[91.733333,41.680000], 'extended': False, 'facet_extended': False, 'reg': '', 'reg_facet': ''}
phasecentre = [90.833333,42.233333] # toorhbrush

##########################################################################################

import sys, os, glob, re
from lofar import bdsm
import pyrap.tables as pt
import numpy as np
import lsmtool
from lib_pipeline import *

set_logger()

#########################################################################################
# Clear
os.info('Cleaning...')
rm_check('peel*MS') 
rm_check('facet*MS') 
rm_check('*log')
rm_check('peel/'+dd['name'])
os.makedirs('peel/'+dd['name'])
os.makedirs('peel/'+dd['name']+'/models')
os.makedirs('peel/'+dd['name']+'/images')
os.makedirs('peel/'+dd['name']+'/instruments')

groups = []
tcs = []

for ms in sorted(glob.glob('group*_TC*.MS')):
    g = re.findall(r'\d+', ms)[0] # group number
    tc = re.findall(r'\d+', ms)[1] # time chunk number
    groups.append(g)
    tcs.append(tc)

#################################################################################################
# Blank unwanted part of models
logging.info('Splitting skymodels...')
for i, model in enumerate(sorted(glob.glob('self/models/wide-g*.model*'))):
    g = re.findall(r'\d+', skymodel)[0] # group number
    peelmodel = os.path.basename(model).replace('wide','peel')
    os.system('cp -r '+model+' peel/'+dd['name']+'/models/'+peelmodel)
    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_blank.py '+dd['reg']+' '+peelmodel)
    # redo the same but leaving the entire facet in the model
    peelmodel = os.path.basename(model).replace('wide','peel_facet')
    os.system('cp -r '+model+' peel/'+dd['name']+'/models/'+peelmodel)
    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_blank.py '+dd['reg_facet']+' '+peelmodel)

###############################################################################################################################
# [PARALLEL] Add DD cal model - group*_TC*.MS:MODEL_DATA (high+low resolution model)
logging.info('Add DD calibrator...')
cmds = []
for tc in tcs:
    # models are opened one at the time, no problems
    for ms in glob.glob('group*_TC'+tc+'.MS'):
        g = re.findall(r'\d+', ms)[0] # group number
        model = 'peel/'+dd['name']+'/models/peel-g'+g+'.model'
        modellr = 'peel/'+dd['name']+'/models/peel-lr-g'+g+'.model'
        cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_ft.py '+ms+' '+model+' > '+ms+'_init-ft.log 2>&1; \
            casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_ft.py '+ms+' '+modellr+' incremental >> '+ms+'_init-ft.log 2>&1')
    thread_cmd(cmds)

# [PARALLEL] ADD and corrupt group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + DD cal from model, cirular, beam correcred)
logging.info('Corrupting model...')
cmds = []
for ms in sorted(glob.glob('group*_TC*.MS')):
    cmds.append('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-init_add.parset > '+ms+'_init-add.log 2>&1')
thread_cmd(cmds)

# [PARALLEL] concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:CORRECTED_DATA -> peel-avg_TC*.MS:DATA (empty+DD, avg, phase shifted)
logging.info('Shifting+averaging...')
cmds = []
for tc in tcs:
    mss = glob.glob('group*_TC'+tc+'.MS')
    msout = 'peel-avg_TC'+tc+'.MS'
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' \
            shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\] > '+msout+'_init-shiftavg.log 2>&1')
thread_cmd(cmds)

# [PARALLEL] create a fake parmdb to be used later for merging slow-amp and fast-phase parmdbs
logging.info('Creating fake parmdb...')
cmds = []
for ms in glob.glob('peel-avg_TC*.MS'):
    cmds.append('calibrate-stand-alone -f --parmdb-name instrument_empty '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-fakeparmdb.parset > '+ms+'_init-fakeparmdb.log 2>&1')
thread_cmd(cmds)

###################################################################################################################
# self-cal cycle
for i in xrange(4):
    logging.info('Start peel cycle: '+str(i))
    check_rm('concat.MS')

    # [PARALLEL] averaging before cleaning peel-avg_TC*.MS:CORRECTED_DATA -> peel-avgavg_TC*.MS:DATA
    # TODO: check that in the first cycle the CORRECTED_DATA is ==DATA (where initial phaseshifted data are!)
    logging.info('Averaging before cleaning...')
    cmds = []
    for ms in glob.glob('peel-avg_TC*.MS'):
        msout = ms.reaplce('avg','avgavg')
        if i == 0: cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-avg.parset msin='+ms+' msin.datacolumn=DATA msout='+msout+' avg.freqstep=4 avg.timestep=10')
        else: cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-avg.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA msout='+msout+' avg.freqstep=4 avg.timestep=10')
    thread_cmd(cmds)

    # Concatenating (in time) before imaging peel-noavg_TC*.MS:DATA -> concat.MS:DATA (beam corrected, only source to peel in the data, all chan)
    # TODO: can we do virtual?
    logging.info('Concatenating TCs...')
    os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_concat.py peel-avgavg_TC*.MS concat-avg.MS > concat-c'+str(i)+'.log 2>&1')

    # Clean mask clean
    logging.info('Cleaning...')
    imagename = 'peel/'+dd['name']'/images/peel-'+str(i)
    if dd['extended']:
        os.system('casa --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_clean_peel-ext.py concat-avg.MS '+imagename+' > casaclean1-c'+str(i)+'.log 2>&1')
        os.system('make_mask.py -t -p 7 -i 5 -r 30,10 -m '+imagename+'.mask '+imagename+'.image.tt0 > make_mask-c'+str(i)+'.log  2>&1')
        os.system('casa --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_clean_peel-ext.py concat-avg.MS '+imagename+'-masked '+imagename+'.mask > casaclean2-c'+str(i)+'.log 2>&1')
    else: 
        os.system('casa --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_clean_peel.py concat-avg.MS '+imagename+' > casaclean1-c'+str(i)+'.log 2>&1')
        os.system('make_mask.py -p 7 -i 5 -r 30,10 -m '+imagename+'.mask '+imagename+'.image.tt0 > make_mask-c'+str(i)+'.log  2>&1')
        os.system('casa --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_clean_peel.py concat-avg.MS '+imagename+'-masked '+imagename+'.mask > casaclean2-c'+str(i)+'.log 2>&1')

    # ft model - peel-avg_TC*.MS:MODEL_DATA (best available model)
    logging.info('FT model...')
    for ms in glob.glob('peel-avg_TC*.MS'):
        os.system('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_ft.py '+ms+' '+imagename+'-masked.model > '+ms+'_ft-c'+str(i)+'.log 2>&1')

    if i < 2:
        ################################################################################################
        # [PARALLEL] calibrate phase-only - peel-avg_TC*.MS:DATA -> peel-avg_TC*.MS:CORRECTED_DATA
        logging.info('Calibrating+correcting phase...')
        cmds = []
        for ms in glob.glob('peel-avg_TC*.MS'):
            cmds.append('calibrate-stand-alone -f '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-solcor.parset peel-group'+groups[0]+'.skymodel > '+ms+'_solcor-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds)

    else:
        ##################################################################################################
        # [PARALLEL] calibrate phase-only - peel-avg_TC*.MS:DATA -> peel-avg_TC*.MS:CORRECTED_DATA_PHASE
        logging.info('Calibrating phase...')
        cmds = []
        for ms in glob.glob('peel-avg_TC*.MS'):
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_csp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-solcor_csp.parset peel-group'+groups[0]+'.skymodel > '+ms+'_calpreamp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds)

        # [PARALLEL] calibrate amplitude 5 min timescale - peel-avg_TC*.MS:CORRECTED_DATA_PHASE (no correction)
        logging.info('Calibrating amplitude...')
        cmds = []
        for ms in glob.glob('peel-avg_TC*.MS'):
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_amp  '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-sol_amp.parset peel-group'+groups[0]+'.skymodel > '+ms+'_calamp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds)

        ########################################################
        # LoSoTo Amp rescaling
        # TODO: add plotting to losoto
        for ms in glob.glob('peel-avg_TC*.MS'):
            h5parm = ms.replace('.MS','.h5')
            check_rm(h5parm)
            os.system('H5parm_importer.py -i instrument_amp -v '+h5parm+' '+ms)
            os.system('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/losoto.parset > '+ms+'_losoto-c'+str(i)+'.log 2>&1')
            os.system('H5parm_exporter.py -i instrument_amp -v -c '+h5parm+' '+ms)
            check_rm(ms+'/instrument_amp')
            os.system('mv '+ms+'/sol000_instrument_amp'+' '+ms+'/instrument_amp')

        # merge parmdbs
        logging.info('Merging instrument tables...')
        for ms in glob.glob('peel-avg_TC*.MS'):
            check_rm(ms+'/instrument')
            merge_parmdb(ms+'/instrument_amp', ms+'/instrument_csp', ms+'/instrument_empty', ms+'/instrument')

        # [PARALLEL] correct amplitude - peel-avg_TC*.MS:DATA -> peel_avg_TC*.MS:CORRECTED_DATA (selfcal phase+amp corrected, beam corrected)
        logging.info('Correcting amplitude...')
        cmds = []
        for ms in glob.glob('peel-avg_TC*.MS'):
            cmds.append('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-cor_ampcsp.parset peel-group'+groups[0]+'.skymodel > '+ms+'_corampcsp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds)

# backup instrument tables
logging.info('Backup instrument tables...')
for ms in glob.glob('peel-avg_TC*.MS'):
    os.system('cp -r '+ms+'/instrument peel/'+dd['name']+'/instruments/'+ms.replace('MS','parmdb'))

# now do as for the DDcal but for the entire facet to obtain a complete image of the facet and do a final subtraction
##############################################################################################################################
# [PARALLEL] Add rest of the facet - group*_TC*.MS:MODEL_DATA (high+low resolution model)
logging.info('Add facet model...')
cmds = []
for tc in tcs:
    # models are opened one at the time, no problems
    for ms in glob.glob('group*_TC'+tc+'.MS'):
        g = re.findall(r'\d+', ms)[0] # group number
        model = 'peel/'+dd['name']+'/models/peel_facet-g'+g+'.model'
        modellr = 'peel/'+dd['name']+'/models/peel_facet-lr-g'+g+'.model'
        cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_ft.py '+ms+' '+model+' > '+ms+'_facet-ft.log 2>&1; \
            casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_ft.py '+ms+' '+modellr+' incremental >> '+ms+'_facet-ft.log 2>&1')
    thread_cmd(cmds)

# [PARALLEL] ADD and corrupt group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + facet from model, cirular, beam correcred)
logging.info('Corrupting facet model...')
cmds = []
for ms in sorted(glob.glob('group*_TC*.MS')):
    cmds.append('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-init_add.parset > '+ms+'_facet-add.log 2>&1')
thread_cmd(cmds)

# [PARALLEL] Concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:CORRECTED_DATA -> facet-avg_TC*.MS:DATA (selfcal corrected, field subtracted but facet, avg, phase shifted, beam corrected)
logging.info('Shifting+averaging facet...')
cmds = []
for tc in tcs:
    mss = glob.glob('group*_TC'+tc+'.MS')
    msout = 'facet-avg_TC'+tc+'.MS'
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' \
            shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\] > '+msout+'_facet-shiftavg.log 2>&1')
thread_cmd(cmds)

# [PARALLEL] correct amp+ph - facet-avg_TC*.MS:DATA -> facet_avg_TC*.MS:CORRECTED_DATA (selfcal phase+amp corrected)
# copy instrument table in facet dataset
for tc in tcs:
    msDD = 'peel-avg_TC'+tc+'.MS'
    msFacet = 'facet-avg_TC'+tc+'.MS'
    os.system('cp -r '+msDD+'/instrument '+msFacet+'/instrument')
logging.info('Correcting facet amplitude+phase...')
cmds = []
for ms in glob.glob('facet-avg_TC*.MS'):
    cmds.append('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-cor_ampcsp.parset peel-group'+groups[0]+'.skymodel > '+ms+'_facet-corampcsp.log 2>&1')
thread_cmd(cmds)

# [PARALLEL] averaging before cleaning peel-avg_TC*.MS:CORRECTED_DATA -> peel-avgavg_TC*.MS:DATA
logging.info('Averaging before cleaning...')
cmds = []
for ms in glob.glob('facet-avg_TC*.MS'):
    msout = ms.reaplce('avg','avgavg')
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-avg.parset msin='+ms+' msin.datacolumn=DATA msout='+msout+' avg.freqstep=4 avg.timestep=10')
thread_cmd(cmds)

# Concatenating (in time) before imaging peel-noavg_TC*.MS:DATA -> concat.MS:DATA (beam corrected, only source to peel in the data, all chan)
# TODO: can we do virtual?
check_rm('concat.MS')
logging.info('Concatenating TCs...')
os.system('casapy --nogui --log2term --nologger -c  /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_concat.py facet-avgavg_TC*.MS concat-avg.MS > concat-facet.log 2>&1')

# Image of the entire facet
logging.info('Cleaning facet...')
imagename = 'peel/'+dd['name']'/images/facet'
if dd['facet_extended']:
    os.system('casa --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_clean_peel-ext.py concat-avg.MS '+imagename+' > casaclean1-facet.log 2>&1')
    os.system('make_mask.py -t -p 7 -i 5 -r 30,10 -m '+imagename+'.mask '+imagename+'.image.tt0 > make_mask-c'+str(i)+'.log  2>&1')
    os.system('casa --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_clean_peel-ext.py concat-avg.MS '+imagename+'-masked '+imagename+'.mask > casaclean2-facet.log 2>&1')
else: 
    os.system('casa --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_clean_peel.py concat-avg.MS '+imagename+' > casaclean1-facet.log 2>&1')
    os.system('make_mask.py -p 7 -i 5 -r 30,10 -m '+imagename+'.mask '+imagename+'.image.tt0 > make_mask-c'+str(i)+'.log  2>&1')
    os.system('casa --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_clean_peel.py concat-avg.MS '+imagename+'-masked '+imagename+'.mask > casaclean2-facet.log 2>&1')

############################################################################################################################
# in the corrected_data there's still the old facet model properly corrupted
# [PARALLEL] shift original dataset -  group*_TC*.MS:CORRECTED_DATA -> group*_TC*-shifted.MS:DATA (empty, phase shifted)
logging.info('Shifting original dataset...')
cmds = []
for ms in sorted(glob.glob('group*_TC*.MS')):
    msout = ms.replace('.MS','-shift.MS')
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shift.parset msin='ms+' msin.datacolumn=CORRECTED_DATA msout='+msout+' \
            msout.datacolumn=DATA shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\] > '+msout+'_final-shift.log 2>&1')
thread_cmd(cmds)
# copy instrument table from the DD calibration inside the phase shifted full-res MS
for ms in sorted(glob.glob('peel-avg_TC*.MS')):
    for g in groups:
        msout = ms.replace('peel-avg','group'+g).replace('.MS','-shift.MS')
        os.system('cp -r '+ms+'/instrument '+msout+'')

# TODO: remove cycle
# [PARALLEL] Add new facet model - group*_TC*-shift.MS:MODEL_DATA (new facet model)
logging.info('Add new facet model...')
cmds = []
for tc in tcs:
    # models are opened one at the time, no problems
    for ms in glob.glob('group*_TC'+tc+'-shift.MS'):
        g = re.findall(r'\d+', ms)[0] # group number
        model = 'peel/'+dd['name']+'/images/facet.model'
        cmds.append('casapy --nogui --log2term --nologger -c /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/casa_ft.py '+ms+' '+model+' > '+ms+'_final-ft.log 2>&1')
    thread_cmd(cmds)

# here the best facet model is subtracted after corruption with DD solution
# [PARALLEL] SUB corrupted facet model group*_TC*-shift.MS:DATA - MODEL_DATA -> group*_TC*-shift.MS:CORRECTED_DATA (empty data + facet from model)
# TODO: check that a parmdb for an avg freq works fine
logging.info('Subtracting new facet model...')
cmds = []
for ms in sorted(glob.glob('group*_TC*-shift.MS')):
    cmds.append('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-final_sub.parset > '+ms+'_final-add.log 2>&1')
thread_cmd(cmds)

# [PARALLEL] shift original dataset -  group*_TC*.MS:CORRECTED_DATA -> group*_TC*-shifted.MS:DATA (empty, phase shifted)
logging.info('Shifting back original dataset...')
cmds = []
for ms in sorted(glob.glob('group*_TC*-shift.MS')):
    msout = ms.replace('-shift.MS','.MS')
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shift.parset msin='ms+' msin.datacolumn=SUBTRACTED_DATA msout='+msout+' \
            msout.datacolumn=SUBTRACTED_DATA shift.phasecenter=\['+str(phasecentre[0])+'deg,'+str(phasecentre[1])+'deg\] > '+msout+'_final-shift2.log 2>&1')
thread_cmd(cmds)

# TODO: make inspection image 

logging.info("Done.")
