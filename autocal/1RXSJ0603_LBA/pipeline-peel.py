#!/usr/bin/python
# perform peeling on a group of TCs. Script should be run in the directory with MSs inside.
# widefield-group#.skymodel is required for each group
# Input:
# * group#_TC###.MS must have instrument table with amp+phase+TEC and beam uncorrected "empty" data in SUBTRACTED_DATA
# it should be possible to collect MSs after the run of pipeline-self.py
# * .model images are in the "self/models" dir made by pipeline-self
# * list of regions
# Output:
# DD-calibrated and subtracted data in SUBTRACTED_DATA
# Image of the facet


# coordinate of the region to find the DD
#coord = [90.833333,42.233333] # toorhbrush
#coord = [91.733333,41.680000] # strong pts
dd = {'name': 'src1', 'coord':[91.733333,41.680000], 'extended': False, 'facet_extended': False, 'mask':'', 'reg': 'sou1.crtf', 'reg_facet': 'facet1.crtf'}
phasecentre = [90.833333,42.233333] # toorhbrush
max_threads = 40
skymodel = '/home/fdg/scripts/autocal/1RXSJ0603_LBA/toothbrush.GMRT150.skymodel' # used only to run bbs, not important the content

##########################################################################################

import sys, os, glob, re
from lofar import bdsm
import pyrap.tables as pt
import numpy as np
import lsmtool
from lib_pipeline import *
from make_mask import make_mask

set_logger()

# TODO: iterate on DD calibrators

#########################################################################################
# Clear
logging.info('Cleaning...')
check_rm('peel*MS') 
check_rm('facet*MS')
check_rm('*shift.MS') 
check_rm('concat*') 
check_rm('*log *last')
check_rm('plot')
check_rm('peel/'+dd['name'])
os.makedirs('peel/'+dd['name'])
os.makedirs('peel/'+dd['name']+'/models')
os.makedirs('peel/'+dd['name']+'/images')
os.makedirs('peel/'+dd['name']+'/instruments')
os.makedirs('peel/'+dd['name']+'/plots')
os.makedirs('peel/'+dd['name']+'/h5')

groups = []
tcs = []
for ms in sorted(glob.glob('group*_TC*.MS')):
    g = re.findall(r'\d+', ms)[0] # group number
    tc = re.findall(r'\d+', ms)[1] # time chunk number
    groups.append(g)
    tcs.append(tc)
groups = list(set(groups))
tcs = list(set(tcs))

#################################################################################################
# Blank unwanted part of models
logging.info('Splitting skymodels...')
for i, model in enumerate(sorted(glob.glob('self/models/wide*-g*model*'))):
    g = re.findall(r'\d+', model)[0] # group number
    logging.debug('Splitting group: '+g+' ('+model+')')
    peelmodel = os.path.basename(model).replace('wide','peel')
    os.system('cp -r '+model+' peel/'+dd['name']+'/models/'+peelmodel)
    run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_blank.py', params={'imgs':'peel/'+dd['name']+'/models/'+peelmodel, 'region':dd['reg'], 'inverse':True}, log='split_skymodels1-g'+g+'.log')
    # redo the same but leaving the entire facet in the model
    peelmodel = os.path.basename(model).replace('wide','peel_facet')
    os.system('cp -r '+model+' peel/'+dd['name']+'/models/'+peelmodel)
    run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_blank.py', params={'imgs':'peel/'+dd['name']+'/models/'+peelmodel, 'region':dd['reg_facet'], 'inverse':True}, log='split_skymodels2-g'+g+'.log')

###############################################################################################################################
# Add DD cal model - group*_TC*.MS:MODEL_DATA (high+low resolution model)
logging.info('Add DD calibrator...')
for g in groups:
    model = 'peel/'+dd['name']+'/models/peel-g'+g+'.model'
    modellr = 'peel/'+dd['name']+'/models/peel-lr-g'+g+'.model'
    logging.debug('Concatenating group: '+g)
    check_rm('concat.MS*')
    pt.msutil.msconcat(sorted(glob.glob('group'+g+'_TC*.MS')), 'concat.MS', concatTime=False)
    run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':model}, log='init-g'+g+'-ft1.log')
    run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':modellr, 'incr':True}, log='init-g'+g+'-ft2.log')

# [PARALLEL] ADD (corrupt) group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + DD cal from model, cirular, beam correcred)
logging.info('Add and corrupt model...')
cmds = []
for ms in sorted(glob.glob('group*_TC*.MS')):
    cmds.append('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-init_add.parset > '+ms+'_init-addcor.log 2>&1')
thread_cmd(cmds, max_threads)

# [PARALLEL] concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:CORRECTED_DATA -> peel-avg_TC*.MS:DATA (empty+DD, avg, phase shifted)
# [PARALLEL] concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:MODEL_DATA -> peel-avg-model_TC*.MS:DATA (DD model, avg, phase shifted)
logging.info('Shifting+averaging...')
cmds = []
for tc in tcs:
    logging.debug('Time chunk (DATA): '+tc)
    mss = glob.glob('group*_TC'+tc+'.MS')
    msout = 'peel-avg_TC'+tc+'.MS'
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' msin.datacolumn=CORRECTED_DATA \
            shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\] > '+msout+'_init-shiftavg.log 2>&1')
for tc in tcs:
    logging.debug('Time chunk (MODEL_DATA): '+tc)
    mss = glob.glob('group*_TC'+tc+'.MS')
    msout = 'peel-avg-model_TC'+tc+'.MS'
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' msin.datacolumn=MODEL_DATA \
            shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\] > '+msout+'_init-shiftavg.log 2>&1')
thread_cmd(cmds, max_threads)

# Copy the phase-shifted MODEL_DATA - peel-avg-model_TC*.MS':DATA -> peel-avg_TC*.MS':MODEL_DATA
logging.info('Copy MODEL_DATA...')
for ms in glob.glob('peel-avg-model_TC*.MS'):
    msout = ms.replace('peel-avg-model','peel-avg')
    logging.debug(ms+':DATA -> '+msout+':MODEL_DATA')
    os.system('addcol2ms.py -i '+msout+' -o MODEL_DATA')
    os.system('taql \'update '+msout+', '+ms+' as model set MODEL_DATA=model.DATA\'')
check_rm('peel-avg-model_TC*.MS')

# [PARALLEL] create a fake parmdb to be used later for merging slow-amp and fast-phase parmdbs
logging.info('Creating fake parmdb...')
cmds = []
for ms in glob.glob('peel-avg_TC*.MS'):
    cmds.append('calibrate-stand-alone -f --parmdb-name instrument_empty '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-fakeparmdb.parset '+skymodel+' > '+ms+'_init-fakeparmdb.log 2>&1')
thread_cmd(cmds, max_threads)

###################################################################################################################
# self-cal cycle
for i in xrange(4):
    logging.info('Start peel cycle: '+str(i))

    # ft model - peel-avg_TC*.MS:MODEL_DATA (best available model)
    # on cycle0 it is already into the MODEL_DATA and it comes from the selfcal image
    if i > 0:
        logging.info('FT model...')
        for ms in glob.glob('peel-avg_TC*.MS'):
            run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':ms, 'model':imagename+'-masked.model'}, log=ms+'_ft-c'+str(i)+'.log')

    if i < 2:
        ################################################################################################
        # [PARALLEL] calibrate phase-only - peel-avg_TC*.MS:DATA -> peel-avg_TC*.MS:CORRECTED_DATA
        logging.info('Calibrating+correcting phase...')
        cmds = []
        for ms in glob.glob('peel-avg_TC*.MS'):
            cmds.append('calibrate-stand-alone -f '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-solcor.parset '+skymodel+' > '+ms+'_solcor-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads)

    else:
        ##################################################################################################
        # [PARALLEL] calibrate phase-only - peel-avg_TC*.MS:DATA -> peel-avg_TC*.MS:CORRECTED_DATA_PHASE
        logging.info('Calibrating phase...')
        cmds = []
        for ms in glob.glob('peel-avg_TC*.MS'):
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_csp '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-solcor_csp.parset '+skymodel+' > '+ms+'_calpreamp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads)

        # [PARALLEL] calibrate amplitude 5 min timescale - peel-avg_TC*.MS:CORRECTED_DATA_PHASE (no correction)
        logging.info('Calibrating amplitude...')
        cmds = []
        for ms in glob.glob('peel-avg_TC*.MS'):
            cmds.append('calibrate-stand-alone -f --parmdb-name instrument_amp  '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-sol_amp.parset '+skymodel+' > '+ms+'_calamp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads)

        # merge parmdbs
        logging.info('Merging instrument tables...')
        for ms in glob.glob('peel-avg_TC*.MS'):
            check_rm(ms+'/instrument')
            merge_parmdb(ms+'/instrument_amp', ms+'/instrument_csp', ms+'/instrument_empty', ms+'/instrument')

        ########################################################
        # LoSoTo Amp rescaling
        os.makedirs('plot')
        check_rm('globaldb')
        os.makedirs('globaldb')
        for num, ms in enumerate(glob.glob('peel-avg_TC*.MS')):
            os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
            if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
        h5parm = 'global-c'+str(i)+'.h5'
        os.system('H5parm_importer.py -v '+h5parm+' globaldb > losoto-c'+str(i)+'.log 2>&1')
        os.system('losoto.py -v '+h5parm+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/losoto.parset >> losoto-c'+str(i)+'.log 2>&1')
        os.system('H5parm_exporter.py -v -c '+h5parm+' globaldb >> losoto-c'+str(i)+'.log 2>&1')
        for num, ms in enumerate(glob.glob('peel-avg_TC*.MS')):
            check_rm(ms+'/instrument')
            os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
        os.system('mv plot peel/'+dd['name']+'/plots/plot-c'+str(i))
        os.system('mv '+h5parm+' peel/'+dd['name']+'/h5')

        # [PARALLEL] correct phase + amplitude - peel-avg_TC*.MS:DATA -> peel_avg_TC*.MS:CORRECTED_DATA (selfcal phase+amp corrected, beam corrected)
        logging.info('Correcting phase+amplitude...')
        cmds = []
        for ms in glob.glob('peel-avg_TC*.MS'):
            cmds.append('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-cor_ampcsp.parset '+skymodel+' > '+ms+'_corampcsp-c'+str(i)+'.log 2>&1')
        thread_cmd(cmds, max_threads)

    # [PARALLEL] averaging before cleaning peel-avg_TC*.MS:CORRECTED_DATA -> peel-avgavg_TC*.MS:DATA
    check_rm('peel-avgavg_TC*.MS')
    logging.info('Averaging before cleaning...')
    cmds = []
    for ms in glob.glob('peel-avg_TC*.MS'):
        msout = ms.replace('avg','avgavg')
        cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-avg.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA msout='+msout+' avg.freqstep=4 avg.timestep=10 > '+ms+'_avgclean-c'+str(i)+'.log 2>&1')
    thread_cmd(cmds, max_threads)

    # Concatenating (in time) before imaging peel-avgavg_TC*.MS:DATA -> concat.MS:DATA (beam corrected, only source to peel in the data, all chan)
    check_rm('concat.MS*')
    logging.info('Concatenating TCs...')
    pt.msutil.msconcat(glob.glob('peel-avgavg_TC*.MS'), 'concat.MS', concatTime=False)

    # Clean mask clean
    imagename = 'peel/'+dd['name']+'/images/peel-'+str(i)
    imsize = size_from_facet('peel/'+dd['name']+'/models/peel-g'+groups[0]+'.model', dd['coord'], 3)
    logging.debug('Image size set to '+str(imsize))
    if dd['extended']: multiscale = [0,3,9,18]
    else: multiscale = []
    logging.info('Cleaning...')
    run_casa(command='/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', params={'msfile':'concat.MS', 'imagename':imagename, 'imsize':imsize, 'niter':1000, 'multiscale':multiscale, 'wproj':128}, log='casaclean1-c'+str(i)+'.log')
    make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask', atrous_do=dd['extended'], threshisl=5)
    logging.info('Cleaning (with mask)...')
    run_casa(command='/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', params={'msfile':'concat.MS', 'imagename':imagename+'-masked', 'imsize':imsize, 'niter':1000, 'multiscale':multiscale, 'wproj':128, 'mask':imagename+'.newmask'}, log='casaclean2-c'+str(i)+'.log')

    # TODO: flag residuals

# backup instrument tables
logging.info('Back up instrument tables...')
for ms in glob.glob('peel-avg_TC*.MS'):
    logging.debug('Creating: peel/'+dd['name']+'/instruments/'+ms.replace('MS','parmdb'))
    os.system('cp -r '+ms+'/instrument peel/'+dd['name']+'/instruments/'+ms.replace('MS','parmdb'))

# now do as for the DDcal but for the entire facet to obtain a complete image of the facet and do a final subtraction
##############################################################################################################################
# Add rest of the facet - group*_TC*.MS:MODEL_DATA (high+low resolution model)
logging.info('Add facet model...')
for g in groups: 
    model = 'peel/'+dd['name']+'/models/peel_facet-g'+g+'.model'
    modellr = 'peel/'+dd['name']+'/models/peel_facet-lr-g'+g+'.model'
    logging.debug('Concatenating group: '+g)
    check_rm('concat.MS*')
    pt.msutil.msconcat(sorted(glob.glob('group'+g+'_TC*.MS')), 'concat.MS', concatTime=False)
    logging.debug('Ft high res model')
    run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':model, 'wproj':512}, log='facet-g'+g+'-ft1.log')
    logging.debug('Ft low res model')
    run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':modellr, 'incr':True, 'wproj':512}, log='facet-g'+g+'-ft2.log')

# [PARALLEL] ADD and corrupt group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + facet from model, cirular, beam correcred)
logging.info('Add and corrupt facet model...')
cmds = []
for ms in sorted(glob.glob('group*_TC*.MS')):
    cmds.append('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-init_add.parset > '+ms+'_facet-add.log 2>&1')
thread_cmd(cmds, max_threads)

# [PARALLEL] Concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:CORRECTED_DATA -> facet-avg_TC*.MS:DATA (selfcal corrected, field subtracted but facet, avg, phase shifted, beam corrected)
logging.info('Shifting+averaging facet...')
cmds = []
for tc in tcs:
    mss = glob.glob('group*_TC'+tc+'.MS')
    msout = 'facet-avg_TC'+tc+'.MS'
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' msin.datacolumn=CORRECTED_DATA \
            shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\] > '+msout+'_facet-shiftavg.log 2>&1')
thread_cmd(cmds, max_threads)

# [PARALLEL] correct amp+ph - facet-avg_TC*.MS:DATA -> facet_avg_TC*.MS:CORRECTED_DATA (selfcal phase+amp corrected)
# copy instrument table in facet dataset
for tc in tcs:
    msDD = 'peel-avg_TC'+tc+'.MS'
    msFacet = 'facet-avg_TC'+tc+'.MS'
    logging.debug(msDD+'/instrument -> '+msFacet+'/instrument')
    os.system('cp -r '+msDD+'/instrument '+msFacet+'/instrument')
logging.info('Correcting facet amplitude+phase...')
cmds = []
for ms in glob.glob('facet-avg_TC*.MS'):
    cmds.append('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-cor_ampcsp.parset '+skymodel+' > '+ms+'_facet-corampcsp.log 2>&1')
thread_cmd(cmds, max_threads)

# [PARALLEL] averaging before cleaning facet-avg_TC*.MS:CORRECTED_DATA -> facet-avgavg_TC*.MS:DATA
logging.info('Averaging before cleaning...')
cmds = []
for ms in glob.glob('facet-avg_TC*.MS'):
    msout = ms.replace('avg','avgavg')
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-avg.parset msin='+ms+' msout='+msout+' avg.freqstep=2 avg.timestep=5 > '+ms+'_facet-avgclean.log  2>&1')
thread_cmd(cmds, max_threads)

# Concatenating (in time) before imaging facet-avgavg_TC*.MS -> concat.MS (beam corrected, only source to peel in the data, all chan)
check_rm('concat.MS*')
logging.info('Concatenating TCs...')
pt.msutil.msconcat(glob.glob('facet-avgavg_TC*.MS'), 'concat.MS', concatTime=False)

# Image of the entire facet
imagename = 'peel/'+dd['name']+'/images/facet'
imsize = size_from_facet('peel/'+dd['name']+'/models/peel_facet-g'+g+'.model', dd['coord'], 3)
logging.debug('Image size set to '+str(imsize))
if dd['extended']: multiscale = [0,3,9,18]
else: multiscale = []
logging.info('Cleaning facet...')
run_casa(command='/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', params={'msfile':'concat.MS', 'imagename':imagename, 'imsize':imsize, 'niter':5000, 'multiscale':multiscale, 'wproj':512}, log='casaclean1-facet.log')
make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask', atrous_do=dd['facet_extended'])
logging.info('Cleaning facet (with mask)...')
run_casa(command='/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', params={'msfile':'concat.MS', 'imagename':imagename+'-masked', 'imsize':imsize, 'niter':5000, 'multiscale':multiscale, 'wproj':512, 'mask':imagename+'.newmask'}, log='casaclean2-facet.log')

############################################################################################################################
# in the corrected_data there's still the old facet model properly corrupted
# [PARALLEL] shift original dataset -  group*_TC*.MS:CORRECTED_DATA -> group*_TC*-shifted.MS:DATA (empty, phase shifted)
logging.info('Shifting original dataset...')
cmds = []
for ms in sorted(glob.glob('group*_TC*.MS')):
    msout = ms.replace('.MS','-shift.MS')
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shift.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA msout='+msout+' \
            msout.datacolumn=DATA shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\] > '+msout+'_final-shift.log 2>&1')
thread_cmd(cmds, max_threads)
# copy instrument table from the DD calibration inside the phase shifted full-res MS
for ms in sorted(glob.glob('peel-avg_TC*.MS')):
    for g in groups:
        logging.debug(ms+'/instrument -> '+msout+'/instrument')
        msout = ms.replace('peel-avg','group'+g).replace('.MS','-shift.MS')
        os.system('cp -r '+ms+'/instrument '+msout+'')

# Add new facet model - group*_TC*-shift.MS:MODEL_DATA (new facet model)
logging.info('Add new facet model...')
check_rm('concat.MS*')
pt.msutil.msconcat(sorted(glob.glob('group*_TC*-shift.MS')), 'concat.MS', concatTime=False)
model = 'peel/'+dd['name']+'/images/facet.model'
run_casa(command='/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':model, 'wproj':512}, log='final-ft.log')

# here the best facet model is subtracted after corruption with DD solution
# [PARALLEL] SUB corrupted facet model group*_TC*-shift.MS:DATA - MODEL_DATA -> group*_TC*-shift.MS:CORRECTED_DATA (empty data + facet from model)
# TODO: check that a parmdb for an avg freq works fine
logging.info('Subtracting new facet model...')
cmds = []
for ms in sorted(glob.glob('group*_TC*-shift.MS')):
    cmds.append('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-final_sub.parset > '+ms+'_final-add.log 2>&1')
thread_cmd(cmds, max_threads)

# [PARALLEL] shift original dataset -  group*_TC*.MS:CORRECTED_DATA -> group*_TC*-shifted.MS:DATA (empty, phase shifted)
logging.info('Shifting back original dataset...')
cmds = []
for ms in sorted(glob.glob('group*_TC*-shift.MS')):
    msout = ms.replace('-shift.MS','.MS')
    cmds.append('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shift.parset msin='+ms+' msin.datacolumn=SUBTRACTED_DATA msout='+msout+' \
            msout.datacolumn=SUBTRACTED_DATA shift.phasecenter=\['+str(phasecentre[0])+'deg,'+str(phasecentre[1])+'deg\] > '+msout+'_final-shift2.log 2>&1')
thread_cmd(cmds, max_threads)

# Make inspection image 
logging.info('Inspection image...')
check_rm('concat.MS*')
pt.msutil.msconcat(sorted(glob.glob('group*_TC*.MS')), 'concat.MS', concatTime=False)
imagename = 'peel/'+dd['name']+'/images/inspection'
os.system('/opt/cep/WSClean/wsclean-1.7/build/wsclean -datacolumn SUBTRACTED_DATA -reorder -name ' + imagename + ' -size 4000 4000 \
           -scale 15arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -no-update-model-required -maxuv-l 2500 concat.MS > final-wsclean.log 2>&1')

logging.info("Done.")
