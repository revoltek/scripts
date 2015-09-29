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
# TODO: extract coords from ms or models
dd = {'name': 'src1', 'coord':[91.733333,41.680000], 'extended': False, 'facet_extended': False, 'mask':'', 'reg': 'sou1.crtf', 'reg_facet': 'facet1.crtf'}
phasecentre = [90.833333,42.233333] # toorhbrush
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
s = Scheduler(qsub=True, max_threads=25, dry=False, max_processors=6)

# TODO: iterate on DD calibrators

#########################################################################################
# Clear
logging.info('Cleaning...')
check_rm('peel*MS') 
check_rm('facet*MS')
check_rm('*shift.MS') 
check_rm('concat*') 
check_rm('*log *last *pickle')
check_rm('plot')
check_rm('tmpCASA_*')

logging.info('Creating dirs...')
check_rm('peel/'+dd['name'])
os.makedirs('peel/'+dd['name'])
os.makedirs('peel/'+dd['name']+'/models')
os.makedirs('peel/'+dd['name']+'/images')
os.makedirs('peel/'+dd['name']+'/instruments')
os.makedirs('peel/'+dd['name']+'/plots')
os.makedirs('peel/'+dd['name']+'/h5')

logging.info('Indexing...')
allmss = sorted(glob.glob('group*_TC*.MS'))

groups = []
tcs = []
for ms in allmss:
    g = re.findall(r'\d+', ms)[0] # group number
    tc = re.findall(r'\d+', ms)[1] # time chunk number
    groups.append(g)
    tcs.append(tc)
groups = list(set(groups))
tcs = list(set(tcs))

#################################################################################################
# [PARALLEL] Blank unwanted part of models
logging.info('Splitting skymodels...')
modeldir = 'peel/'+dd['name']+'/models'
for model in sorted(glob.glob('self/models/wide*-g*model*')):
    g = re.findall(r'\d+', model)[0] # group number
    logging.debug('Splitting group: '+g+' ('+model+')')
    peelmodel = os.path.basename(model).replace('wide','peel')
    # tmp directory are created to run CASA inside
    os.system('mkdir '+modeldir+'/'+peelmodel+'-tmp')
    os.system('cp -r '+model+' '+modeldir+'/'+peelmodel+'-tmp/'+peelmodel)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
            params={'imgs':peelmodel, 'region':os.getcwd()+'/'+dd['reg'], 'inverse':True}, wkd=modeldir+'/'+peelmodel+'-tmp', log='split_skymodels1-g'+g+'.log')
    # re-do per facet
    peelmodel = os.path.basename(model).replace('wide','peel_facet')
    os.system('mkdir '+modeldir+'/'+peelmodel+'-tmp')
    os.system('cp -r '+model+' '+modeldir+'/'+peelmodel+'-tmp/'+peelmodel)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
            params={'imgs':peelmodel, 'region':os.getcwd()+'/'+dd['reg_facet'], 'inverse':True}, wkd=modeldir+'/'+peelmodel+'-tmp', log='split_skymodels2-g'+g+'.log')
s.run(check=True) # TEST parallel

# Add DD cal model - group*_TC*.MS:MODEL_DATA (high+low resolution model)
logging.info('Add DD calibrator...')
for g in groups:
    tmpdir = os.getcwd()+'/'+'peel/'+dd['name']+'/models/peel-g'+g+'.model-tmp'
    model = tmpdir+'/peel-g'+g+'.model'
    concat = tmpdir+'/concat.MS'
    check_rm(concat+'*')
    pt.msutil.msconcat(sorted(glob.glob('group'+g+'_TC*.MS')), concat, concatTime=False)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat, 'model':model}, wkd=tmpdir, log='init-g'+g+'-ft1.log')
s.run(check=True) # TEST parallel

logging.info('Add DD calibrator (lr)...')
for g in groups:
    tmpdir = os.getcwd()+'/'+'peel/'+dd['name']+'/models/peel-lr-g'+g+'.model-tmp'
    model = tmpdir+'/peel-lr-g'+g+'.model'
    concat = tmpdir+'/concat.MS'
    check_rm(concat+'*')
    pt.msutil.msconcat(sorted(glob.glob('group'+g+'_TC*.MS')), concat, concatTime=False)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat, 'model':model, 'incr':True}, wkd=tmpdir, log='init-g'+g+'-ft2.log')
s.run(check=True) # TEST parallel

# cleanup the tmp dirs
for d in sorted(glob.glob(modeldir+'/*-tmp')):
    os.system('mv '+d+'/*model '+modeldir)
check_rm(modeldir+'/*-tmp')

###########################################################################################################
# [PARALLEL] ADD (corrupt) group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + DD cal from model, cirular, beam correcred)
# TODO: modify to account for changes in pipeline-self
logging.info('Add and corrupt model...')
for ms in allmss:
    s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-init_add.parset', \
            log=ms+'_init-addcor.log', cmd_type='BBS')
s.run(check=True)

# [PARALLEL] concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:CORRECTED_DATA -> peel_TC*.MS:DATA (empty+DD, avg, phase shifted)
# [PARALLEL] concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:MODEL_DATA -> peel-model_TC*.MS:DATA (DD model, avg, phase shifted)
logging.info('Shifting+averaging...')
for tc in tcs:
    logging.debug('Time chunk (DATA): '+tc)
    mss = glob.glob('group*_TC'+tc+'.MS')
    msout = 'peel_TC'+tc+'.MS'
    s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' msin.datacolumn=CORRECTED_DATA \
            shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\]', log=msout+'_init-shiftavg.log', cmd_type='NDPPP')
s.run(check=True)
for tc in tcs:
    logging.debug('Time chunk (MODEL_DATA): '+tc)
    mss = glob.glob('group*_TC'+tc+'.MS')
    msout = 'peel-model_TC'+tc+'.MS'
    s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' msin.datacolumn=MODEL_DATA \
            shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\]', log=msout+'_init-shiftavg.log', cmd_type='NDPPP')
s.run(check=True)

peelmss = sorted(glob.glob('peel_TC*.MS'))

# Copy the phase-shifted MODEL_DATA - peel-model_TC*.MS':DATA -> peel_TC*.MS':MODEL_DATA
logging.info('Copy MODEL_DATA...')
for ms in peelmss:
    s.add('addcol2ms.py -i '+ms+' -o MODEL_DATA', log=ms+'_init-addcol.log', cmd_type='python', processors='max')
s.run(check=True)
for ms in peelmss:
    msmodel = ms.replace('peel', 'peel-model')
    logging.debug(msmodel+':DATA -> '+ms+':MODEL_DATA')
    s.add('taql "update '+ms+', '+msmodel+' as model set MODEL_DATA=model.DATA"', log=msout+'_init-taql.log', cmd_type='general')
s.run(check=True)
check_rm('peel-model_TC*.MS')

# BL avg 
logging.info('BL-based averaging...')
for ms in peelmss:
    s.add('BLavg.py -m '+ms, log=ms+'_smooth.log', cmd_type='python')
s.run(check=True)
BLavgpeelmss = sorted(glob.glob('peel_TC*-BLavg.MS'))
sys.exit(1)

###################################################################################################################
# self-cal cycle
for i in xrange(4):
    logging.info('Start peel cycle: '+str(i))

    # ft model - peel_TC*.MS:MODEL_DATA (best available model)
    # on cycle0 it is already into the MODEL_DATA and it comes from the selfcal image
    if i > 0:
        logging.info('FT model...')
        for ms in BLavgpeelmss:
            s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':ms, 'model':imagename+'-masked.model'}, log=ms+'_ft-c'+str(i)+'.log')
            s.run(check=True) # no parallel (problem multiple accesses to model file)

    if i < 2:
        ################################################################################################
        # [PARALLEL] calibrate phase-only - peel_TC*.MS:DATA -> peel_TC*.MS:CORRECTED_DATA
        logging.info('Calibrating TEC...')
        for ms in BLavgpeelmss:
            s.add('calibrate-stand-alone -f '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-sol.parset '+skymodel, \
                    log=ms+'_sol-c'+str(i)+'.log', cmd_type='BBS')
        s.run(check=True)
        for ms in peelmss:
            check_rm(ms+'/instrument')
            os.system('cp -r '+ms.replace('.MS','-BLavg.MS')+'/instrument '+ms+'/instrument')

        logging.info('Correcting TEC...')
        for ms in peelmss:
            s.add('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-cor.parset '+skymodel, \
                    log=ms+'_cor-c'+str(i)+'.log', cmd_type='BBS')
        s.run(check=True)
    else:
        ##################################################################################################
        # [PARALLEL] calibrate phase-only - peel_TC*.MS:DATA -> peel_TC*.MS:CORRECTED_DATA_PHASE
        logging.info('Calibrating phase...')
        for ms in peelmss:
            s.add('calibrate-stand-alone -f --parmdb-name instrument_tec '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-solcor_tec.parset '+skymodel, \
                    log=ms+'_calpreamp-c'+str(i)+'.log', cmd_type='BBS')
        s.run(check=True)

        # [PARALLEL] calibrate amplitude 5 min timescale every 20 SBs - peel_TC*.MS:CORRECTED_DATA_PHASE (no correction)
        logging.info('Calibrating amplitude...')
        for ms in peelmss:
            s.add('calibrate-stand-alone -f --parmdb-name instrument_amp  '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-sol_amp.parset '+skymodel, \
                    log=ms+'_calamp-c'+str(i)+'.log', cmd_type='BBS', processors = 'max')
        s.run(check=True)

        # merge parmdbs
        logging.info('Merging instrument tables...')
        for ms in peelmss:
            merge_parmdb(ms+'/instrument_tec', ms+'/instrument_amp', ms+'/instrument', clobber=True)

        ########################################################
        # LoSoTo Amp rescaling
        os.makedirs('plot')
        check_rm('globaldb')
        os.makedirs('globaldb')
        for num, ms in enumerate(peelmss):
            os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
            if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
        h5parm = 'global-c'+str(i)+'.h5'

        s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
        s.run(check=False)
        s.add('losoto -v '+h5parm+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/losoto.parset', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
        s.run(check=False)
        s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(i)+'.log', log_append=True, cmd_type='python')
        s.run(check=True)

        for num, ms in enumerate(peelmss):
            check_rm(ms+'/instrument')
            os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/instrument')
        os.system('mv plot peel/'+dd['name']+'/plots/plot-c'+str(i))
        os.system('mv '+h5parm+' peel/'+dd['name']+'/h5')

        # [PARALLEL] correct phase + amplitude - peel_TC*.MS:DATA -> peel_TC*.MS:CORRECTED_DATA (selfcal TEC+ph+amp corrected)
        logging.info('Correcting phase+amplitude...')
        for ms in peelmss:
            s.add('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-cor_amptec.parset '+skymodel, \
                    log=ms+'_coramptec-c'+str(i)+'.log', cmd_type='BBS')
        s.run(check=True)

    # [PARALLEL] averaging before cleaning peel_TC*.MS:CORRECTED_DATA -> peel_TC*-avg.MS:DATA
    check_rm('peel_TC*-avg.MS')
    logging.info('Averaging before cleaning...')
    nchan = find_nchan(peelmss[0])
    for ms in peelmss:
        msout = ms.replace('.MS','-avg.MS')
        s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-avg.parset msin='+ms+' msin.nchan='+str(nchan-nchan%4)+' msin.datacolumn=CORRECTED_DATA \
                msout='+msout+' avg.freqstep=4 avg.timestep=10', log=ms+'_avgclean-c'+str(i)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # Concatenating (in time) before imaging peel_TC*-avg.MS:DATA -> concat.MS:DATA (beam corrected, only source to peel in the data, all chan)
    check_rm('concat.MS*')
    logging.info('Concatenating TCs...')
    pt.msutil.msconcat(sorted(glob.glob('peel_TC*-avg.MS')), 'concat.MS', concatTime=False)

    # Clean mask clean
    imagename = 'peel/'+dd['name']+'/images/peel-'+str(i)
    imsize = size_from_facet('peel/'+dd['name']+'/models/peel-g'+groups[0]+'.model', dd['coord'], 3)
    logging.debug('Image size set to '+str(imsize))
    if dd['extended']: multiscale = [0,3,9,18]
    else: multiscale = []
    logging.info('Cleaning...')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', \
            params={'msfile':'concat.MS', 'imagename':imagename, 'imsize':imsize, 'niter':1000, 'multiscale':multiscale, 'wproj':128}, log='casaclean1-c'+str(i)+'.log')
    s.run(check=True)
    make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask', atrous_do=dd['extended'], threshisl=5)
    logging.info('Cleaning (with mask)...')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', \
            params={'msfile':'concat.MS', 'imagename':imagename+'-masked', 'imsize':imsize, 'niter':1000, 'multiscale':multiscale, 'wproj':128, 'mask':imagename+'.newmask'}, log='casaclean2-c'+str(i)+'.log')
    s.run(check=True)

    # TODO: flag residuals

# backup instrument tables
logging.info('Back up instrument tables...')
for ms in peelmss:
    logging.debug('Creating: peel/'+dd['name']+'/instruments/'+ms.replace('MS','parmdb'))
    os.system('cp -r '+ms+'/instrument peel/'+dd['name']+'/instruments/'+ms.replace('MS','parmdb'))

# now do the same but for the entire facet to obtain a complete image of the facet and do a final subtraction
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
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':model, 'wproj':512}, log='facet-g'+g+'-ft1.log')
    s.run(check=True) # no parallel
    logging.debug('Ft low res model')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':modellr, 'incr':True, 'wproj':512}, log='facet-g'+g+'-ft2.log')
    s.run(check=True) # no parallel

# [PARALLEL] ADD and corrupt group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + facet from model, cirular, beam correcred)
logging.info('Add and corrupt facet model...')
for ms in allmss:
    s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-init_add.parset '+skymodel, \
        log=ms+'_facet-add.log', cmd_type='BBS')
s.run(check=True)

# [PARALLEL] Concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:CORRECTED_DATA -> facet_TC*.MS:DATA (not corrected, field subtracted but facet, avg, phase shifted)
logging.info('Shifting+averaging facet...')
for tc in tcs:
    mss = glob.glob('group*_TC'+tc+'.MS')
    msout = 'facet_TC'+tc+'.MS'
    s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' msin.datacolumn=CORRECTED_DATA \
            shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\]', \
            log=msout+'_facet-shiftavg.log', cmd_type='NDPPP')
s.run(check=True)

facetmss = sorted(glob.glob('facet_TC*.MS'))

# [PARALLEL] correct amp+ph - facet_TC*.MS:DATA -> facet_TC*.MS:CORRECTED_DATA (selfcal phase+amp corrected)
# copy instrument table in facet dataset
for tc in tcs:
    msDD = 'peel_TC'+tc+'.MS'
    msFacet = 'facet_TC'+tc+'.MS'
    logging.debug(msDD+'/instrument -> '+msFacet+'/instrument')
    os.system('cp -r '+msDD+'/instrument '+msFacet+'/instrument')
logging.info('Correcting facet amplitude+phase...')
for ms in facetmss:
    s.add('calibrate-stand-alone '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-cor_amptec-facet.parset '+skymodel, \
            log=ms+'_facet-coramptec.log', cmd_type='BBS')
s.run(check=True)

# [PARALLEL] averaging before cleaning facet_TC*.MS:CORRECTED_DATA -> facet_TC*-avg.MS:DATA
logging.info('Averaging before cleaning...')
nchan = find_nchan(facetmss[0])
for ms in facetmss:
    msout = ms.replace('.MS','-avg.MS')
    s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-avg.parset msin='+ms+' msin.nchan='+str(nchan-nchan%4)+' msout='+msout+'\
            avg.freqstep=2 avg.timestep=5', log=ms+'_facet-avgclean.log', cmd_type='NDPPP')
s.run(check=True)

# Concatenating (in time) before imaging facet_TC*-avg.MS -> concat.MS (beam corrected, only source to peel in the data, all chan)
check_rm('concat.MS*')
logging.info('Concatenating TCs...')
pt.msutil.msconcat(sorted(glob.glob('facet_TC*-avg.MS')), 'concat.MS', concatTime=False)

# Image of the entire facet
imagename = 'peel/'+dd['name']+'/images/facet'
imsize = size_from_facet('peel/'+dd['name']+'/models/peel_facet-g'+g+'.model', dd['coord'], 3)
logging.debug('Image size set to '+str(imsize))
if dd['extended']: multiscale = [0,3,9,18]
else: multiscale = []
logging.info('Cleaning facet...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', \
        params={'msfile':'concat.MS', 'imagename':imagename, 'imsize':imsize, 'niter':5000, 'multiscale':multiscale, 'wproj':512}, log='casaclean1-facet.log')
s.run(check=True)
make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask', atrous_do=dd['facet_extended'])
logging.info('Cleaning facet (with mask)...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', \
        params={'msfile':'concat.MS', 'imagename':imagename+'-masked', 'imsize':imsize, 'niter':5000, 'multiscale':multiscale, 'wproj':512, 'mask':imagename+'.newmask'}, log='casaclean2-facet.log')
s.run(check=True)

sys.exit(1)

############################################################################################################################
# in the corrected_data there's still the old facet model properly corrupted
# [PARALLEL] shift original dataset -  group*_TC*.MS:CORRECTED_DATA -> group*_TC*-shifted.MS:DATA (empty+facet, phase shifted)
logging.info('Shifting original dataset...')
for ms in allmss:
    msout = ms.replace('.MS','-shift.MS')
    s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shift.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA msout='+msout+' \
            msout.datacolumn=DATA shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\]', \
            log=msout+'_final-shift.log', cmd_type='NDPPP')
s.run(check=True)
# copy instrument table from the DD calibration inside the phase shifted full-res MS
for ms in peelmss:
    for g in groups:
        msout = ms.replace('peel','group'+g).replace('.MS','-shift.MS')
        logging.debug(ms+'/instrument -> '+msout+'/instrument')
        os.system('cp -r '+ms+'/instrument '+msout+'')

# Add new facet model - group*_TC*-shift.MS:MODEL_DATA (new facet model)
logging.info('Add new facet model...')
check_rm('concat.MS*')
pt.msutil.msconcat(sorted(glob.glob('group*_TC*-shift.MS')), 'concat.MS', concatTime=False)
model = 'peel/'+dd['name']+'/images/facet.model'
s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':model, 'wproj':512}, log='final-ft.log')
s.run(check=True)

# here the best facet model is subtracted after corruption with DD solution
# [PARALLEL] SUB corrupted facet model group*_TC*-shift.MS:DATA - MODEL_DATA -> group*_TC*-shift.MS:CORRECTED_DATA (empty data + facet from model)
# TODO: check that a parmdb for an avg freq works fine
logging.info('Subtracting new facet model...')
for ms in sorted(glob.glob('group*_TC*-shift.MS')):
    s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-final_sub.parset', \
            log=ms+'_final-add.log', cmd_type='BBS')
s.run(check=True)

# [PARALLEL] shift back dataset -  group*_TC*-shifted.MS:CORRECTED_DATA -> group*_TC*.MS:DATA (empty, phase shifted)
logging.info('Shifting back original dataset...')
for ms in sorted(glob.glob('group*_TC*-shift.MS')):
    msout = ms.replace('-shift.MS','.MS')
    s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-shift.parset msin='+ms+' msin.datacolumn=SUBTRACTED_DATA msout='+msout+' \
            msout.datacolumn=SUBTRACTED_DATA shift.phasecenter=\['+str(phasecentre[0])+'deg,'+str(phasecentre[1])+'deg\]', \
            log=msout+'_final-shift2.log', cmd_type='NDPPP')
s.run(check=True)
check_rm('group*_TC*-shift.MS') # otherwise next wildcard select them

# Make inspection image 
logging.info('Inspection image...')
check_rm('concat.MS*')
pt.msutil.msconcat(allmss, 'concat.MS', concatTime=False)
imagename = 'peel/'+dd['name']+'/images/inspection'
s.add('/opt/cep/WSClean/wsclean-1.7/build/wsclean -datacolumn SUBTRACTED_DATA -reorder -name ' + imagename + ' -size 4000 4000 \
           -scale 15arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -no-update-model-required -maxuv-l 2500 concat.MS', \
           log='final-wsclean.log', cmd_type='wsclean')
s.run(check=True)

logging.info("Done.")
