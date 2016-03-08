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
ddset = [{'name': 'src1', 'coord':[91.733333,41.680000], 'extended': False, 'facet_extended': False, 'mask':'', 'reg': 'src1.crtf', 'reg_facet': 'facet1.crtf', 'faint': False},
         {'name': 'src2', 'coord':[91.391897,41.530003], 'extended': False, 'facet_extended': False, 'mask':'', 'reg': 'src2.crtf', 'reg_facet': 'facet2.crtf', 'faint': False},
         {'name': 'tooth', 'coord':[90.833333,42.233333], 'extended': False, 'facet_extended': False, 'mask':'tooth_mask.crtf', 'reg': 'src3.crtf', 'reg_facet': 'facet3.crtf', 'faint': True}]
skymodel = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.GMRT150.skymodel' # used only to run bbs, not important the content
parset_dir = '/home/fdg/scripts/autocal/LBAsurvey/parset_peel'
niter = 10

##########################################################################################

import sys, os, glob, re, time
import pyrap.tables as pt
import numpy as np
from lib_pipeline import *
from make_mask import make_mask

set_logger()
check_rm('logs')
s = Scheduler(dry=False)

def clean(c, mss, dd, avgfreq=4, avgtime=10, facet=False, skip_mask=False):
    """
    c = cycle/name
    mss = list of mss to avg/clean
    """
    # averaging before cleaning *.MS:CORRECTED_DATA -> *-avg.MS:DATA
    logging.info('Averaging before cleaning...')
    nchan = find_nchan(mss[0])
    for ms in mss:
        msout = ms.replace('.MS','-avg.MS')
        check_rm(msout)
        s.add('NDPPP '+parset_dir+'/NDPPP-avg.parset msin='+ms+' msin.nchan='+str(nchan-nchan%4)+' msin.datacolumn=CORRECTED_DATA \
                msout='+msout+' avg.freqstep='+str(avgfreq)+' avg.timestep='+str(avgtime), log=ms+'_cleanavg-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)
    mssavg = [ms.replace('.MS','-avg.MS') for ms in mss]

    # TEST
    os.system('ls -d *MS')
    while True:
        if all([os.path.exists(ms) for ms in mssavg]): break
        time.sleep(5)
    os.system('ls -d *MS')

    # Concatenating (in time) before imaging *-avg.MS:DATA -> concat.MS:DATA
    check_rm('concat-avg.MS*')
    logging.info('Concatenating TCs...')
    pt.msutil.msconcat(mssavg, 'concat-avg.MS', concatTime=False)

    # set image name
    imagename = 'peel/'+dd['name']+'/images/peel-'+str(c)

    # set imsize and niter
    if facet:
        imsize = size_from_facet('peel/'+dd['name']+'/models/facet.cut', dd['coord'], 3)
        niter = 50000 # TODO: check x10
    else:
        imsize = size_from_facet('peel/'+dd['name']+'/models/dd.cut', dd['coord'], 3)
        niter = 20000 # TODO: check x10

    # get wproj scaled with pixels
    wproj = 1
    if imsize > 512:
        wproj = 64
    if imsize > 799:
        wproj = 96
    if imsize > 1023:
        wproj = 128
    if imsize > 1599:
        wproj = 256
    if imsize > 2047:
        wproj = 384
    if imsize > 3000:
        wproj = 448 
    if imsize > 4095:
        wproj = 512

    # Clean mask clean
    logging.debug('Image size set to '+str(imsize))
    if dd['extended']: multiscale = [0,3,9,18]
    else: multiscale = []
    logging.info('Cleaning (cycle: '+str(c)+')...')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/LBAsurvey/casa_clean_peel.py', \
            params={'msfile':'concat-avg.MS', 'imagename':imagename, 'imsize':imsize, 'niter':niter, 'multiscale':multiscale, 'wproj':wproj}, log='casaclean1-c'+str(c)+'.log')
    s.run(check=True)

    if skip_mask: return imagename+'.model'

    if dd['faint'] or facet:
        make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask', atrous_do=dd['extended'], threshisl=5)
    else:
        make_mask(image_name = imagename+'.image.tt0', mask_name = imagename+'.newmask', atrous_do=dd['extended'], threshisl=20)

    # if dd['mask'] is set then add it to the new mask
    if dd['mask'] != '':
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
            params={'imgs':imagename+'.newmask', 'region':dd['mask'], 'setTo':1}, log='casablank-c'+str(c)+'.log', log_append=True)
        s.run(check=True)

    logging.info('Cleaning with mask (cycle: '+str(c)+')...')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/LBAsurvey/casa_clean_peel.py', \
            params={'msfile':'concat-avg.MS', 'imagename':imagename+'-masked', 'imsize':imsize, 'niter':int(niter/2.), 'multiscale':multiscale, 'wproj':wproj, 'mask':imagename+'.newmask'}, log='casaclean2-c'+str(c)+'.log')
    s.run(check=True)

    return imagename+'-masked.model'


def losoto(c, mss, dd, parset):
    logging.info('Running LoSoTo...')
    check_rm('plots')
    os.makedirs('plots')
    check_rm('globaldb')
    os.makedirs('globaldb')

    for num, ms in enumerate(mss):
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
    os.system('mv plots peel/'+dd['name']+'/plots/plots-c'+str(c))
    os.system('mv '+h5parm+' peel/'+dd['name']+'/h5')


def peel(dd):

    logging.info('Peeling: '+dd['name'])

    #########################################################################################
    # Clear
    logging.info('Cleaning...')
    check_rm('peel*MS') 
    check_rm('facet*MS')
    check_rm('*shift.MS') 
    check_rm('concat*') 
    check_rm('*last *pickle')
    check_rm('plot')
    
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
    phasecentre = get_phase_centre(allmss[0])
    
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
    # Blank unwanted part of models
    modeldir = 'peel/'+dd['name']+'/models/'
    
    for model in glob.glob('self/models/wide*_g*model*'):
        os.system('cp -r '+model+' '+modeldir+'/'+os.path.basename(model).replace('wide','peel_dd'))
        os.system('cp -r '+model+' '+modeldir+'/'+os.path.basename(model).replace('wide','peel_facet'))
    
    logging.info('Splitting skymodels...')
    models = glob.glob(modeldir+'/peel_dd*')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
                params={'imgs':models, 'region':dd['reg'], 'inverse':True}, log='split_skymodels.log')
    s.run(check=True)
    
    logging.info('Splitting skymodels (low-resolution)...')
    models = glob.glob(modeldir+'/peel_facet*')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
                params={'imgs':models, 'region':dd['reg_facet'], 'inverse':True}, log='split_skymodels.log', log_append=True)
    s.run(check=True)

    # Ugly solution to make images of the minimal size just to find out how many pixels to use
    logging.info('Making image cuts...')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_cutimg.py', \
                params={'img':glob.glob(modeldir+'/peel_dd_g*.model.ms')[0], 'region':dd['reg'], 'out':modeldir+'/dd.cut'}, log='cut_skymodels.log')
    s.run(check=True)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_cutimg.py', \
                params={'img':glob.glob(modeldir+'/peel_facet_g*.model.ms')[0], 'region':dd['reg_facet'], 'out':modeldir+'/facet.cut'}, log='cut_skymodels.log', log_append=True)
    s.run(check=True)
    
    # Add DD cal model - group*_TC*.MS:MODEL_DATA (high+low resolution model)
    logging.info('Ft DD calibrator model...')
    for g in groups:
        # tmp directory are created to run CASA inside and prevent CASA bug when multiple istances are run in the same dir
        tmpdir = os.getcwd()+'/'+modeldir+'/tmp_'+g
        model = os.getcwd()+'/'+modeldir+'/peel_dd_g'+g+'.model.ms'
        os.makedirs(tmpdir)
        concat = tmpdir+'/concat.MS'
        check_rm(concat+'*')
        pt.msutil.msconcat(sorted(glob.glob('group'+g+'_TC*.MS')), concat, concatTime=False)
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat, 'model':model, 'wproj':512}, wkd=tmpdir, log='init-g'+g+'-ft.log')
    s.run(check=True)
    
    logging.info('Ft DD calibrator (lr) model...')
    for g in groups:
        tmpdir = os.getcwd()+'/'+modeldir+'/tmp_'+g
        model = os.getcwd()+'/'+modeldir+'/peel_dd_lr_g'+g+'.model.ms'
        concat = tmpdir+'/concat.MS'
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat, 'model':model, 'wproj':512, 'incr':True}, wkd=tmpdir, log='init-g'+g+'-ft.log', log_append=True)
    s.run(check=True)
    
    # cleanup the tmp dirs
    check_rm(modeldir+'/tmp*')
    
    ###########################################################################################################
    # ADD model group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + DD cal from model, cirular, beam correcred)
    logging.info('Add model...')
    for ms in allmss:
        s.add('taql "update '+ms+' set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log=ms+'_init-taql.log', cmd_type='general')
    s.run(check=True)
    
    # concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:CORRECTED_DATA -> peel_TC*.MS:DATA (empty+DD, avg, phase shifted)
    # concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:MODEL_DATA -> peel-model_TC*.MS:DATA (DD model, avg, phase shifted)
    logging.info('Shifting+averaging (CORRECTED_DATA)...')
    for tc in tcs:
        mss = glob.glob('group*_TC'+tc+'.MS')
        msout = 'peel_TC'+tc+'.MS'
        s.add('NDPPP '+parset_dir+'/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\]', log=msout+'_init-shiftavg.log', cmd_type='NDPPP')
    s.run(check=True)
    
    peelmss = sorted(glob.glob('peel_TC*.MS'))
    
    # Add CORRECTED_DATA for cleaning
    logging.info('Add CORRECTED_DATA...')
    for ms in peelmss:
        s.add('addcol2ms.py -m '+ms+' -c CORRECTED_DATA -i DATA', log=ms+'_init-addcol.log', cmd_type='python', processors='max')
    s.run(check=True)

    # do a first hi-res clean (CORRECTED_DATA is == DATA now)
    model = clean('init', peelmss, dd)
    #model = clean('initfacet', peelmss, dd, facet=True) # DEBUG
   
    ###################################################################################################################
    # self-cal cycle
    for c in xrange(niter):
        logging.info('Start peel cycle: '+str(c))

        # Smooth
        logging.info('BL-based smoothing...')
        for ms in peelmss:
            s.add('BLavg.py -r -w -i DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True)

        if c == 0:
            # make concat after the smoother to have the WEIGHT_SPECTRUM_ORIG included
            logging.info('Concatenating TCs...')
            check_rm('concat.MS*')
            pt.msutil.msconcat(peelmss, 'concat.MS', concatTime=False)
    
        # ft model - peel_TC*.MS:MODEL_DATA (best available model)
        logging.info('FT model...')
        for ms in peelmss:
            s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':ms, 'model':model}, log=ms+'_ft-c'+str(c)+'.log')
            s.run(check=True) # no parallel (problem multiple accesses to model file)
    
        # calibrate phase-only (solve only) - peel_TC*.MS:DATA
        logging.info('Calibrating phase...')
        for ms in peelmss:
            s.add('calibrate-stand-alone -f --parmdb-name instrument_tec '+ms+' '+parset_dir+'/bbs-sol_tec.parset '+skymodel, \
                    log=ms+'_calpreamp-c'+str(c)+'.log', cmd_type='BBS')
        s.run(check=True)
    
        # calibrate phase-only - group*_TC.MS:DATA -> group*_TC.MS:CORRECTED_DATA (selfcal phase corrected, beam corrected)
        logging.info('Correcting phase...')
        for ms in peelmss:
            s.add('calibrate-stand-alone --parmdb-name instrument_tec '+ms+' '+parset_dir+'/bbs-cor_tec.parset '+skymodel, \
            log=ms+'_corpreamp-c'+str(c)+'.log', cmd_type='BBS')
        s.run(check=True)

        # Smooth
        logging.info('Smoothing...')
        for ms in peelmss:
            s.add('BLavg.py -r -w -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth-preamp-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True)

        # calibrate amplitude (solve only) - peel_TC*.MS:SMOOTH_DATA
        logging.info('Calibrating amplitude...')
        for ms in peelmss:
            s.add('calibrate-stand-alone -f --parmdb-name instrument_amp '+ms+' '+parset_dir+'/bbs-sol_amp.parset '+skymodel, \
                    log=ms+'_calamp-c'+str(c)+'.log', cmd_type='BBS', processors = 'max')
        s.run(check=True)
    
        # merge parmdbs
        logging.info('Merging instrument tables...')
        for ms in peelmss:
            merge_parmdb(ms+'/instrument_tec', ms+'/instrument_amp', ms+'/instrument', clobber=True)
    
        # LoSoTo Amp rescaling
        losoto(c, peelmss, dd, parset_dir+'/losoto.parset')
    
        # correct phase + amplitude - peel_TC*.MS:DATA -> peel_TC*.MS:CORRECTED_DATA (selfcal TEC+ph+amp corrected)
        logging.info('Correcting phase+amplitude...')
        for ms in peelmss:
            s.add('calibrate-stand-alone '+ms+' '+parset_dir+'/bbs-cor_amptec.parset '+skymodel, \
                    log=ms+'_coramptec-c'+str(c)+'.log', cmd_type='BBS')
        s.run(check=True)

        logging.info('Restoring WEIGHT_SPECTRUM...')
        s.add('taql "update concat.MS set WEIGHT_SPECTRUM = WEIGHT_SPECTRUM_ORIG"', log='taql-resetweights-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
    
        ######################################################################################################################
        # clean
        model = clean(c, peelmss, dd)
    
    #clean('emptyfacet', peelmss, dd, avgfreq=2, avgtime=5, facet=True, skip_mask=True) # DEBUG

    # backup instrument tables
    logging.info('Back up instrument tables...')
    for ms in peelmss:
        logging.debug('Creating: peel/'+dd['name']+'/instruments/'+ms.replace('MS','parmdb'))
        os.system('cp -r '+ms+'/instrument peel/'+dd['name']+'/instruments/'+ms.replace('MS','parmdb'))
    
    # now do the same but for the entire facet to obtain a complete image of the facet and do a final subtraction
    ##############################################################################################################################
    # Add rest of the facet - group*_TC*.MS:MODEL_DATA (high+low resolution model)
    logging.info('Ft facet model...')
    for g in groups:
        # tmp directory are created to run CASA inside and prevent CASA bug when multiple istances are run in the same dir
        tmpdir = os.getcwd()+'/'+modeldir+'/tmp_'+g
        model = os.getcwd()+'/'+modeldir+'/peel_facet_g'+g+'.model.ms'
        os.makedirs(tmpdir)
        concat = tmpdir+'/concat.MS'
        check_rm(concat+'*')
        pt.msutil.msconcat(sorted(glob.glob('group'+g+'_TC*.MS')), concat, concatTime=False)
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat, 'model':model, 'wproj':512}, wkd=tmpdir, log='facet-g'+g+'-ft.log')
    s.run(check=True)

    logging.info('Ft facet model (lr)...')
    for g in groups:
        tmpdir = os.getcwd()+'/'+modeldir+'/tmp_'+g
        model = os.getcwd()+'/'+modeldir+'/peel_facet_lr_g'+g+'.model.ms'
        concat = tmpdir+'/concat.MS'
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat, 'model':model, 'wproj':512, 'incr':True}, wkd=tmpdir, log='facet-g'+g+'-ft.log', log_append=True)
    s.run(check=True)
    
    # cleanup the tmp dirs
    check_rm(modeldir+'/tmp*')
    
    # ADD group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + facet from model, cirular, beam correcred)
    logging.info('Add facet model...')
    for ms in allmss:
        s.add('taql "update '+ms+' set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log=ms+'_facet-taql.log', cmd_type='general')
    s.run(check=True)
    
    # Concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:CORRECTED_DATA -> facet_TC*.MS:DATA (not corrected, field subtracted but facet, avg, phase shifted)
    logging.info('Shifting+averaging (CORRECTED_DATA)...')
    for tc in tcs:
        mss = glob.glob('group*_TC'+tc+'.MS')
        msout = 'facet_TC'+tc+'.MS'
        s.add('NDPPP '+parset_dir+'/NDPPP-shiftavg.parset msin="['+','.join(mss)+']" msout='+msout+' msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\]', \
                log=msout+'_facet-shiftavg.log', cmd_type='NDPPP')
    s.run(check=True)
    
    facetmss = sorted(glob.glob('facet_TC*.MS'))

    # DEBUG
    for ms in facetmss:
        s.add('addcol2ms.py -m '+ms+' -c CORRECTED_DATA -i DATA', log=ms+'_facet-addcol.log', cmd_type='python', processors='max', log_append=True)
    s.run(check=True)
    clean('precalfacet', facetmss, dd, avgfreq=2, avgtime=5, facet=True, skip_mask=True) # DEBUG
    
    # Correct amp+ph - facet_TC*.MS:DATA -> facet_TC*.MS:CORRECTED_DATA (selfcal phase+amp corrected)
    # Copy instrument table in facet dataset
    for tc in tcs:
        msDD = 'peel_TC'+tc+'.MS'
        msFacet = 'facet_TC'+tc+'.MS'
        logging.debug(msDD+'/instrument -> '+msFacet+'/instrument')
        os.system('cp -r '+msDD+'/instrument '+msFacet+'/instrument')
    logging.info('Correcting facet amplitude+phase...')
    for ms in facetmss:
        s.add('calibrate-stand-alone '+ms+' '+parset_dir+'/bbs-cor_amptec.parset '+skymodel, \
                log=ms+'_facet-coramptec.log', cmd_type='BBS')
    s.run(check=True)
    
    # Cleaning facet
    facetmodel = clean('facet', facetmss, dd, avgfreq=2, avgtime=5, facet=True)

    # Blank pixels outside facet, new foccussed sources are cleaned (so they don't interfere) but we don't want to subtract them
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
            params={'imgs':glob.glob(facetmodel+'*'), 'region':dd['reg_facet'], 'inverse':True}, log='final_casablank.log')
    s.run(check=True)

    ############################################################################################################################
    # in CORRECTED_DATA there's still SUBTRACTED_DATA + (old) MODEL_DATA
    # shift original dataset -  group*_TC*.MS:CORRECTED_DATA -> group*_TC*-shifted.MS:DATA (empty+facet, phase shifted)
    logging.info('Shifting original dataset...')
    for ms in allmss:
        msout = ms.replace('.MS','-shift.MS')
        s.add('NDPPP '+parset_dir+'/NDPPP-shift.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA msout='+msout+' \
                msout.datacolumn=DATA shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\]', \
                log=msout+'_final-shift.log', cmd_type='NDPPP')
    s.run(check=True)
    # copy instrument table from the DD calibration inside the phase shifted full-res MS
    for ms in peelmss:
        for g in groups:
            msout = ms.replace('peel','group'+g).replace('.MS','-shift.MS')
            logging.debug(ms+'/instrument -> '+msout+'/instrument')
            os.system('cp -r '+ms+'/instrument '+msout)

    allmssshifted = sorted(glob.glob('group*_TC*-shift.MS'))

    # add columns that will be used to do ft() in concat mode
    for ms in allmssshifted:
        s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA', log=ms+'_facet-addcol.log', cmd_type='python')
    s.run(check=True)
    
    #################################################################
    # here the new best facet model is subtracted after corruption with DD solution
    # ft model - group*_TC*-shift.MS:MODEL_DATA (best available model)
    logging.info('FT new facet model...')
    for g in groups:
        check_rm('concat.MS*')
        pt.msutil.msconcat(sorted(glob.glob('group'+g+'_TC*-shift.MS')), 'concat.MS', concatTime=False)
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':facetmodel, 'wproj':512}, log='final_ft-g'+g+'.log')
        s.run(check=True) # no parallel (problem multiple accesses to model file)
 
    # SUB corrupted facet model group*_TC*-shift.MS:DATA - MODEL_DATA -> group*_TC*-shift.MS:SUBTRACTED_DATA (empty data + facet from model)
    logging.info('Subtracting new facet model...')
    for ms in allmssshifted:
        s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' '+parset_dir+'/bbs-final_sub.parset', \
                log=ms+'_final-sub.log', cmd_type='BBS')
    s.run(check=True)

    # Shift back dataset -  group*_TC*-shifted.MS:SUBTRACTED_DATA -> group*_TC*-shiftback.MS:DATA (empty, phase shifted)
    logging.info('Shifting back original dataset...')
    for ms in sorted(glob.glob('group*_TC*-shift.MS')):
        msout = ms.replace('-shift.MS','-shiftback.MS')
        s.add('NDPPP '+parset_dir+'/NDPPP-shift.parset msin='+ms+' msin.datacolumn=SUBTRACTED_DATA msout='+msout+' \
                msout.datacolumn=DATA shift.phasecenter=\['+str(phasecentre[0])+'deg,'+str(phasecentre[1])+'deg\]', \
                log=msout+'_final-shift2.log', cmd_type='NDPPP')
    s.run(check=True)

    for ms in sorted(glob.glob('group*_TC*-shiftback.MS')):
        msorig = ms.replace('-shiftback.MS','.MS')
        s.add('taql "update '+msorig+', '+ms+' as shiftback set SUBTRACTED_DATA=shiftback.DATA"', log=ms+'_final-taql.log', cmd_type='general')
    s.run(check=True)

    check_rm('group*_TC*-shift*.MS') # otherwise next wildcard selects them
    
    # Make inspection image 
    logging.info('Inspection image...')
    check_rm('concat.MS*')
    pt.msutil.msconcat(allmss, 'concat.MS', concatTime=False)
    imagename = 'peel/'+dd['name']+'/images/inspection'
    s.add('wsclean_1.8 -datacolumn SUBTRACTED_DATA -reorder -name ' + imagename + ' -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
            -scale 5arcsec -weight briggs 0.0 -niter 1 -no-update-model-required -maxuv-l 8000 -mgain 0.85 concat.MS', \
            log='wsclean-empty.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    
    # backup logs
    os.system('mv logs peel/'+dd['name']+'/')
# end peeling function

for dd in ddset: peel(dd)
logging.info("Done.")
