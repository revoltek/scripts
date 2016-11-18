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
ddset = [{'name': 'src1', 'extended': False, 'facet_extended': False, 'mask':'', 'reg': 'src1.reg', 'reg_facet': 'facet1.reg', 'faint': False, 'coord':[]},
        {'name': 'src2', 'extended': False, 'facet_extended': False, 'mask':'', 'reg': 'src2.reg', 'reg_facet': 'facet2.reg', 'faint': False, 'coord':[]},
        {'name': 'tooth', 'extended': False, 'facet_extended': False, 'mask':'tooth_mask.crtf', 'reg': 'src3.reg', 'reg_facet': 'facet3.reg', 'faint': True, 'coord':[]}]
parset_dir = '/home/fdg/scripts/autocal/LBAsurvey/parset_peel'
niter = 10

##########################################################################################

import sys, os, glob, re
import numpy as np
from lofar import bdsm
import pyrap.tables as pt
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
    check_rm('mss_imgavg')
    os.mkdir('mss_imgavg')
    nchan = find_nchan(mss[0])
    for ms in mss:
        msout = 'mss_imgavg/'+os.path.basename(ms)
        check_rm(msout)
        s.add('NDPPP '+parset_dir+'/NDPPP-avg.parset msin='+ms+' msin.nchan='+str(nchan-nchan%4)+' msin.datacolumn=CORRECTED_DATA \
                msout='+msout+' avg.freqstep='+str(avgfreq)+' avg.timestep='+str(avgtime), log=ms+'_cleanavg-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)
    mssavg = [ms for ms in glob.glob('mss_imgavg/*MS')]

    # set image name
    imagename = 'peel/'+dd['name']+'/images/peel-'+str(c)

    # set pixscale, imsize and niter
    pixscale = scale_from_ms(mss[0])
    if facet:
        imsize = size_from_reg('peel/'+dd['name']+'/models/peel_facet-0000-model.fits', 'regions/'+dd['reg_facet'], dd['coord'], pixscale)
#        niter = 10000
    else:
        imsize = size_from_reg('peel/'+dd['name']+'/models/peel_dd-0000-model.fits', 'regions/'+dd['reg'], dd['coord'], pixscale)
#        niter = 2000

    if imsize < 512:
        trim = 512
        imsize = 1024
    elif imsize < 1024:
        trim = 1024
        imsize = 1024
    else:
        trim = imsize

    logging.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    # TODO: add multiscale
    if dd['extended']: multiscale = [0,4,16,64]
    else: multiscale = []

    # Clean mask clean
    logging.info('Cleaning (cycle: '+str(c)+')...')
    s.add('wsclean -reorder -name ' + imagename + ' -size '+str(imsize)+' '+str(imsize)+' -trim '+str(trim)+' '+str(trim)+' \
            -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale '+str(pixscale)+'arcsec -weight briggs 0.0 -threshold 0.01 -niter 10000 -no-update-model-required -mgain 0.6 \
            -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 '+' '.join(mss), \
            log='wsclean-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    if skip_mask:
        check_rm('mss_imgavg')
        return imagename, trim, pixscale

    if dd['faint'] or facet:
        make_mask(image_name = imagename+'-MFS-image.fits', mask_name = imagename+'.newmask', atrous_do=dd['extended'], threshisl=5)
    else:
        make_mask(image_name = imagename+'-MFS-image.fits', mask_name = imagename+'.newmask', atrous_do=dd['extended'], threshisl=20)

    # if dd['mask'] is set then add it to the new mask
    if dd['mask'] != '':
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
            params={'imgs':imagename+'.newmask', 'region':dd['mask'], 'setTo':1}, log='casablank-c'+str(c)+'.log', log_append=True)
        s.run(check=True)

    logging.info('Cleaning with mask (cycle: '+str(c)+')...')
    s.add('wsclean -reorder -name ' + imagename + '-M -size '+str(imsize)+' '+str(imsize)+' -trim '+str(trim)+' '+str(trim)+' \
            -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 -casamask '+imagename+'.newmask \
            -scale '+str(pixscale)+'arcsec -weight briggs 0.0 -threshold 0.01 -niter 5000 -no-update-model-required -mgain 0.6 \
            -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 '+' '.join(mss), \
            log='wscleanM-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    check_rm('mss_imgavg')
    return imagename, trim, pixscale


def losoto(c, mss, dd, parset):
    logging.info('Running LoSoTo...')
    check_rm('plots')
    os.makedirs('plots')
    check_rm('globaldb')
    os.makedirs('globaldb')

    for num, ms in enumerate(mss):
        os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
        if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD globaldb/')

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
    check_rm('mss_peel') 
    check_rm('mss_shift')
    check_rm('mss_facet') 
    check_rm('plot')
    
    logging.info('Creating dirs...')
    os.makedirs('logs/mss')
    os.makedirs('logs/mss_peel')
    os.makedirs('logs/mss_shift')
    os.makedirs('logs/mss_facet')
    os.makedirs('mss_peel')
    os.makedirs('mss_shift')
    os.makedirs('mss_facet')
    check_rm('peel/'+dd['name'])
    os.makedirs('peel/'+dd['name'])
    os.makedirs('peel/'+dd['name']+'/models')
    os.makedirs('peel/'+dd['name']+'/images')
    os.makedirs('peel/'+dd['name']+'/instruments')
    os.makedirs('peel/'+dd['name']+'/plots')
    os.makedirs('peel/'+dd['name']+'/h5')
    
    logging.info('Indexing...')
    allmss = sorted(glob.glob('mss/TC*.MS'))
    phasecentre = get_phase_centre(allmss[0])
    
    tcs = []
    for ms in allmss:
        tc = re.findall(r'\d+', ms)[0] # time chunk number
        tcs.append(tc)
    tcs = list(set(tcs))

    centroid_ra, centroid_dec = get_coord_centroid(glob.glob('self/models/*.fits')[0], 'regions/'+dd['reg'])
    dd['coord']=[centroid_ra,centroid_dec]

    # preparing concatenated dataset
    concat_ms = 'mss/concat.MS'
    check_rm(concat_ms+'*')
    pt.msutil.msconcat(allmss, 'mss/concat.MS', concatTime=False)
    
    #################################################################################################
    # Blank unwanted part of models
    modeldir = 'peel/'+dd['name']+'/models/'
    
    logging.info('Splitting skymodels...')
    for model in sorted(glob.glob('self/models/*.fits')):
        logging.debug(model)
        outfile = modeldir+'/'+os.path.basename(model).replace('coadd','peel_dd')
        blank_image(model, 'regions/'+dd['reg'], outfile, inverse=True)
        outfile = modeldir+'/'+os.path.basename(model).replace('coadd','peel_facet')
        blank_image(model, 'regions/'+dd['reg_facet'], outfile, inverse=True)

    #####################################################################################################
    # Add DD cal model - group*_TC*.MS:MODEL_DATA (high+low resolution model)
    logging.info('Ft DD calibrator model...')
    s.add('wsclean -predict -name ' + modeldir + 'peel_dd -size 8000 8000 -mem 90 -j '+str(s.max_processors)+' \
            -scale 10arcsec -channelsout 10 '+concat_ms, \
            log='wscleanPRE.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    ###########################################################################################################
    # ADD model group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + DD cal from model, cirular, beam correcred)
    logging.info('Add model...')
    for ms in allmss:
        s.add('taql "update '+ms+' set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log=ms+'_init-taql.log', cmd_type='general')
    s.run(check=True)
    
    # concat all groups (freq) + avg (to 1 chan/SB, 5 sec) -  group*_TC*.MS:CORRECTED_DATA -> peel_TC*.MS:DATA (empty+DD, avg, phase shifted)
    logging.info('Shifting+averaging (CORRECTED_DATA)...')
    for ms in allmss:
        msout = ms.replace('mss','mss_peel')
        s.add('NDPPP '+parset_dir+'/NDPPP-shiftavg.parset msin='+ms+' msout='+msout+' msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=\['+str(dd['coord'][0])+'deg,'+str(dd['coord'][1])+'deg\]', log=msout+'_init-shiftavg.log', cmd_type='NDPPP')
    s.run(check=True)
    
    peelmss = sorted(glob.glob('mss_peel/TC*.MS'))
    
    # Add MODEL_DATA and CORRECTED_DATA for cleaning
    logging.info('Add MODEL_DATA and CORRECTED_DATA...')
    for ms in peelmss:
        s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA,CORRECTED_DATA -i DATA', log=ms+'_init-addcol.log', cmd_type='python', processors='max')
    s.run(check=True, max_threads=1)

    # do a first hi-res clean (CORRECTED_DATA is == DATA now)
    # TODO: use available model? it's a mess with phase centers: needs reprojecting
    model, imsize, pixscale = clean('init', peelmss, dd)
    #clean('initfacet', peelmss, dd, avgfreq=2, avgtime=5, facet=True, skip_mask=True) # DEBUG
   
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
            concat_ms_peel = 'mss_peel/concat.MS'
            check_rm(concat_ms_peel+'*')
            pt.msutil.msconcat(peelmss, concat_ms_peel, concatTime=False)
    
        # ft model - peel_TC*.MS:MODEL_DATA (best available model)
        logging.info('FT model...')
        s.add('wsclean -predict -name ' + model + ' -size '+str(imsize)+' '+str(imsize)+' -mem 90 -j '+str(s.max_processors)+' \
                -scale '+str(pixscale)+'asec -channelsout 10 '+concat_ms_peel, \
                log='wscleanPRE-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
    
        # solve+correct TEC - group*_TC.MS:SMOOTHED_DATA -> group*_TC.MS:CORRECTED_DATA
        logging.info('Solving TEC...')
        for ms in peelmss:
            check_rm(ms+'/instrument-tec')
            s.add('NDPPP '+parset_dir+'/NDPPP-solTEC.parset msin='+ms+' cal.parmdb='+ms+'/instrument-tec', \
                log=ms+'_sol-tec-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)
        logging.info('Correcting TEC...')
        for ms in peelmss:
            s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' cor1.parmdb='+ms+'/instrument-tec', \
                log=ms+'_cor-tec-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # calibrate amplitude (solve only) - peel_TC*.MS:CORRECTED_DATA
        logging.info('Calibrating amplitude...')
        for ms in peelmss:
            check_rm(ms+'/instrument-amp')
            s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' cal.parmdb='+ms+'/instrument-amp cal.solint=60', \
                log=ms+'_sol-g-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)
    
        # merge parmdbs
        logging.info('Merging instrument tables...')
        for ms in peelmss:
            merge_parmdb(ms+'/instrument-tec', ms+'/instrument-amp', ms+'/instrument', clobber=True)
    
        # LoSoTo Amp rescaling + plotting
        losoto(c, peelmss, dd, parset_dir+'/losoto.parset')
        sys.exit(1)
    
        # correct TEC+amplitude - peel_TC*.MS:DATA -> peel_TC*.MS:CORRECTED_DATA
        logging.info('Correcting phase+amplitude...')
        for ms in peelmss:
            s.add('NDPPP '+parset_dir+'/NDPPP-corTECG.parset msin='+ms+' cor1.parmdb='+ms+'/instrument cor2.parmdb='+ms+'/instrument cor3.parmdb='+ms+'/instrument', \
                log=ms+'_cor-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        logging.info('Restoring WEIGHT_SPECTRUM...')
        s.add('taql "update '+concat_ms_peel+' set WEIGHT_SPECTRUM = WEIGHT_SPECTRUM_ORIG"', log='taql-resetweights-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
    
        ######################################################################################################################
        # clean
        model, imsize, pixscale = clean(c, peelmss, dd)
    
    clean('emptyfacet', peelmss, dd, avgfreq=2, avgtime=5, facet=True, skip_mask=True) # DEBUG

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
        msout = ms.replace('mss','mss_shift')
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

    allmssshifted = sorted(glob.glob('mss_shift/TC*.MS'))

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
    s.add('wsclean_1.8 -datacolumn SUBTRACTED_DATA -reorder -name ' + imagename + ' -size 5000 5000 -mem 90 -j '+str(s.max_processors)+' \
            -scale 5arcsec -weight briggs 0.0 -niter 1 -no-update-model-required -maxuv-l 8000 -mgain 0.85 concat.MS', \
            log='wsclean-empty.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    
    # backup logs
    os.system('mv logs peel/'+dd['name']+'/')
# end peeling function

for dd in ddset: peel(dd)
logging.info("Done.")
