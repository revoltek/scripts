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

parset_dir = '/home/fdg/scripts/autocal/LBAsurvey/parset_peel'
maxniter = 10 # max iteration if source not converged
pb_cut = 5 # degree to cut faceting

##########################################################################################

import sys, os, glob, re
import numpy as np
from autocal.lib_pipeline import *
from autocal.lib_pipeline_dd import *
from autocal.lib_pipeline_img import *
from make_mask import make_mask
from lofar import bdsm
import pyrap.tables as pt

set_logger('pipeline-peel.logging')
check_rm('logs')
s = Scheduler(dry=False)
allmss = sorted(glob.glob('mss/TC*.MS'))
phasecentre = get_phase_centre(allmss[0])

def clean(c, mss, dd, avgfreq=8, avgtime=10, facet=False):
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
    mssavg = [ms for ms in sorted(glob.glob('mss_imgavg/*MS'))]

    # set pixscale and imsize
    pixscale = scale_from_ms(mss[0])
    if facet:
        imsize = int((dd['facet_size']/(pixscale/3600.))*1.5)
    else:
        imsize = int((dd['dd_size']/(pixscale/3600.))*1.5)

    if imsize < 512:
        imsize = 512

    trim = int(imsize*0.7)

    if imsize % 2 == 1: imsize += 1 # make even
    if trim % 2 == 1: trim += 1 # make even

    logging.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    # -trim '+str(trim)+' '+str(trim)+'
    # -auto-mask 5 -auto-threshold 1 -rms-background -rms-background-window 25 \
    # -multiscale

    # clean 1
    logging.info('Cleaning (cycle: '+str(c)+')...')
    if facet: imagename = 'img/facet-'+str(c)
    else: imagename = 'img/ddcal-'+str(c)
    s.add('/home/fdg/opt/src/wsclean-2.2.9/build/wsclean -reorder -name ' + imagename + ' -size '+str(imsize)+' '+str(imsize)+' -trim '+str(trim)+' '+str(trim)+' \
            -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale '+str(pixscale)+'arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -mgain 0.7 -pol I \
            -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 \
            -auto-threshold 20 '+' '.join(mss), \
            log='wsclean-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    os.system('cat logs/wsclean-c'+str(c)+'.log | grep Jy')

    # make mask
    maskname = imagename+'-mask.fits'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 3)

    # clean 2
    logging.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    if facet: imagename = 'img/facetM-'+str(c)
    else: imagename = 'img/ddcalM-'+str(c)
    s.add('/home/fdg/opt/src/wsclean-2.2.9/build/wsclean -reorder -name ' + imagename + ' -size '+str(imsize)+' '+str(imsize)+' -trim '+str(trim)+' '+str(trim)+' \
            -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale '+str(pixscale)+'arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -mgain 0.7 -pol I \
            -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 \
            -auto-threshold 1 -fitsmask '+maskname+' '+' '.join(mss), \
            log='wscleanM-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    os.system('cat logs/wscleanM-c'+str(c)+'.log | grep Jy')

    # remove CC not in mask
    maskname = imagename+'-mask.fits'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 5)
    for modelname in sorted(glob.glob(imagename+'*model.fits')):
        blank_image_fits(modelname, maskname, inverse=True)

    check_rm('mss_imgavg')
    return imagename


def losoto(c, mss, dd, parset, instrument='instrument', putback=True):
    """
    c : name for plot dir
    mss : input set of mss
    dd : direction table
    parset : losoto parset
    instrument : parmdb name
    putback : if False skip the exporter
    """
    logging.info('Running LoSoTo...')
    check_rm('plots')
    os.makedirs('plots')
    check_rm('globaldb')
    os.makedirs('globaldb')

    for num, ms in enumerate(mss):
        os.system('cp -r '+ms+'/'+instrument+' globaldb/instrument-'+str(num))
        if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD globaldb/')

    h5parm = 'global-c'+str(c)+'.h5'

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('losoto -v '+h5parm+' '+parset, log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    if putback:
        s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
        s.run(check=True)

        for num, ms in enumerate(mss):
            check_rm(ms+'/'+instrument)
            os.system('mv globaldb/sol000_instrument-'+str(num)+' '+ms+'/'+instrument)

    os.system('mv plots peel/'+dd['name']+'/plots/plots-c'+str(c))
    os.system('mv '+h5parm+' peel/'+dd['name']+'/h5')


def peel(dd):

    logging.info('######################## Peeling: '+dd['name']+' #################################')

    with open('pipeline-peel.status', 'r') as status_file:
        if any(dd['name'] == line.rstrip('\r\n') for line in status_file):
            logging.warning('Direction %s already done: skip.' % dd['name'])
            return

    #########################################################################################
    # Clear
    logging.info('Cleaning...')
    check_rm('mss_peel') 
    check_rm('mss_facet') 
    check_rm('mss_shiftback') 
    check_rm('plot')
    check_rm('img')
    
    logging.info('Creating dirs...')
    os.makedirs('logs/mss')
    os.makedirs('logs/mss_peel')
    os.makedirs('logs/mss_facet')
    os.makedirs('logs/mss_shiftback')
    os.makedirs('mss_peel')
    os.makedirs('mss_shiftback')
    check_rm('peel/'+dd['name'])
    os.makedirs('peel/'+dd['name'])
    os.makedirs('peel/'+dd['name']+'/models')
    os.makedirs('peel/'+dd['name']+'/images')
    os.makedirs('peel/'+dd['name']+'/masks')
    os.makedirs('peel/'+dd['name']+'/plots')
    os.makedirs('peel/'+dd['name']+'/h5')
    os.makedirs('img')
    
    logging.info('Indexing...')
    modeldir = 'peel/'+dd['name']+'/models/'
   
    #clean('emptybefore', allmss, dd, avgfreq=1, avgtime=5, facet=True) # DEBUG

    #################################################################################################
    # Blank unwanted part of models + intersect with beam
    logging.info('Splitting skymodels...')
    for model in sorted(glob.glob('self/models/*.fits')):
        logging.debug(model)
        outfile = modeldir+'/'+os.path.basename(model).replace('coadd','large_peel_dd')
        blank_image_reg(model, ['regions/'+dd['name']+'.reg', 'regions/beam.reg'], outfile, inverse=True, op='AND')
        #blank_image_reg(outfile, 'regions/beam.reg', outfile, inverse=True)
        if dd['facet_size'] > 0:
            outfile = modeldir+'/'+os.path.basename(model).replace('coadd','large_peel_facet')
            blank_image_reg(model, ['regions/'+dd['name']+'-facet.reg', 'regions/beam.reg'], outfile, inverse=True, op='AND')
            #blank_image_reg(outfile, 'regions/beam.reg', outfile, inverse=True)

    ##############################################################
    # reproject + cut model image to speed up prediction
    logging.info("Reprojecting models...")
    s.add('mHdr -p 10 "%f %f" %f %s/dd.hdr' % (dd['RA'], dd['DEC'], dd['dd_size']*5, modeldir), log='reproject.log', log_append=True, cmd_type="general")
    if dd['facet_size'] > 0:
        s.add('mHdr -p 10 "%f %f" %f %s/facet.hdr' % (dd['RA'], dd['DEC'], dd['facet_size']*1.2, modeldir), log='reproject.log', log_append=True, cmd_type="general")
    s.run(check=True)
    os.system('sed -i \'s/TAN/SIN/\' '+modeldir+'/*.hdr') # wsclean wants SIN projection
    for model in sorted(glob.glob(modeldir+'large_peel_dd*.fits')):
        outmodel = model.replace('large_','')
        s.add("mProjectPP "+model+" "+outmodel+" "+modeldir+"dd.hdr", log='reproject.log', log_append=True, cmd_type="general")
    if dd['facet_size'] > 0:
        for model in sorted(glob.glob(modeldir+'large_peel_facet*.fits')):
            outmodel = model.replace('large_','')
            s.add("mProjectPP "+model+" "+outmodel+" "+modeldir+"facet.hdr", log='reproject.log', log_append=True, cmd_type="general")
    s.run(check=True)
    check_rm(modeldir+'*hdr')
    check_rm(modeldir+'*area.fits')
    check_rm(modeldir+'large*')
    # remove NaNs that mProject can create
    for model in glob.glob(modeldir+'/*fits'):
        nan2zeros(model)

    ###########################################################
    # keep SUBTRACTED_DATA as a working columns so we can re-start each time
    logging.info('Add SUBTRACTED_DATA (only first direction)...')
    for ms in allmss:
        s.add('addcol2ms.py -m '+ms+' -c SUBTRACTED_DATA -i CORRECTED_DATA', log=ms+'_init-addcol.log', cmd_type='python', processors='max', log_append=True)
    s.run(check=True)
 
    ###################################################################
    # avg and ph-shift (to 4 chan/SB, 4 sec) -  mss/TC*.MS:SUBTRACTED_DATA -> mss_peel/TC*.MS:DATA (empty+DD, phase shifted)
    logging.info('Shifting (SUBTRACTED_DATA)...')
    for ms in allmss:
        msout = ms.replace('mss','mss_peel')
        s.add('NDPPP '+parset_dir+'/NDPPP-shift.parset msin='+ms+' msout='+msout+' msin.datacolumn=SUBTRACTED_DATA \
                shift.phasecenter=\['+str(dd['RA'])+'deg,'+str(dd['DEC'])+'deg\]', log=msout+'_init-shift.log', cmd_type='NDPPP')
    s.run(check=True)

    peelmss = sorted(glob.glob('mss_peel/TC*.MS'))

    #####################################################################################################
    # BKP empty DATA for faceting
    logging.info('Add EMTPY_DATA...')
    for ms in peelmss:
        s.add('addcol2ms.py -m '+ms+' -c EMPTY_DATA -i DATA', log=ms+'_init-addcol.log', cmd_type='python', processors='max', log_append=True)
    s.run(check=True)
 
    ######################################################################################################
    # Add DD cal model - peel_mss/TC*.MS:MODEL_DATA
    logging.info('Add MODEL_DATA...')
    for ms in peelmss:
        s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA', log=ms+'_init-addcol.log', cmd_type='python', processors='max', log_append=True)
    s.run(check=True)
    logging.info('Ft DD calibrator model...')
    s.add('wsclean -predict -name ' + modeldir + 'peel_dd -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+' '.join(peelmss), \
            log='wscleanPRE-dd.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    # ADD model peel_mss/TC*.MS:DATA + MODEL_DATA -> peel_mss/TC*.MS:DATA (empty data + DD cal from model)
    logging.info('Set DATA = DATA + MODEL_DATA...')
    for ms in peelmss:
        s.add('taql "update '+ms+' set DATA = DATA + MODEL_DATA"', log=ms+'_init-taql.log', cmd_type='general')
    s.run(check=True)

    ###########################################################################################################
    # Add CORRECTED_DATA for cleaning
    logging.info('Add CORRECTED_DATA...')
    for ms in peelmss:
        s.add('addcol2ms.py -m '+ms+' -c CORRECTED_DATA -i DATA', log=ms+'_init-addcol.log', cmd_type='python', processors='max', log_append=True)
    s.run(check=True)
    # do a first clean to get the starting model
    model = clean('init', peelmss, dd)
    rms_noise = get_noise_img(model+'-MFS-residual.fits')

    ###################################################################################################################
    # self-cal cycle
    rms_noise_pre = np.inf
    for c in xrange(maxniter):
        logging.info('Start peel cycle: '+str(c))

        # ft model - mss_peel/TC*.MS:MODEL_DATA (best available model)
        logging.info('FT model...')
        s.add('wsclean -predict -name ' + model + ' -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+' '.join(peelmss), \
                log='wscleanPRE-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
 
        # Smooth peel_mss/TC*.MS:DATA -> peel_mss/TC*.MS:CORRECTED_DATA (smoothed data)
        logging.info('BL-based smoothing...')
        for ms in peelmss:
            s.add('BLsmooth.py -r -i DATA -o CORRECTED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True)
   
        ####################################
        # solve+correct TEC - mss_peel/TC*.MS:CORRECTED_DATA (smoothed!) -> mss_peel/TC*.MS:CORRECTED_DATA (smoothed corrected TEC)

        # sol.solint depends on peak flux
        if dd['Peak_flux'] > 5: solint = 1
        elif dd['Peak_flux'] > 4: solint = 2
        elif dd['Peak_flux'] > 3: solint = 4
        elif dd['Peak_flux'] > 2: solint = 8
        else: solint = 16

        logging.info('Solving TEC...')
        for ms in peelmss:
            check_rm(ms+'/instrument-tec')
            s.add('NDPPP '+parset_dir+'/NDPPP-solTEC.parset msin='+ms+' sol.parmdb='+ms+'/instrument-tec sol.solint='+str(solint), \
                log=ms+'_sol-tec-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)
        
        losoto(str(c)+'-tec', peelmss, dd, parset_dir+'/losoto-tec.parset', instrument='instrument-tec', putback=False)

        # correct on smoothed data only when solve also amp
        if c > 0 and dd['Total_flux'] > 3: incol = 'CORRECTED_DATA' # <- smoothed
        else: incol = 'DATA'

        logging.info('Correcting TEC...')
        for ms in peelmss:
            s.add('NDPPP '+parset_dir+'/NDPPP-corTEC.parset msin='+ms+' msin.datacolumn='+incol+' cor1.parmdb='+ms+'/instrument-tec cor2.parmdb='+ms+'/instrument-tec', \
                log=ms+'_cor-tec-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        #######################################
        # calibrate amplitude (solve only) - mss_peel/TC*.MS:CORRECTED_DATA
        if c > 0 and dd['Peak_flux'] > 3:

            # TODO: if strong sources add multichan amp solve, wait for NDPPP to divide chans evenly
            #nchan = find_nchan(peelmss[0])
            #nchan = min(np.rint(nchan/10.), np.rint(nchan/dd['Total_flux']/10.))

            logging.info('Calibrating amplitude...')
            for ms in peelmss:
                check_rm(ms+'/instrument-amp')
                s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' sol.parmdb='+ms+'/instrument-amp sol.solint=100 sol.nchan=0', \
                    log=ms+'_sol-g-c'+str(c)+'.log', cmd_type='NDPPP')
            s.run(check=True)
    
            # LoSoTo Amp rescaling + plotting
            losoto(str(c)+'-amp', peelmss, dd, parset_dir+'/losoto-g.parset', instrument='instrument-amp')
    
            # correct TEC+amplitude - mss_peel/TC*.MS:DATA -> mss_peel/TC*.MS:CORRECTED_DATA (corrected TEC+G)
            logging.info('Correcting phase+amplitude...')
            for ms in peelmss:
                s.add('NDPPP '+parset_dir+'/NDPPP-corTECG.parset msin='+ms+' cor1.parmdb='+ms+'/instrument-tec cor2.parmdb='+ms+'/instrument-tec cor3.parmdb='+ms+'/instrument-amp', \
                    log=ms+'_cor-tecg-c'+str(c)+'.log', cmd_type='NDPPP')
            s.run(check=True)

        # clean
        model = clean(c, peelmss, dd)

        # get noise, if larger than 95% of prev cycle: break
        rms_noise = get_noise_img(model+'-MFS-residual.fits')
        if rms_noise > 0.95 * rms_noise_pre: break
        rms_noise_pre = rms_noise

    if dd['facet_size'] > 0:
        logging.info('Doing facet...')

        # now do the same but for the entire facet to obtain a complete image of the facet and do a final subtraction
        ##############################################################################################################################
        # Cannot avg since the same dataset has to be shifted back and used for other facets
        # Phase shift -  mss/TC*.MS:SUBTRACTED_DATA -> mss_facet/TC*.MS:DATA (not corrected, field subtracted but facet, phase shifted)
        #logging.info('Shifting (SUBTRACTED_DATA)...')
        #for ms in allmss:
        #    msout = ms.replace('mss','mss_facet')
        #    s.add('NDPPP '+parset_dir+'/NDPPP-shift.parset msin='+ms+' msout='+msout+' msin.datacolumn=SUBTRACTED_DATA \
        #            shift.phasecenter=\['+str(dd['RA'])+'deg,'+str(dd['DEC'])+'deg\]', log=msout+'_facet-shift.log', cmd_type='NDPPP')
        #s.run(check=True)
        #
        #facetmss = sorted(glob.glob('mss_facet/TC*.MS'))

        # mss_peel/TC*.MS:DATA = EMPTY_DATA (empty data)
        logging.info('Copy back EMPTY_DATA...')
        for ms in peelmss:
            s.add('taql "update '+ms+' set DATA = EMPTY_DATA"', log=ms+'_facet-taql.log', cmd_type='general')
        s.run(check=True)

        # Add rest of the facet - mss_peel/TC*.MS:MODEL_DATA (high+low resolution facet model)
        logging.info('Ft facet model...')
        s.add('wsclean -predict -name ' + modeldir + 'peel_facet -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+' '.join(peelmss), \
                log='wscleanPRE-facet1.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
        # ADD mss_peel/TC*.MS:DATA + MODEL_DATA -> mss_facet/TC*.MS:DATA (empty data + facet from model)
        logging.info('Add facet model...')
        for ms in peelmss:
            s.add('taql "update '+ms+' set DATA = DATA + MODEL_DATA"', log=ms+'_facet-taql.log', cmd_type='general', log_append=True)
        s.run(check=True)
    
        ### DEBUG
        #logging.info('Set CORRECTED_DATA = DATA')
        #for ms in facetmss:
        #    s.add('addcol2ms.py -m '+ms+' -c CORRECTED_DATA -i DATA', log=ms+'_facet-addcolDEBUG.log', cmd_type='python', processors='max', log_append=True)
        #s.run(check=True)
        #clean('initfacet', facetmss, dd, avgfreq=4, avgtime=5, facet=True) # DEBUG
        
        # Correct amp+ph - mss_facet/TC*.MS:DATA -> mss_facet/TC*.MS:CORRECTED_DATA (selfcal tec+amp corrected)
        # Copy instrument table in facet dataset
        #for msFacet in facetmss:
        #    msDD = msFacet.replace('mss_facet','mss_peel')
        #    logging.debug(msDD+'/instrument -> '+msFacet+'/instrument')
        #    check_rm(msFacet+'/instrument')
        #    os.system('cp -r '+msDD+'/instrument '+msFacet+'/instrument')
        logging.info('Correcting facet amplitude+phase...')
        for ms in peelmss:
            s.add('NDPPP '+parset_dir+'/NDPPP-corTECG.parset msin='+ms+' \
                   cor1.parmdb='+ms+'/instrument-tec cor2.parmdb='+ms+'/instrument-tec cor3.parmdb='+ms+'/instrument-amp', \
                log=ms+'_facet-cor.log', cmd_type='NDPPP')
        s.run(check=True)
        
        # Cleaning facet
        model = clean('facet', peelmss, dd, avgtime=5, facet=True)
        # DEBUG
        #model = 'peel/src1/images/peel-facet'
       
        logging.info('Add MODEL_DATA')
        for ms in peelmss:
            s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA', log=ms+'_facet-addcol.log', cmd_type='python', processors='max', log_append=True)
        s.run(check=True)
    
        # Blank pixels outside facet, new foccussed sources are cleaned (so they don't interfere) but we don't want to subtract them
        logging.info('Blank pixels outside facet...')
        for modelfits in glob.glob(model+'*model.fits'):
            blank_image_reg(modelfits, 'regions/'+dd['name']+'-facet.reg', inverse=True, blankval=0.)
 
    # for ddcal without associated facet
    else:
        logging.info('This DD cal does not have an associate facet, just subtract it...')

        # Blank pixels outside facet, new foccussed sources are cleaned (so they don't interfere) but we don't want to subtract them
        logging.info('Blank pixels outside dd region...')
        for modelfits in glob.glob(model+'*model.fits'):
            blank_image_reg(modelfits, 'regions/'+dd['name']+'.reg', inverse=True, blankval=0.)

    # ft model - mss_peel/TC*.MS:MODEL_DATA (best available model)
    logging.info('FT model...')
    s.add('wsclean -predict -name ' + model + ' -mem 90 -j '+str(s.max_processors)+' -channelsout 10 '+' '.join(peelmss), \
                log='wscleanPRE-facet2.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    for ms in peelmss:
        s.add('taql "update '+ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log=ms+'_facet-taql.log', cmd_type='general', log_append=True)
    s.run(check=True)

    # Corrupt empty data amp+ph - mss_peel/TC*.MS:CORRECTED_DATA -> mss_peel/TC*.MS:CORRECTED_DATA (selfcal empty)
    # TODO: do not corrupt on first calibrator to propagate its solutions
    logging.info('Corrupting facet amplitude+phase...')
    for ms in peelmss:
        s.add('NDPPP '+parset_dir+'/NDPPP-corTECG.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA \
                cor1.parmdb='+ms+'/instrument-tec cor1.invert=false cor2.parmdb='+ms+'/instrument-tec cor2.invert=false cor3.parmdb='+ms+'/instrument-amp cor3.invert=false', \
                log=ms+'_facet-corrupt.log', cmd_type='NDPPP')
    s.run(check=True)

    logging.info('Shifting back...')
    for ms in peelmss:
        msout = ms.replace('mss_peel','mss_shiftback')
        s.add('NDPPP '+parset_dir+'/NDPPP-shift.parset msin='+ms+' msout='+msout+' msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=\['+str(phasecentre[0])+'deg,'+str(phasecentre[1])+'deg\]', log=msout+'_facet-shiftback.log', cmd_type='NDPPP')
    s.run(check=True)

    shiftbackmss = sorted(glob.glob('mss_shiftback/TC*.MS'))
    for ms in shiftbackmss:
        msorig = ms.replace('mss_shiftback','mss')
        s.add('taql "update '+msorig+', '+ms+' as shiftback set SUBTRACTED_DATA=shiftback.DATA"', log=ms+'_taql.log', cmd_type='general')
    s.run(check=True)

    logging.info('Add SUBTRACTED_DATA_DDCAL## (backup)...')
    for ms in allmss:
        colout = 'SUBTRACTED_DATA_'+dd['name'].upper()
        s.add('addcol2ms.py -m '+ms+' -c '+colout+' -i SUBTRACTED_DATA', log=ms+'_bkp-addcol.log', cmd_type='python', processors='max')
    s.run(check=True)

    #clean('emptyafter', allmss, dd, avgfreq=4, avgtime=5, facet=True) # DEBUG

    # backup logs
    os.system('mv logs peel/'+dd['name']+'/')
    # save images
    os.system('mv img/*MFS-image.fits peel/'+dd['name']+'/images/')
    # save masks
    os.system('mv img/*mask.fits peel/'+dd['name']+'/masks/')
    # direction completed, write status
    with open('pipeline-peel.status', 'a') as status_file:
        status_file.write(dd['name']+'\n')

# end peeling function

# Run pyBDSM to create a model used to find good DD-calibrator and tassellate the sky
logging.info('Finding directions...')
imagename = sorted(glob.glob('self/images/wide-[0-9]-MFS-image.fits'))[-1]
if not os.path.exists('regions/DIEcatalog.fits'):
    bdsm_img = bdsm.process_image(imagename, rms_box=(55,12), \
        thresh_pix=5, thresh_isl=3, atrous_do=False, atrous_jmax=3, \
        adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(80,20), \
        quiet=True)
    check_rm('regions')
    os.makedirs('regions')
    bdsm_img.write_catalog(outfile='regions/DIEcatalog.fits', catalog_type='srl', format='fits')

ddset = make_directions_from_skymodel('regions/DIEcatalog.fits', outdir='regions', flux_min_Jy=1.0, size_max_arcmin=3.0,
        directions_separation_max_arcmin=5.0, directions_max_num=20, flux_min_for_merging_Jy=0.2)

logging.info('Voronoi tassellation...')
make_beam_reg(phasecentre[0], phasecentre[1], pb_cut, 'regions/beam.reg')
ddset = make_tassellation(ddset, imagename, beam_reg='regions/beam.reg')

print ddset
ddset.write('ddset.txt', format='ascii')

if not os.path.exists("pipeline-peel.status"):
    os.mknod("pipeline-peel.status")

for dd in ddset: peel(dd)
logging.info("Done.")
