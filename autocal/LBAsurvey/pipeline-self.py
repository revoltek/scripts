#!/usr/bin/python
# perform self-calibration on a group of SBs concatenated in TCs. Script must be run in dir with MS.
# number/chan in MS are flexible but the must be concatenable (same chans/freq!)
# Input:
# TCs are blocks of SBs should have calibrator corrected (a+p) data in DATA (beam not applied).
# file format of TCs is: group#_TC###.MS.
# Output:
# TCs with selfcal NOT corrected source subtracted data in SUBTRACTED_DATA
# TCs with selfcal corrected source subtracted data in CORRECTED_DATA
# instrument tables contain gain (slow) + fast (scalarphase+TEC) solutions
# last high/low resolution models are copied in the "self/models" dir
# last high/low resolution images + masks + empty images (CORRECTED_DATA) are copied in the "self/images" dir
# h5parm solutions and plots are copied in the "self/solutions" dir

parset_dir = '/home/fdg/scripts/autocal/LBAsurvey/parset_self/'
skymodel = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.GMRT150_field.skymodel'
sourcedb = '/home/fdg/scripts/autocal/LBAsurvey/toothbrush.GMRT150_field.skydb'
niter = 2

#######################################################################################

import sys, os, glob, re
import numpy as np
from lofar import bdsm
import pyrap.tables as pt
from lib_pipeline import *
from make_mask import make_mask

set_logger()
check_rm('logs')
s = Scheduler(dry=False)

def losoto(c, mss, g, parset, instrument_in='instrument', instrument_out='instrument'):
    """
    c = cycle
    mss = list of mss
    g = group number
    parset = losoto parset
    """
    logging.info('Running LoSoTo...')
    check_rm('plots')
    os.makedirs('plots')
    check_rm('globaldb')
    os.makedirs('globaldb')

    for num, ms in enumerate(mss):
        if num == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
        os.system('cp -r '+ms+'/'+instrument_in+' globaldb/'+instrument_in+'-'+str(num))

    h5parm = 'global-c'+str(c)+'.h5'
    check_rm(h5parm)

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('losoto -v '+h5parm+' '+parset, log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -t scalaramplitude000 -c '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    for num, ms in enumerate(mss):
        check_rm(ms+'/'+instrument_out)
        os.system('mv globaldb/sol000_'+instrument_in+'-'+str(num)+' '+ms+'/'+instrument_out)
    os.system('mv plots self/solutions/g'+g+'/plots-c'+str(c))
    os.system('mv '+h5parm+' self/solutions/g'+g)

#################################################
## Clear
logging.info('Cleaning...')

check_rm('*bak *last *pickle')

# here images, models, solutions for each group will be saved
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/models'): os.makedirs('self/models')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')
check_rm('self/solutions/gall')
os.makedirs('self/solutions/gall')

check_rm('img')
os.makedirs('img')

mss = sorted(glob.glob('all/all_TC*[0-9].MS'))
concat_ms = 'all/concat.MS'
os.makedirs('logs/all')
nchan = find_nchan(mss[0])

##################################################################################################
## Add model to MODEL_DATA
logging.info('Add model to MODEL_DATA...')
# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for ms in mss:
    check_rm(ms+'/'+sourcedb_basename)
    os.system('cp -r '+sourcedb+' '+ms)
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.sourcedb='+ms+'/'+sourcedb_basename, log=ms+'_pre.log', cmd_type='NDPPP')
s.run(check=True)

# 1. find and remove FR

####################################################################################################
# To circular - SB.MS:DATA -> SB.MS:CORRECTED_DATA (circular)
logging.info('Convert to circular...')
for ms in mss:
    s.add('/home/fdg/scripts/mslin2circ.py -s -w -i '+ms+':DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
s.run(check=True, max_threads=1)

################################################################################################
# Smooth CORRECTED_DATA -> SMOOTHED_DATA (circular, smooth)
logging.info('BL-based smoothing...')
for ms in mss:
    s.add('BLavg.py -r -w -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python', processors='max')
s.run(check=True, max_threads=1)

#################################################################################################
# solve+correct TEC - group*_TC.MS:SMOOTHED_DATA -> group*_TC.MS:CORRECTED_DATA (circular, smooth, TEC-calibrated)
# TODO: merge with next step?
logging.info('Calibrating TEC...')
for ms in mss:
    check_rm(ms+'/instrument-tec')
    s.add('NDPPP '+parset_dir+'/NDPPP-solTEC.parset msin='+ms+' msin.datacolumn=SMOOTHED_DATA cal.parmdb='+ms+'/instrument-tec', log=ms+'_sol-tec.log', cmd_type='NDPPP')
s.run(check=True)

# TODO: BBS for correct, move to NDPPP
for ms in mss:
    s.add('calibrate-stand-alone --parmdb-name instrument-tec '+ms+' '+parset_dir+'/bbs-cor_tec.parset '+skymodel, \
              log=ms+'_cor-tec.log', cmd_type='BBS', processors=2)
s.run(check=True)
# TODO: smooth csp + rerun with only tec?

##############################################################################################
# Solve SB.MS:CORRECTED_DATA (only solve)
logging.info('Calibrating for FR...')
for ms in mss:
    check_rm(ms+'/instrument')
    s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA cal.parmdb='+ms+'/instrument cal.solint=30 cal.nchan=4', log=ms+'_sol-g.log', cmd_type='NDPPP')
s.run(check=True)
        
#################################################################################
# Preapre fake FR parmdb
logging.info('Prepare fake FR parmdb...')
for ms in mss:
    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
s.run(check=True)
for ms in mss:
    s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,24)"', log=ms+'_taql.log', cmd_type='general')
s.run(check=True)

################################################
# Prepare and run losoto
check_rm('globaldb')
check_rm('globaldb-fr')
os.system('mkdir globaldb')
os.system('mkdir globaldb-fr')
for i, ms in enumerate(mss):
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb/')
    if i == 0: os.system('cp -r '+ms+'/ANTENNA '+ms+'/FIELD '+ms+'/sky globaldb-fr/')
    num = re.findall(r'\d+', ms)[-1]
    logging.debug('Copy instrument of '+ms+' into globaldb/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument globaldb/instrument-'+str(num))
    logging.debug('Copy instrument-fr of '+ms+' into globaldb-fr/instrument-'+str(num))
    os.system('cp -r '+ms+'/instrument-fr globaldb-fr/instrument-fr-'+str(num))

logging.info('Running LoSoTo...')
check_rm('plots')
check_rm('cal-fr.h5')
s.add('H5parm_importer.py -v cal-fr.h5 globaldb', log='losoto1.log', cmd_type='python', processors='max')
s.run(check=True)
#s.add('losoto -v cal-fr.h5 '+parset_dir+'/losoto-flag.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
s.add('losoto -v cal-fr.h5 '+parset_dir+'/losoto-fr.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('H5parm_exporter.py -v -t rotationmeasure000 cal-fr.h5 globaldb-fr', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
check_rm('plots-fr')
os.system('mv plots plots-fr')

for i, ms in enumerate(mss):
    num = re.findall(r'\d+', ms)[-1]
    check_rm(ms+'/instrument-fr')
    logging.debug('Copy globaldb-fr/sol000_instrument-fr-'+str(num)+' into '+ms+'/instrument-fr')
    os.system('cp -r globaldb-fr/sol000_instrument-fr-'+str(num)+' '+ms+'/instrument-fr')

sys.exit(1)

#####################################################
# Correct FR SB.MS:DATA->DATA_INIT
logging.info('Faraday rotation correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-corFR.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr', log=ms+'_corFR.log', cmd_type='NDPPP')
s.run(check=True)

# 2: recalibrate without FR

#################################################################################################
## Create columns
logging.info('Creating MODEL_DATA, MODEL_DATA_HIGHRES, SUBTRACTED_DATA...')
for ms in mss:
    s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA,MODEL_DATA_HIGHRES,SUBTRACTED_DATA', log=ms+'_addcol.log', cmd_type='python')
s.run(check=True)

###################################################################################################
# Self-cal cycle
for c in xrange(niter):
    logging.info('Start selfcal cycle: '+str(c))

    #################################################################################################
    # Smooth
    logging.info('BL-based smoothing...')
    for ms in mss:
        s.add('BLavg.py -r -w -i DATA_INIT -o SMOOTHED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python', processors='max')
    s.run(check=True, max_threads=4)

    if c == 0:
        # on first cycle concat (need to be done after smoothing)
        logging.info('Concatenating TCs...')
        check_rm(concat_ms+'*')
        pt.msutil.msconcat(mss, concat_ms, concatTime=False)

    sys.exit(1)
    # TODO: merge with next step?
    # solve+correct TEC - group*_TC.MS:SMOOTHED_DATA -> group*_TC.MS:CORRECTED_DATA
    logging.info('Calibrating TEC...')
    for ms in mss:
        check_rm(ms+'/instrument-tec')
        s.add('NDPPP '+parset_dir+'/NDPPP-solTEC.parset msin='+ms+' msin.datacolumn=SMOOTHED_DATA cal.parmdb='+ms+'/instrument-tec', log=ms+'_sol-tec-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)
    # TODO: BBS for correct, move to NDPPP
    for ms in mss:
        s.add('calibrate-stand-alone --parmdb-name instrument-tec '+ms+' '+parset_dir+'/bbs-cor_tec.parset '+skymodel, \
                  log=ms+'_cor-tec-c'+str(c)+'.log', cmd_type='BBS', processors=2)
    s.run(check=True)

    # solve AMP - group*_TC.MS:CORRECTED_DATA
    logging.info('Calibrating fast amp...')
    for ms in mss:
       check_rm(ms+'/instrument-amp')
       s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA cal.parmdb='+ms+'/instrument-amp cal.solint=2 cal.nchan=0 cal.caltype=commonscalaramplitude',\
            log=ms+'_sol-amp-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # merge parmdbs to plot everything
    logging.info('Merging instrument tables...')
    for ms in mss:
        merge_parmdb(ms+'/instrument-tec', ms+'/instrument-amp', ms+'/instrument', clobber=True)

    # LoSoTo Amp rescaling
    losoto(c, mss, 'all', parset_dir+'/losoto-norm.parset', instrument_in='instrument', instrument_out='instrument')

    # correct amp - CORRECTED_DATA -> CORRECTED_DATA
    logging.info('Correcting fast amp...')
    for ms in mss:
       s.add('NDPPP '+parset_dir+'/NDPPP-corG.parset msin='+ms+' cor.parmdb='+ms+'/instrument cor.correction=commonscalaramplitude',\
            log=ms+'_cor-amp-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    logging.info('Restoring WEIGHT_SPECTRUM before imaging...')
    s.add('taql "update '+concat_ms+' set WEIGHT_SPECTRUM = WEIGHT_SPECTRUM_ORIG"', log='taql-resetweights-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    ###################################################################################################################
    # concat all TCs in one MS - group*_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected)

    # clean mask clean (cut at 8k lambda) - MODEL_DATA updated
    logging.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/wide-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
            -scale 5arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 8000 -mgain 0.6 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 '+concat_ms, \
            log='wscleanA-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask')
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
               params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/LBAsurvey/tooth_mask.crtf', 'setTo':1}, log='casablank-c'+str(c)+'.log')
    s.run(check=True)
    logging.info('Cleaning low resolution (cycle: '+str(c)+')...')
    s.add('wsclean -reorder -name ' + imagename + '-masked -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
            -scale 5arcsec -weight briggs 0.0 -niter 20000 -update-model-required -maxuv-l 8000 -mgain 0.6 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -casamask '+imagename+'.newmask '+concat_ms, \
            log='wscleanB-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    # Convert model to ms and ft() it
    #logging.info('Ft() model...')
    #model = imagename+'-masked-model.fits'
    #s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_fits2ms.py', \
    #            params={'imgs':model, 'del_fits':False}, log='casa_fits2ms-c'+str(c)+'.log')
    #s.run(check=True)
    #s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat_ms, 'model':model.replace('.fits','.ms'), 'wproj':512}, log='ft-c'+str(c)+'.log')
    #s.run(check=True)
   
    ####################################################################
    # FAST VERSION (no low-res)
    #continue
    ####################################################################

    logging.info('Moving MODEL_DATA to MODEL_DATA_HIGHRES...')
    s.add('taql "update '+concat_ms+' set MODEL_DATA_HIGHRES = MODEL_DATA"', log='taql1-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    ############################################################################################################
    # Subtract model from all TCs - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected, high-res model subtracted)
    logging.info('Subtracting high-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
    s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql3-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    # reclean low-resolution
    logging.info('Cleaning low resolution (cycle: '+str(c)+')...')
    imagename = 'img/wide-lr-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 4000 4000 -mem 30 -j '+str(s.max_processors)+' \
            -scale 15arcsec -weight briggs 0.0 -niter 50000 -no-update-model-required -maxuv-l 2500 -mgain 0.6 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 '+concat_ms, \
            log='wscleanA-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask', threshpix=6) # a bit higher treshold
    logging.info('Cleaning low resolution with mask (cycle: '+str(c)+')...')
    s.add('wsclean -reorder -name ' + imagename + '-masked -size 4000 4000 -mem 30 -j '+str(s.max_processors)+' \
            -scale 15arcsec -weight briggs 0.0 -niter 10000 -update-model-required -maxuv-l 2500 -mgain 0.6 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -casamask '+imagename+'.newmask '+concat_ms, \
            log='wscleanB-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    # Convert model to ms and ft() it
    #logging.info('Ft() model...')
    #model = imagename+'-masked-model.fits'
    #s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_fits2ms.py', \
    #            params={'imgs':model, 'del_fits':False}, log='casa_fits2ms-lr-c'+str(c)+'.log')
    #s.run(check=True)
    #s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat_ms, 'model':model.replace('.fits','.ms'), 'wproj':512}, log='ft-lr-c'+str(c)+'.log')
    #s.run(check=True)

    ###############################################################################################################
    # Subtract low-res model - concat.MS:CORRECTED_DATA - MODEL_DATA -> concat.MS:CORRECTED_DATA (empty)
    logging.info('Subtracting low-res model (CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA)...')
    s.add('taql "update '+concat_ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='taql4-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

    # Flag on residuals
    logging.info('Flagging residuals...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms, \
                log=ms+'_flag-c'+str(c)+'.log', cmd_type='NDPPP')
    s.run(check=True)

    # Concat models
    logging.info('Adding model data columns (MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA)...')
    s.add('taql "update '+concat_ms+' set MODEL_DATA = MODEL_DATA_HIGHRES + MODEL_DATA"', log='taql5-c'+str(c)+'.log', cmd_type='general')
    s.run(check=True)

# Perform a final clean to create an inspection image of SUBTRACTED_DATA which should be very empty
logging.info('Empty cleaning...')
s.add('taql "update '+concat_ms+' set SUBTRACTED_DATA = CORRECTED_DATA"', log='taql_sub.log', cmd_type='general')
s.run(check=True)
imagename = 'img/empty'
s.add('wsclean -reorder -name ' + imagename + ' -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
        -scale 5arcsec -weight briggs 0.0 -niter 1 -no-update-model-required -maxuv-l 8000 -mgain 0.6 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -datacolumn SUBTRACTED_DATA '+concat_ms, \
        log='wsclean-empty.log', cmd_type='wsclean', processors='max')
s.run(check=True)

# Copy last *model
logging.info('Copying models/images...')
os.system('mv img/wide-'+str(c)+'-masked-model.fits self/models/wide_g'+g+'.model.fits')
os.system('mv img/wide-lr-'+str(c)+'-masked-model.fits self/models/wide_lr_g'+g+'.model.fits')

s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_fits2ms.py', \
                params={'imgs':['self/models/wide_g'+g+'.model.fits', 'self/models/wide_g'+g+'.model.fits'], 'del_fits':True}, log='casa_fits2ms.log')
s.run(check=True)

# Copy images
[ os.system('mv img/wide-'+str(c)+'.newmask self/images/g'+g) for c in xrange(niter) ]
[ os.system('mv img/wide-lr-'+str(c)+'.newmask self/images/g'+g) for c in xrange(niter) ]
[ os.system('mv img/wide-'+str(c)+'-image.fits self/images/g'+g) for c in xrange(niter) ]
[ os.system('mv img/wide-lr-'+str(c)+'-image.fits self/images/g'+g) for c in xrange(niter) ]
[ os.system('mv img/wide-'+str(c)+'-masked-image.fits self/images/g'+g) for c in xrange(niter) ]
[ os.system('mv img/wide-lr-'+str(c)+'-masked-image.fits self/images/g'+g) for c in xrange(niter) ]
os.system('mv img/empty-image.fits self/images/g'+g)
os.system('mv logs '+group)

# TODO: add final imaging with all SBs

logging.info("Done.")
