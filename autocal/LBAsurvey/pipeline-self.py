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
globaldb = 'globaldb'
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

def losoto(c, mss, g, parset, instrument='instrument'):
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
        for instTab in glob.glob(ms+'/'+instrument):
            os.system('cp -r '+instTab+' globaldb/'+instTab.split('/')[-1]+'-'+str(num))

    h5parm = 'global-c'+str(c)+'.h5'

    s.add('H5parm_importer.py -v '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('losoto -v '+h5parm+' '+parset, log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=False)
    s.add('H5parm_exporter.py -v -c '+h5parm+' globaldb', log='losoto-c'+str(c)+'.log', log_append=True, cmd_type='python')
    s.run(check=True)

    for num, ms in enumerate(mss):
        check_rm(ms+'/instrument')
        for instTab in glob.glob('globaldb/sol000_'+instrument+'-'+str(num)):
            os.system('mv '+instTab+' '+ms+'/'+instTab.replace('globaldb/sol000_','').rsplit('-',1)[0])
    os.system('mv plots self/solutions/g'+g+'/plots-c'+str(c))
    os.system('mv '+h5parm+' self/solutions/g'+g)

#################################################
## Clear
logging.info('Cleaning...')

check_rm('*bak *last *pickle')
mss = sorted(glob.glob('*MS'))
nchan = find_nchan(mss[0])

# here images, models, solutions for each group will be saved
if not os.path.exists('self/images'): os.makedirs('self/images')
if not os.path.exists('self/models'): os.makedirs('self/models')
if not os.path.exists('self/solutions'): os.makedirs('self/solutions')
check_rm('self/solutions/gall')
os.makedirs('self/solutions/gall')

# copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
sourcedb_basename = sourcedb.split('/')[-1]
for ms in mss:
    check_rm(ms+'/'+sourcedb_basename)
    os.system('cp -r '+sourcedb+' '+ms)

# 1. initial processing

###############################################
# Initial processing
logging.info('Fixing beam table...')
for ms in mss:
    s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
s.run(check=False)

###################################################
# Beam correction DATA -> CORRECTED_DATA (beam corrected)
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
s.run(check=True)

##########################################################################################
# Copy instrument tables
for ms in mss:
    num = re.findall(r'\d+', ms)[-1]
    check_rm(ms+'/instrument')
    logging.debug('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')
    os.system('cp -r '+globaldb+'/sol000_instrument-'+str(num)+' '+ms+'/instrument')

# Apply cal sol (clock, phase offsetts & bandpass) - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (calibrator corrected data, beam corrected, lin)
logging.info('Apply solutions...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-initcor.parset msin='+ms+' cor1.parmdb='+ms+'/instrument'+' cor2.parmdb='+ms+'/instrument', log=ms+'_cor.log', cmd_type='NDPPP')
s.run(check=True)

##############################################################
# TODO: create a virtual dataset for flags
# Flagging on concatenated dataset
logging.info('Flagging...')
for groupname in groupnames:
    s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+groupname+'/'+groupname+'.MS', \
                log=groupname+'_NDPPP_flag.log', cmd_type='NDPPP')
s.run(check=True)

# 2. find and remove FR

####################################################################################################
# To circular - SB.MS:CORRECTED_DATA -> SB.MS:CIRC_DATA (beam corrected, circular)
logging.info('Convert to circular...')
for ms in mss:
    s.add('/home/fdg/scripts/mslin2circ.py -s -i '+ms+':CORRECTED_DATA -o '+ms+':CIRC_DATA', log=ms+'_circ2lin.log', cmd_type='python')
s.run(check=True)

#################################################################################################
# Smooth CIRC_DATA -> SMOOTHED_DATA
logging.info('BL-based smoothing...')
for ms in mss:
    s.add('BLavg.py -r -w -i CIRC_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python', processors='max')
s.run(check=True)

#################################################################################################
# solve+correct TEC - group*_TC.MS:SMOOTHED_DATA -> group*_TC.MS:CORRECTED_DATA
logging.info('Calibrating TEC...')
for ms in mss:
    s.add('calibrate-stand-alone -f --parmdb-name instrument-tec '+ms+' '+parset_dir+'/bbs-solcor_tec.parset '+skymodel, \
              log=ms+'_solcor_tec.log', cmd_type='BBS', processors=2)
s.run(check=True)

###############################################################################################
#  Solve SB.MS:CORRECTED_DATA (only solve)
logging.debug('Calibrating for FR...')
for ms in mss:
    check_rm(ms+'/instrument')
    s.add('NDPPP '+parset_dir+'/NDPPP-solG.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA cal.parmdb='+ms+'/instrument cal.sourcedb='+ms+'/'+sourcedb 'cal.solint=10 cal.nchan=4', log=ms+'_sol-circ.log', cmd_type='NDPPP')
s.run(check=True)

#################################################################################
# Preapre fake FR parmdb
logging.debug('Prepare fake FR parmdb...')
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
check_rm('cal1.h5')
s.add('H5parm_importer.py -v cal1.h5 globaldb', log='losoto1.log', cmd_type='python', processors='max')
s.run(check=True)
#s.add('losoto -v cal1.h5 '+parset_dir+'/losoto-flag.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
s.add('losoto -v cal1.h5 '+parset_dir+'/losoto-fr.parset', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
s.add('H5parm_exporter.py -v -t rotationmeasure000 cal1.h5 globaldb-fr', log='losoto1.log', log_append=True, cmd_type='python', processors='max')
s.run(check=True)
check_rm('plots-fr')
os.system('mv plots plots-fr')

for i, ms in enumerate(mss):
    num = re.findall(r'\d+', ms)[-1]
    check_rm(ms+'/instrument-fr')
    logging.debug('Copy globaldb-fr/sol000_instrument-fr-'+str(num)+' into '+ms+'/instrument-fr')
    os.system('cp -r globaldb-fr/sol000_instrument-fr-'+str(num)+' '+ms+'/instrument-fr')

######################################################
# Correct FR
logging.info('Faraday rotation correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-corFR.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr', log=ms+'_corFR.log', cmd_type='NDPPP')
s.run(check=True)

# 3: recalibrate without FR

###############################################
# Beam correction CORRECTED_DATA -> CORRECTED_DATA
logging.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA', log=ms+'_beam2.log', cmd_type='NDPPP')
s.run(check=True)

#################################################################################################
# Smooth CORRECTED_DATA -> SMOOTHED_DATA
logging.info('BL-based smoothing...')
for ms in mss:
    s.add('BLavg.py -r -w -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth2.log', cmd_type='python', processors='max')
s.run(check=True)




# TODO: add below to above and rewrite the init to use time-splitted data


##################################################################################################
# calibrate slow (3 min/16 chan) diag ph for FR
chanblock = 4#16 # number of channel to solve at the same time
assert nchan % chanblock == 0
logging.debug('Calibrating FR - iterating on '+str(1.*nchan/chanblock)+' channel blocks.')
# find how many channels per block to put together. If nchan%chanblock != 0, add one to initial runs
chansets = [chanblock for ch in xrange(nchan/chanblock)]
for i in range(nchan%chanblock): chansets[i] = chanblock + 1
logging.debug('Chansets: '+str(chansets)+' - total runs: '+str(len(chansets)))
startchan = 0
for i, chan in enumerate(chansets):
    logging.debug('Channel/16: '+str(i))
    for ms in mss:
        check_rm(ms+'/instrument-ph'+str(i))
        s.add('NDPPP '+parset_dir+'/NDPPP-sol_ph.parset msin='+ms+' msin.startchan='+str(startchan)+' msin.nchan='+str(chan)+' cal.solint=36 cal.parmdb='+ms+'/instrument-ph'+str(i)+' cal.sourcedb='+ms+'/'+sourcedb_basename, log=ms+'chan'+str(i)+'_solph.log', cmd_type='NDPPP')
    s.run(check=True)
    startchan += chan

losoto(0, mss, 'all', parset_dir+'/losoto-fr.parset', 'instrument-ph*')
sys.exit(1)

for group in sorted(glob.glob('group*'))[::-1]:

    mss = sorted(glob.glob(group+'/group*_TC*[0-9].MS'))
    concat_ms = group+'/concat.MS'
    g = str(re.findall(r'\d+', mss[0])[0])
    os.makedirs('logs/'+group)
    logging.info('Working on group: '+g+'...')

    # copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
    for ms in mss:
        check_rm(ms+'/'+sourcedb_basename)
        os.system('cp -r '+sourcedb+' '+ms)
    
    ################################################################################################
    # Clear
    logging.info('Cleaning...')
    check_rm(group+'/plots* plots')
    check_rm(group+'/*h5 *h5 globaldb')
    check_rm(group+'/logs')
    check_rm('img')
    os.makedirs('img')
    check_rm('self/images/g'+g)
    os.makedirs('self/images/g'+g)
    check_rm('self/solutions/g'+g)
    os.makedirs('self/solutions/g'+g)

    #################################################################################################
    # Create columns
    logging.info('Creating MODEL_DATA_HIGHRES, SUBTRACTED_DATA, MODEL_DATA and CORRECTED_DATA...')
    for ms in mss:
        s.add('addcol2ms.py -m '+ms+' -c MODEL_DATA,CORRECTED_DATA,MODEL_DATA_HIGHRES,SUBTRACTED_DATA', log=ms+'_addcol.log', cmd_type='python')
    s.run(check=True)

    ###################################################################################################
    # Create rotationmeasure fakeparmdb
    for ms in mss:
        s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
    s.run(check=True)

    # fill with solutions TODO

    # Correct RM and beam DATA -> DATA_INIT
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-cor_beam-fr.parset msin='+ms+' cal.parmdb='+ms+'/instrument-fr',\
            log=ms+'_cor-beam-fr.log', cmd_type='NDPPP')
    s.run(check=True)

    ####################################################################################################
    # Self-cal cycle
    for c in xrange(niter):
        logging.info('Start selfcal cycle: '+str(c))

        #################################################################################################
        # Smooth
        logging.info('BL-based smoothing...')
        for ms in mss:
            s.add('BLavg.py -r -w -i DATA_INIT -o SMOOTHED_DATA '+ms, log=ms+'_smooth-c'+str(c)+'.log', cmd_type='python')
        s.run(check=True)

        if c == 0:
            # on first cycle concat (need to be done after smoothing)
            logging.info('Concatenating TCs...')
            check_rm(concat_ms+'*')
            pt.msutil.msconcat(mss, concat_ms, concatTime=False)

        # solve+correct TEC - group*_TC.MS:SMOOTHED_DATA -> group*_TC.MS:CORRECTED_DATA
        logging.info('Calibrating TEC...')
        for ms in mss:
            s.add('calibrate-stand-alone -f --parmdb-name instrument-tec '+ms+' '+parset_dir+'/bbs-solcor_tec.parset '+skymodel, \
                      log=ms+'_solcor_tec-c'+str(c)+'.log', cmd_type='BBS', processors=2)
        s.run(check=True)

        sys.exit(1)
 
        # solve+correct AMP - group*_TC.MS:CORRECTED_DATA -> group*_TC.MS:CORRECTED_DATA
        # TODO: update to scalar amplitude add channles if no groups?
        logging.info('Calibrating fast amp...')
        for ms in mss:
           check_rm(ms+'/instrument-amp')
           s.add('NDPPP '+parset_dir+'/NDPPP-sol_amp.parset msin='+ms+' cal.parmdb='+ms+'/instrument-amp cal.sourcedb='+ms+'/'+sourcedb_basename,\
                log=ms+'_solamp-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        # LoSoTo Amp rescaling / Ph=0
        losoto(c, mss, g, parset_dir+'/losoto-amp.parset')

        # merge parmdbs
        logging.info('Merging instrument tables...')
        for ms in mss:
            merge_parmdb(ms+'/instrument-tec', ms+'/instrument-amp', ms+'/instrument', clobber=True)

        # correct TEC phases + amp - DATA_INIT -> CORRECTED_DATA
        # TODO: TAQL to copy CORRECTED_DATA_TMP[chans] back
        for ms in mss:
           s.add('NDPPP '+parset_dir+'/NDPPP-cor_amp.parset msin='+ms+' cor.parmdb='+ms+'/instrument-amp',\
                log=ms+'_coramp-c'+str(c)+'.log', cmd_type='NDPPP')
        s.run(check=True)

        sys.exit(1)
        
        logging.info('Restoring WEIGHT_SPECTRUM before imaging...')
        s.add('taql "update '+concat_ms+' set WEIGHT_SPECTRUM = WEIGHT_SPECTRUM_ORIG"', log='taql-resetweights-c'+str(c)+'.log', cmd_type='general')
        s.run(check=True)
    
        ###################################################################################################################
        # concat all TCs in one MS - group*_TC.MS:CORRECTED_DATA -> concat.MS:CORRECTED_DATA (selfcal corrected, beam corrected)
    
        # clean mask clean (cut at 8k lambda) - MODEL_DATA updated
        logging.info('Cleaning (cycle: '+str(c)+')...')
        imagename = 'img/wide-'+str(c)
        s.add('wsclean -reorder -name ' + imagename + ' -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
                -scale 5arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -maxuv-l 8000 -mgain 0.85 '+concat_ms, \
                log='wscleanA-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
        make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask')
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_blank.py', \
                   params={'imgs':imagename+'.newmask', 'region':'/home/fdg/scripts/autocal/LBAsurvey/tooth_mask.crtf', 'setTo':1}, log='casablank-c'+str(c)+'.log')
        s.run(check=True)
        logging.info('Cleaning low resolution (cycle: '+str(c)+')...')
        s.add('wsclean -reorder -name ' + imagename + '-masked -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
                -scale 5arcsec -weight briggs 0.0 -niter 20000 -no-update-model-required -maxuv-l 8000 -mgain 0.85 -casamask '+imagename+'.newmask '+concat_ms, \
                log='wscleanB-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)

        # Convert model to ms and ft() it
        logging.info('Ft() model...')
        model = imagename+'-masked-model.fits'
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_fits2ms.py', \
                    params={'imgs':model, 'del_fits':False}, log='casa_fits2ms-c'+str(c)+'.log')
        s.run(check=True)
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat_ms, 'model':model.replace('.fits','.ms'), 'wproj':512}, log='ft-c'+str(c)+'.log')
        s.run(check=True)
       
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
        s.add('wsclean_1.8 -reorder -name ' + imagename + ' -size 4000 4000 -mem 30 -j '+str(s.max_processors)+'\
                -scale 15arcsec -weight briggs 0.0 -niter 50000 -no-update-model-required -maxuv-l 2500 -mgain 0.85 '+concat_ms, \
                log='wscleanA-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)
        make_mask(image_name = imagename+'-image.fits', mask_name = imagename+'.newmask', threshpix=6) # a bit higher treshold
        logging.info('Cleaning low resolution with mask (cycle: '+str(c)+')...')
        s.add('wsclean_1.8 -reorder -name ' + imagename + '-masked -size 4000 4000 -mem 30 -j '+str(s.max_processors)+' \
                -scale 15arcsec -weight briggs 0.0 -niter 10000 -no-update-model-required -maxuv-l 2500 -mgain 0.85 -casamask '+imagename+'.newmask '+concat_ms, \
                log='wscleanB-lr-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
        s.run(check=True)

        # Convert model to ms and ft() it
        logging.info('Ft() model...')
        model = imagename+'-masked-model.fits'
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_fits2ms.py', \
                    params={'imgs':model, 'del_fits':False}, log='casa_fits2ms-lr-c'+str(c)+'.log')
        s.run(check=True)
        s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':concat_ms, 'model':model.replace('.fits','.ms'), 'wproj':512}, log='ft-lr-c'+str(c)+'.log')
        s.run(check=True)
 
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
    s.add('wsclean_1.8 -reorder -name ' + imagename + ' -size 5000 5000 -mem 30 -j '+str(s.max_processors)+' \
            -scale 5arcsec -weight briggs 0.0 -niter 1 -no-update-model-required -maxuv-l 8000 -mgain 0.85 -datacolumn SUBTRACTED_DATA '+concat_ms, \
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
