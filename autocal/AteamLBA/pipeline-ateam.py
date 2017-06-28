#!/usr/bin/python
# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation
import sys, os, glob, re
import numpy as np
from autocal.lib_pipeline import *
import pyrap.tables as pt

sourcedb = '/home/fdg/scripts/model/A-team_4_CC.skydb'
skymodel = '/home/fdg/scripts/model/A-team_4_CC.skymodel'
patch = 'VirA'
datadir = '/home/fdg/lofar2/LOFAR/Ateam_LBA/VirA/tgts-cycle1-bkp'
bl2flag = 'CS013LBA\;CS031LBA\;RS409LBA\;RS310LBA' # virgo
parset_dir = '/home/fdg/scripts/autocal/AteamLBA/parset_ateam/'
user_mask = '/home/fdg/scripts/autocal/AteamLBA/VirA.reg'

########################################################
logger = set_logger('pipeline-ateam.logger')
check_rm('logs')
check_rm('img')
os.makedirs('img')
s = Scheduler(dry=False)
mss = sorted(glob.glob(datadir+'/*MS'))

############################################################
# Avg to 4 chan and 4 sec
# Remove internationals
# TODO: move to download pipeline
nchan = find_nchan(mss[0])
timeint = find_timeint(mss[0])
if nchan % 4 != 0 and nchan != 1:
    logger.error('Channels should be a multiple of 4.')
    sys.exit(1)

avg_factor_f = nchan / 4 # to 4 ch/SB
if avg_factor_f < 1: avg_factor_f = 1
avg_factor_t = int(np.round(4/timeint)) # to 4 sec
if avg_factor_t < 1: avg_factor_t = 1

if avg_factor_f != 1 or avg_factor_t != 1:
    logger.info('Average in freq (factor of %i) and time (factor of %i)...' % (avg_factor_f, avg_factor_t))
    for ms in mss:
        msout = ms.replace('.MS','-avg.MS').split('/')[-1]
        if os.path.exists(msout): continue
        s.add('NDPPP '+parset_dir+'/NDPPP-avg.parset msin='+ms+' msout='+msout+' msin.datacolumn=DATA avg.timestep='+str(avg_factor_t)+' avg.freqstep='+str(avg_factor_f), \
                log=msout+'_avg.log', cmd_type='NDPPP')
    s.run(check=True)
    nchan = nchan / avg_factor_f
    timeint = timeint * avg_factor_t
else:
    logger.info('Copy data - no averaging...')
    for ms in mss:
        msout = ms.replace('.MS','-avg.MS').split('/')[-1]
        if os.path.exists(msout): continue
        s.add('NDPPP '+parset_dir+'/NDPPP-avg.parset msin='+ms+' msout='+msout+' msin.datacolumn=DATA avg.timestep=1 avg.freqstep=1', \
                log=msout+'_cp.log', cmd_type='NDPPP') # better than cp as activates dysco
    s.run(check=True)

mss = sorted(glob.glob('*.MS'))

#############################################################   
## flag bad stations, and low-elev
logger.info('Flagging...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms+' msout=. flag1.baseline='+bl2flag+' msin.datacolumn=DATA', \
            log=ms+'_flag.log', cmd_type='NDPPP')
s.run(check=True)
 
# Prepare output parmdb
# TODO: remove as soon as losoto has the proper exporter
logger.info('Creating fake parmdb...')
for ms in mss:
    if os.path.exists(ms+'/instrument-clock'): continue
    s.add('calibrate-stand-alone -f --parmdb-name instrument-clock '+ms+' '+parset_dir+'/bbs-fakeparmdb-clock.parset '+skymodel, log=ms+'_fakeparmdb-clock.log', cmd_type='BBS')
s.run(check=True)
for ms in mss:
    if os.path.exists(ms+'/instrument-fr'): continue
    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
s.run(check=True)
for ms in mss:
    s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,24)"', log=ms+'_taql.log', cmd_type='general')
s.run(check=True)

# predict to save time ms:MODEL_DATA
logger.info('Predict...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.sourcedb='+sourcedb+' pre.sources='+patch, log=ms+'_pre.log', cmd_type='NDPPP')
s.run(check=True)

for c in xrange(10):

    #################################################
    # 1: find the FR and remve it
    
    # Beam correction DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
    s.run(check=True)
    
    # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Converting to circular...')
    for ms in mss:
        s.add('mslin2circ.py -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
    s.run(check=True)
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    for ms in mss:
        s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth1.log', cmd_type='python')
    s.run(check=True)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    for ms in mss:
        check_rm(ms+'/instrument')
        s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms, log=ms+'_sol1.log', cmd_type='NDPPP')
    s.run(check=True)
    
    run_losoto(s, 'fr-c'+str(c), mss, [parset_dir+'/losoto-fr.parset'], outtab='rotationmeasure000', \
        inglobaldb='globaldb', outglobaldb='globaldb-fr', ininstrument='instrument', outinstrument='instrument-fr', putback=True)
    
    #####################################################
    # 2: find CROSS DELAY and remve it
    
    # Beam correction DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
    s.run(check=True)
    
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', log=ms+'_corFR.log', cmd_type='NDPPP')
    s.run(check=True)
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    for ms in mss:
        s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth2.log', cmd_type='python')
    s.run(check=True)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    for ms in mss:
        check_rm(ms+'/instrument')
        s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms, log=ms+'_sol2.log', cmd_type='NDPPP')
    s.run(check=True)
    
    run_losoto(s, 'cd-c'+str(c), mss, [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-cd.parset'], outtab='amplitude000,crossdelay', \
        inglobaldb='globaldb', outglobaldb='globaldb', ininstrument='instrument', outinstrument='instrument-cd', putback=True)
    
    #################################################
    # 3: recalibrate without FR
    
    # Beam correction (and update weight in case of imaging) DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms+' corrbeam.updateweight=True', log=ms+'_beam2.log', cmd_type='NDPPP')
    s.run(check=True)
    
    # Correct DELAY CORRECTED_DATA (beam corrected) -> CORRECTED_DATA
    logger.info('Cross delay correction...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-cd cor.correction=gain', log=ms+'_corCD.log', cmd_type='NDPPP')
    s.run(check=True)
    
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', log=ms+'_corFR.log', cmd_type='NDPPP')
    s.run(check=True)
    
    # Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
    logger.info('BL-smooth...')
    for ms in mss:
        s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth3.log', cmd_type='python')
    s.run(check=True)
    
    # Solve cal_SB.MS:SMOOTHED_DATA (only solve)
    logger.info('Calibrating...')
    for ms in mss:
        check_rm(ms+'/instrument')
        s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms, log=ms+'_sol3.log', cmd_type='NDPPP')
    s.run(check=True)
    
    run_losoto(s, 'final-c'+str(c), mss, [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-amp.parset',parset_dir+'/losoto-ph.parset'], outtab='amplitudeOrig000,phaseOrig000', \
               inglobaldb='globaldb', outglobaldb='globaldb', ininstrument='instrument', outinstrument='instrument', putback=True)
    
    from make_mask import make_mask
    
    # Correct all CORRECTED_DATA (beam, CD, FR corrected) -> CORRECTED_DATA
    logger.info('Amp/ph correction...')
    for ms in mss:
        if c == 0:
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.updateweights=True cor.parmdb='+ms+'/instrument cor.correction=gain', log=ms+'_corG.log', cmd_type='NDPPP')
        else:
            # update weight only first time, it should be at first order correct
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.updateweights=False cor.parmdb='+ms+'/instrument cor.correction=gain', log=ms+'_corG.log', cmd_type='NDPPP')
    s.run(check=True)
    
#    logger.info('Cleaning...')
#    imagename = 'img/wide-c'+str(c)
#    s.add('wsclean -reorder -name ' + imagename + ' -size 3500 3500 -trim 3000 3000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
#            -scale 6arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -mgain 0.9 \
#            -pol I -joinchannels -fit-spectral-pol 3 -channelsout 10 -auto-threshold 20 '+' '.join(mss), \
#            log='wscleanA-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
#    s.run(check=True)
#    
#    # make mask
#    maskname = imagename+'-mask.fits'
#    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 3, atrous_do=True)
#    if user_mask is not None:
#        blank_image_reg(maskname, user_mask, inverse=False, blankval=1)
    
    logger.info('Cleaning w/ mask')
    imagename = 'img/wideM-c'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 1200 1200 -trim 800 800 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 1.5 \
            -scale 4arcsec -weight briggs -0.5 -niter 100000 -no-update-model-required -mgain 0.7 \
            -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,4,8,16,32 -auto-mask 5\
            -pol I -joinchannels -fit-spectral-pol 3 -channelsout 15 -threshold 0.005 '+' '.join(mss), \
            log='wscleanB-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)
    sys.exit(1)

    logger.info('Predict (ft)...')
    s.add('wsclean -predict -name ' + imagename + ' -mem 90 -j '+str(s.max_processors)+' -channelsout 15 '+' '.join(mss), \
            log='wscleanPRE-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

logger.info("Done.")
