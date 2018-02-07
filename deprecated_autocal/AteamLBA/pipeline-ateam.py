#!/usr/bin/python
# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation
import sys, os, glob, re
import numpy as np
from autocal.lib_pipeline import *
import pyrap.tables as pt

if 'VirA2013' in os.getcwd():
    patch = 'VirA'
    datadir = '/home/fdg/lofar2/LOFAR/Ateam_LBA/VirA/tgts2013-bkp'
    bl2flag = 'CS013LBA\;CS031LBA'
    blrange = '[0,1e30]'
elif 'VirA2015' in os.getcwd():
    patch = 'VirA'
    datadir = '/home/fdg/lofar2/LOFAR/Ateam_LBA/VirA/tgts2015-bkp'
    bl2flag = 'CS017LBA\;RS407LBA'
    blrange = '[0,1e30]'
elif 'VirA2017' in os.getcwd():
    patch = 'VirA'
    datadir = '/home/fdg/lofar2/LOFAR/Ateam_LBA/VirA/tgts2017-bkp'
    bl2flag = ''
    blrange = '[0,1e30]'
elif 'TauA' in os.getcwd():
    patch = 'TauA'
    datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/TauA/tgts-bkp'
    bl2flag = 'RS310LBA\;RS210LBA\;RS407LBA\;RS409LBA'
    blrange = '[0,1000,5000,1e30]'
elif 'CasA' in os.getcwd():
    patch = 'CasA'
    datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/CasA/tgts1-bkp'
    #datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/CasA/tgts2-bkp'
    bl2flag = 'CS031LBA'
    blrange = '[0,30000]'
elif 'CygA' in os.getcwd():
    patch = 'CygA'
    datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/CygA/tgts1-bkp'
    #datadir='/home/fdg/lofar2/LOFAR/Ateam_LBA/CygA/tgts2-bkp'
    bl2flag = 'CS031LBA'
    blrange = '[0,30000]'

parset_dir = '/home/fdg/scripts/autocal/AteamLBA/parset_ateam/'

########################################################
logger = set_logger('pipeline-ateam.logger')
check_rm('logs')
check_rm('img')
os.makedirs('img')
s = Scheduler(dry=False)
mss = sorted(glob.glob(datadir+'/*MS'))
mss = mss[len(mss)*2/5:len(mss)*4/5] # use only upper half of the band

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
    if os.path.exists(ms+'/instrument-fr'): continue
    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset /home/fdg/scripts/model/calib-simple.skymodel', log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
s.run(check=True)
for ms in mss:
    s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,24)"', log=ms+'_taql.log', cmd_type='general')
s.run(check=True)

# predict to save time ms:MODEL_DATA
if os.path.exists('/home/fdg/scripts/model/AteamLBA/'+patch+'/wideM-MFS-model.fits'):
    logger.info('Predict (wsclean)...')
    s.add('wsclean -predict -name /home/fdg/scripts/model/AteamLBA/'+patch+'/wideM -mem 90 -j '+str(s.max_processors)+' -channelsout 15 '+' '.join(mss), \
          log='wscleanPRE-init.log', cmd_type='wsclean', processors='max')
else:
    logger.info('Predict (NDPPP)...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.sourcedb=/home/fdg/scripts/model/A-team_4_CC.skydb pre.sources='+patch, log=ms+'_pre.log', cmd_type='NDPPP')
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
        s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' filter.blrange='+blrange, log=ms+'_sol1.log', cmd_type='NDPPP')
    s.run(check=True)
    
    run_losoto(s, 'fr-c'+str(c), mss, [parset_dir+'/losoto-fr.parset'], outtab='rotationmeasure000', \
        inglobaldb='globaldb', outglobaldb='globaldb-fr', ininstrument='instrument', outinstrument='instrument-fr', putback=True)
    
    #####################################################
    # 2: find CROSS DELAY
    
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
        s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' filter.blrange='+blrange, log=ms+'_sol2.log', cmd_type='NDPPP')
    s.run(check=True)
    
    run_losoto(s, 'cd-c'+str(c), mss, [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-amp.parset',parset_dir+'/losoto-cd.parset'], outtab='amplitudeSmooth000,crossdelay', \
        inglobaldb='globaldb', outglobaldb='globaldb', ininstrument='instrument', outinstrument='instrument-cd', putback=True)
    
    #################################################
    # 3: recalibrate without FR

    # Correct DELAY + ampBP DATA (beam corrected) -> CORRECTED_DATA
    logger.info('Cross delay+ampBP correction...')
    for ms in mss:
        if c == 0:
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' msin.datacolumn=DATA cor.updateweights=True cor.parmdb='+ms+'/instrument-cd cor.correction=gain', log=ms+'_corCD.log', cmd_type='NDPPP')
        else:
            s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' msin.datacolumn=DATA cor.updateweights=False cor.parmdb='+ms+'/instrument-cd cor.correction=gain', log=ms+'_corCD.log', cmd_type='NDPPP')
    s.run(check=True)
 
    # Beam correction (and update weight in case of imaging) CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Beam correction...')
    for ms in mss:
        if c == 0:
            s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA corrbeam.updateweights=True', log=ms+'_beam2.log', cmd_type='NDPPP')
        else:
            s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA corrbeam.updateweights=False', log=ms+'_beam2.log', cmd_type='NDPPP')
    s.run(check=True)
       
    # Correct FR CORRECTED_DATA -> CORRECTED_DATA
    logger.info('Faraday rotation correction...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', log=ms+'_corFR.log', cmd_type='NDPPP')
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
        s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms+' filter.blrange='+blrange, log=ms+'_sol3.log', cmd_type='NDPPP')
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
    
    # briggs: -1.2 for virgo
    logger.info('Cleaning (cycle %i)...' % c)
    imagename = 'img/wideM-c'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size 1700 1700 -trim 1500 1500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 1.5 \
            -scale 2arcsec -weight briggs -1.5 -niter 100000 -no-update-model-required -mgain 0.7 \
            -multiscale -multiscale-scale-bias 0.5 -multiscale-scales 0,4,8,16,32 -auto-mask 5\
            -pol I -joinchannels -fit-spectral-pol 3 -channelsout 15 -threshold 0.005 '+' '.join(mss), \
            log='wscleanB-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    logger.info('Predict (ft)...')
    s.add('wsclean -predict -name ' + imagename + ' -mem 90 -j '+str(s.max_processors)+' -channelsout 15 '+' '.join(mss), \
            log='wscleanPRE-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

#    logger.info('Sub model...')
#    for ms in mss:
#        s.add('taql "update '+ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log=ms+'_taql.log', cmd_type='general')
#    s.run(check=True)
#
#    logger.info('Cleaning sub (cycle %i)...' % c)
#    imagename = 'img/wideMsub-c'+str(c)
#    s.add('wsclean -reorder -name ' + imagename + ' -size 500 500 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 1.5 \
#            -scale 10arcsec -weight briggs 0.0 -niter 10000 -no-update-model-required -mgain 0.7 -taper-gaussian 45arcsec \
#            -pol I -joinchannels -fit-spectral-pol 3 -channelsout 15 '+' '.join(mss), \
#            log='wscleanB-c'+str(c)+'.log', cmd_type='wsclean', processors='max')
#    s.run(check=True)

logger.info("Done.")
