#!/usr/bin/python
# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation
import sys, os, glob, re
import numpy as np
from astropy.time import Time
from autocal.lib_pipeline import *
import pyrap.tables as pt

parset_dir = '/home/fdg/scripts/autocal/parset_cal/'
skymodel = '/home/fdg/scripts/model/calib-simple.skymodel'
imaging = True
clock = False

if 'tooth' in os.getcwd(): # tooth 2013
    datadir = '../cals-bkp/'
    bl2flag = 'CS031LBA\;RS409LBA\;RS310LBA'
elif 'bootes' in os.getcwd(): # bootes 2013
    datadir = '../cals-bkp/'
    bl2flag = 'CS013LBA\;CS031LBA\;RS409LBA'
elif 'survey' in os.getcwd():
    obs = os.getcwd().split('/')[-2] # assumes .../c??-o??/3c196
    calname = os.getcwd().split('/')[-1] # assumes .../c??-o??/3c196
    datadir = '../../download/%s/%s' % (obs, calname)
    bl2flag = 'CS031LBA\;RS310LBA\;RS210LBA\;RS409LBA'
    if 'c07-o00' in os.getcwd() or 'c07-o01' in os.getcwd() or 'c07-o02' in os.getcwd() or 'c07-o03' in os.getcwd() or 'c07-o04' in os.getcwd() or 'c07-o05' in os.getcwd() or 'c07-o06' in os.getcwd():
        bl2flag = 'CS031LBA\;RS310LBA\;RS210LBA\;RS409LBA\;RS407LBA'
else:
    datadir = '../cals-bkp/'
    bl2flag = 'RS310HBA'

########################################################
logger = set_logger('pipeline-cal.logger')
check_rm('logs')
s = Scheduler(dry=False)
mss = sorted(glob.glob(datadir+'/*MS'))
calname = mss[0].split('/')[-1].split('_')[0].lower()
logger.info("Calibrator name: %s." % calname)

sourcedb = '/home/fdg/scripts/model/calib-simple.skydb'

if calname == '3c196':
    patch = '3C196'
elif calname == '3c380':
    patch = '3C380'
elif calname == '3c295':
    patch = '3C295'
elif calname == 'CygA':
    sourcedb = '/home/fdg/scripts/model/A-team_4_CC.skydb'
    patch = 'CygA'
else:
    logger.error("Calibrator not recognised.")
    sys.exit(1)

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
## flag bad stations, flags will propagate
#logger.info('Flagging...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-flag.parset msin='+ms+' msout=. flag1.baseline='+bl2flag+' msin.datacolumn=DATA', \
#            log=ms+'_flag.log', cmd_type='NDPPP')
#s.run(check=True)
# 
## Prepare output parmdb
## TODO: remove as soon as losoto has the proper exporter
#logger.info('Creating fake parmdb...')
#for ms in mss:
#    if os.path.exists(ms+'/instrument-clock'): continue
#    s.add('calibrate-stand-alone -f --parmdb-name instrument-clock '+ms+' '+parset_dir+'/bbs-fakeparmdb-clock.parset '+skymodel, log=ms+'_fakeparmdb-clock.log', cmd_type='BBS')
#s.run(check=True)
#for ms in mss:
#    if os.path.exists(ms+'/instrument-fr'): continue
#    s.add('calibrate-stand-alone -f --parmdb-name instrument-fr '+ms+' '+parset_dir+'/bbs-fakeparmdb-fr.parset '+skymodel, log=ms+'_fakeparmdb-fr.log', cmd_type='BBS')
#s.run(check=True)
#
## are we doing HBA?
#tab = pt.table(mss[0]+'/instrument-fr/NAMES/', ack=False)
#HBA = 'HBA' in tab.getcol('NAME')[0]
#tab.close()
#
#for ms in mss:
#    if HBA:
#        s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,25)"', log=ms+'_taql.log', cmd_type='general')
#    else:
#        s.add('taql "update '+ms+'/instrument-fr::NAMES set NAME=substr(NAME,0,24)"', log=ms+'_taql.log', cmd_type='general')
#s.run(check=True)
#
## predict to save time ms:MODEL_DATA
#logger.info('Predict...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-predict.parset msin='+ms+' pre.sourcedb='+sourcedb+' pre.sources='+patch, log=ms+'_pre.log', cmd_type='NDPPP')
#s.run(check=True)
#
##################################################
## 1: find the FR and remve it
#
## Beam correction DATA -> CORRECTED_DATA
#logger.info('Beam correction...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
#s.run(check=True)
#
## Convert to circular CORRECTED_DATA -> CORRECTED_DATA
#logger.info('Converting to circular...')
#for ms in mss:
#    s.add('mslin2circ.py -i '+ms+':CORRECTED_DATA -o '+ms+':CORRECTED_DATA', log=ms+'_circ2lin.log', cmd_type='python')
#s.run(check=True)
#
## Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
#logger.info('BL-smooth...')
#for ms in mss:
#    s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth1.log', cmd_type='python')
#s.run(check=True)
#
## Solve cal_SB.MS:SMOOTHED_DATA (only solve)
#logger.info('Calibrating...')
#for ms in mss:
#    check_rm(ms+'/instrument')
#    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms, log=ms+'_sol1.log', cmd_type='NDPPP')
#s.run(check=True)
#
#run_losoto(s, 'fr', mss, [parset_dir+'/losoto-fr.parset'], outtab='rotationmeasure000', \
#    inglobaldb='globaldb', outglobaldb='globaldb-fr', ininstrument='instrument', outinstrument='instrument-fr', putback=True)
#
######################################################
## 2: find CROSS DELAY and amplitude
#
## Beam correction DATA -> CORRECTED_DATA
#logger.info('Beam correction...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms, log=ms+'_beam.log', cmd_type='NDPPP')
#s.run(check=True)
#
## Correct FR CORRECTED_DATA -> CORRECTED_DATA
#logger.info('Faraday rotation correction...')
#for ms in mss:
#    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.parmdb='+ms+'/instrument-fr cor.correction=RotationMeasure', log=ms+'_corFR.log', cmd_type='NDPPP')
#s.run(check=True)
#
## Smooth data CORRECTED_DATA -> SMOOTHED_DATA (BL-based smoothing)
#logger.info('BL-smooth...')
#for ms in mss:
#    s.add('BLsmooth.py -r -i CORRECTED_DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth2.log', cmd_type='python')
#s.run(check=True)
#
## Solve cal_SB.MS:SMOOTHED_DATA (only solve)
#logger.info('Calibrating...')
#for ms in mss:
#    check_rm(ms+'/instrument')
#    s.add('NDPPP '+parset_dir+'/NDPPP-sol.parset msin='+ms, log=ms+'_sol2.log', cmd_type='NDPPP')
#s.run(check=True)
#
#run_losoto(s, 'cd', mss, [parset_dir+'/losoto-flag.parset',parset_dir+'/losoto-amp.parset',parset_dir+'/losoto-cd.parset'], outtab='amplitudeSmooth000,crossdelay', \
#    inglobaldb='globaldb', outglobaldb='globaldb', ininstrument='instrument', outinstrument='instrument-cd', putback=True)

#################################################
# 3: recalibrate without FR

# Correct cd+amp DATA -> CORRECTED_DATA
logger.info('Cross delay+ampBP correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' msin.datacolumn=DATA cor.parmdb='+ms+'/instrument-cd cor.correction=gain', log=ms+'_corCD.log', cmd_type='NDPPP')
s.run(check=True)

# Beam correction (and update weight in case of imaging) CORRECTED_DATA -> CORRECTED_DATA
# TEST: wegith -> true
logger.info('Beam correction...')
for ms in mss:
    s.add('NDPPP '+parset_dir+'/NDPPP-beam.parset msin='+ms+' msin.datacolumn=CORRECTED_DATA corrbeam.updateweights=False', log=ms+'_beam2.log', cmd_type='NDPPP')
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

if clock:
    run_losoto(s, 'iono', mss, [parset_dir+'/losoto-iono.parset'], outtab='amplitudeSmooth000,phase000,clock000', \
           inglobaldb='globaldb', outglobaldb='globaldb-clock', ininstrument='instrument', outinstrument='instrument-clock', putback=False)
else:
    run_losoto(s, 'iono', mss, [parset_dir+'/losoto-iono.parset'], outtab='amplitude000,phaseOrig000', \
           inglobaldb='globaldb', outglobaldb='globaldb', ininstrument='instrument', outinstrument='instrument', putback=True)

# copy amp from cd h5parm to iono h5parm
#logger.info('Prepare final globaldb...')
#s.add('H5parm_merge.py cal-cd.h5:sol000 cal-iono.h5:solcd', log='losoto-iono.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
#
#s.add('losoto -v cal-iono.h5 '+parset_dir+'/losoto-dup.parset', log='losoto-iono.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)
#
#s.add('H5parm_exporter.py -v -c --soltab amplitudeSmooth000,phaseOrig000 cal-iono.h5 globaldb', log='losoto-iono.log', log_append=True, cmd_type='python', processors='max')
#s.run(check=True)

if 'survey' in os.getcwd():
    check_rm('globaldb/instrument*') # keep only filled instrument tables
    newglobaldb = 'globaldb_'+os.getcwd().split('/')[-2]
    logger.info('Copy: globaldb -> dsk:/disks/paradata/fdg/LBAsurvey/%s' % newglobaldb)
    os.system('ssh dsk "rm -f /disks/paradata/fdg/LBAsurvey/%s"' % newglobaldb)
    os.system('scp -q -r globaldb dsk:/disks/paradata/fdg/LBAsurvey/%s' % newglobaldb)

# a debug image
if imaging:

    from make_mask import make_mask
    if not 'survey' in os.getcwd():
        mss = mss[int(len(mss)/2.):] # keep only upper band

    # Correct all CORRECTED_DATA (beam, CD, FR, BP corrected) -> CORRECTED_DATA
    # TEST: wegith -> true
    logger.info('Amp/ph correction...')
    for ms in mss:
        s.add('NDPPP '+parset_dir+'/NDPPP-cor.parset msin='+ms+' cor.updateweights=False cor.parmdb='+ms+'/instrument cor.correction=gain', log=ms+'_corG.log', cmd_type='NDPPP')
    s.run(check=True)

#    logger.info('Subtract model...')
#    for ms in mss:
#        s.add('taql "update '+ms+' set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log=ms+'_taql.log', cmd_type='general')
#    s.run(check=True)

    logger.info('Cleaning...')
    check_rm('img')
    os.makedirs('img')
    imagename = 'img/wide'
    s.add('wsclean -reorder -name ' + imagename + ' -size 5000 5000 -trim 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 6arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -mgain 0.9 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 20 '+' '.join(mss), \
            log='wscleanA.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

    # make mask
    maskname = imagename+'-mask.fits'
    make_mask(image_name = imagename+'-MFS-image.fits', mask_name = maskname, threshisl = 3, atrous_do=True)
    # remove CC not in mask
    for modelname in sorted(glob.glob(imagename+'*model.fits')):
        blank_image_fits(modelname, maskname, inverse=True)

    logger.info('Cleaning w/ mask')
    imagename = 'img/wideM'
    s.add('wsclean -reorder -name ' + imagename + ' -size 5000 5000 -trim 4000 4000 -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale 6arcsec -weight briggs 0.0 -niter 100000 -no-update-model-required -mgain 0.8 -minuv-l 100 \
            -pol I -joinchannels -fit-spectral-pol 2 -channelsout 10 -auto-threshold 0.1 -fitsmask '+maskname+' '+' '.join(mss), \
            log='wscleanB.log', cmd_type='wsclean', processors='max')
    s.run(check=True)

logger.info("Done.")
