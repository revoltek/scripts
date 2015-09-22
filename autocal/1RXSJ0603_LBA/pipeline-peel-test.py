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
s = Scheduler(qsub=True, max_threads=50, dry=False, max_processors=6)

# TODO: iterate on DD calibrators

#########################################################################################
# Clear
logging.info('Cleaning...')
check_rm('concat*') 
check_rm('*log *last *pickle')
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
# Add DD cal model - group*_TC*.MS:MODEL_DATA (high+low resolution model)
logging.info('Add DD calibrator...')
for g in groups:
    logging.debug('Working group: '+g)
    model = 'self/models/wide-g'+g+'.model'
    modellr = 'self/models/wide-lr-g'+g+'.model'
    check_rm('concat.MS*')
    pt.msutil.msconcat(sorted(glob.glob('group'+g+'_TC*.MS')), 'concat.MS', concatTime=False)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':model, 'wproj':1024}, log='init-g'+g+'-ft1.log')
    s.run(check=True) # no parallel
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':modellr, 'wproj':1024, 'incr':True}, log='init-g'+g+'-ft2.log')
    s.run(check=True) # no parallel

####################################################################################################
# [PARALLEL] ADD (corrupt) group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + DD cal from model, cirular, beam correcred)
#logging.info('Add and corrupt model...')
#for ms in sorted(glob.glob('group*_TC*.MS')):
#    s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-init_add-test.parset', \
#            log=ms+'_init-addcor.log', cmd_type='BBS')
#s.run(check=True)

#####################################################################################################
logging.info('Correct...')
for ms in sorted(glob.glob('group*_TC*.MS')):
    s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' bbs-subadd.parset', \
            log=ms+'_subadd.log', cmd_type='BBS')
s.run(check=True)

# pre-concat
#logging.info('Averaging...')
#check_rm('peel*MS') 
#for tc in tcs:
#    logging.debug('Time chunk (DATA): '+tc)
#    mss = glob.glob('group*_TC'+tc+'.MS')
#    msout = 'peel-avg_TC'+tc+'.MS'
#    s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-avg.parset msin="['+','.join(mss)+']" msout='+msout+' msin.datacolumn=CORRECTED_DATA \
#           avg.freqstep=4 avg.timestep=10', log=msout+'_init-shiftavg.log', cmd_type='NDPPP')
#s.run(check=True)

# Concatenating (in time) before imaging peel-avgavg_TC*.MS:DATA -> concat.MS:DATA (beam corrected, only source to peel in the data, all chan)
check_rm('concat.MS*')
logging.info('Concatenating TCs...')
pt.msutil.msconcat(sorted(glob.glob('group*_TC*.MS')), 'concat.MS', concatTime=False)
#pt.msutil.msconcat(glob.glob('peel-avg_TC*.MS'), 'concat.MS', concatTime=False)

# Clean mask clean
imagename = 'peel/'+dd['name']+'/images/peel-test'
logging.info('Cleaning 1...')
s.add('wsclean -reorder -name ' + imagename + ' -size 5000 5000 -mem 90 \
        -scale 5arcsec -weight briggs 0.0 -niter 100000 -mgain 0.75 -no-update-model-required -maxuv-l 8000 concat.MS', \
        log='wsclean.log', cmd_type='wsclean')
s.run(check=True)
#s.add_casa('/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', \
#            params={'msfile':'concat.MS', 'imagename':imagename, 'imsize':4096, 'niter':5000, 'multiscale':[0,3,9,18], 'wproj':512}, log='casaclean1.log')
#s.run(check=True)
#logging.info('Cleaning 2...')
#s.add_casa('/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', \
#            params={'msfile':'concat.MS', 'imagename':imagename+'-lr', 'imsize':1024, 'cell':'10arcsec', 'niter':5000, 'multiscale':[0], 'wproj':512, 'uvtaper':True, 'outertaper':'60arcsec'}, log='casaclean2.log')
#s.run(check=True)
logging.info('Done.')
