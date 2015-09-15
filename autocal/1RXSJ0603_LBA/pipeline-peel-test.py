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
    model = 'self/models/wide-g'+g+'model*'
    model = 'self/models/wide-lr-g'+g+'model*'
    check_rm('concat.MS*')
    pt.msutil.msconcat(sorted(glob.glob('group'+g+'_TC*.MS')), 'concat.MS', concatTime=False)
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':model}, log='init-g'+g+'-ft1.log')
    s.run(check=True) # no parallel
    s.add_casa('/home/fdg/scripts/autocal/casa_comm/casa_ft.py', params={'msfile':'concat.MS', 'model':modellr, 'incr':True}, log='init-g'+g+'-ft2.log')
    s.run(check=True) # no parallel

# [PARALLEL] ADD (corrupt) group*_TC*.MS:SUBTRACTED_DATA + MODEL_DATA -> group*_TC*.MS:CORRECTED_DATA (empty data + DD cal from model, cirular, beam correcred)
logging.info('Add and corrupt model...')
for ms in sorted(glob.glob('group*_TC*.MS')):
    s.add('calibrate-stand-alone --parmdb-name instrument '+ms+' /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/bbs-init_add-test.parset', \
            log=ms+'_init-addcor.log', cmd_type='BBS')
s.run(check=True)

# [PARALLEL] averaging before cleaning peel-avg_TC*.MS:CORRECTED_DATA -> peel-avgavg_TC*.MS:DATA
check_rm('peel-avgavg_TC*.MS')
logging.info('Averaging before cleaning...')
nchan = find_nchan('group8_TC0.MS')
for ms in glob.glob('group*_TC*.MS'):
    msout = ms.replace('.MS','avg.MS')
    s.add('NDPPP /home/fdg/scripts/autocal/1RXSJ0603_LBA/parset_peel/NDPPP-avg.parset msin='+ms+' msin.nchan='+str(nchan-nchan%4)+' msin.datacolumn=CORRECTED_DATA \
                msout='+msout+' avg.freqstep=4 avg.timestep=10', log=ms+'_avgclean-c'+str(i)+'.log', cmd_type='NDPPP')
s.run(check=True)

# Concatenating (in time) before imaging peel-avgavg_TC*.MS:DATA -> concat.MS:DATA (beam corrected, only source to peel in the data, all chan)
check_rm('concat.MS*')
logging.info('Concatenating TCs...')
pt.msutil.msconcat(glob.glob('group*_TC*.MS'), 'concat.MS', concatTime=False)

# Clean mask clean
imagename = 'peel/'+dd['name']+'/images/peel-'+str(i)
imsize = 4096
logging.debug('Image size set to '+str(imsize))
logging.info('Cleaning...')
s.add_casa('/home/fdg/scripts/autocal/casa_comm/1RXSJ0603_LBA/casa_clean_peel.py', \
            params={'msfile':'concat.MS', 'imagename':imagename, 'imsize':imsize, 'niter':50000, 'multiscale':[0], 'wproj':512}, log='casaclean1-c'+str(i)+'.log')
s.run(check=True)
