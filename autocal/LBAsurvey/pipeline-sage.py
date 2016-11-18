#!/usr/bin/python
# perform DD-calibration with sagecal and imaging
# Input:
# all DIE-calibrated MSs splitted in SBs
# a sky model extrated during DIE calibration
# Output:
#

parset_dir = '/home/fdg/scripts/autocal/LBAsurvey/parset_sage'
concat_mss = '/data/scratch/fdg/tooth/tgts10h-selfcal-half/mss'
model = '~/scripts/autocal/LBAsurvey/toothbrush.LBA.skymodel'
modelname = 'LBA'
#model = '~/scripts/autocal/LBAsurvey/toothbrush.HBA150.skymodel' # add beam (-B 1)
#modelname = 'HBA'
number_of_sbs = 61 # 61 to use all
number_of_chan_per_sb = 8 # 8 to use all
number_of_clusters = 10
restoreid = 2 # this is the cluster number that will be resotred
inttime = 40 # default 120

###############################################################
import os, sys, glob, re
import numpy as np
from lib_pipeline import *

set_logger()
check_rm('logs')
s = Scheduler(dry=False)

img_extname = '_model'+modelname+'_nClusters'+str(number_of_clusters)+'_nSB'+str(number_of_sbs)+'_time'+str(inttime)

##############################################################
# Split mss
# TODO: reduce SB size, freq variation should be mild
logging.info('Splitting MSs...')
for ms in sorted(glob.glob(concat_mss+'/TC*MS')):
    tc = int(re.findall(r'\d+', ms)[-1])
    logging.info('Working on time chunk: %i' % tc)
    tcdir = 'split/tc%02i/' % tc
    if not os.path.exists(tcdir): os.makedirs(tcdir)
    if not os.path.exists('logs/'+tcdir): os.makedirs('logs/'+tcdir)
    nch = find_nchan(ms)
    assert nch % number_of_chan_per_sb == 0 # channels must be divisible by number_of_chan_per_sb
    for ch in xrange(int(nch/number_of_chan_per_sb)):
        msout = tcdir+'/'+ms.split('/')[-1].replace('.MS','-c%03i.MS' % ch)
        if os.path.exists(msout): continue
        s.add('NDPPP '+parset_dir+'/NDPPP-split.parset msin='+ms+' msin.startchan='+str(ch*number_of_chan_per_sb)+' msin.nchan='+str(number_of_chan_per_sb)+' msout='+msout, \
                log=msout+'_split.log', cmd_type='NDPPP')
    s.run(check=True, max_threads=10)

##############################################################
# prepare model
logging.info('Prepare model file...')
check_rm('sagecal.model')
s.add('sagecal_convert_skymodel.py -i '+model+' -o sagecal.model -b', log='convert_skymodel.log', cmd_type='python')
s.run(check=True)

logging.info('Prepare cluster file...')
check_rm('sagecal.clusters')
s.add('sagecal_create_clusters.py -s sagecal.model -c '+str(number_of_clusters)+' -o sagecal.clusters -i 10', log='create_clusters.log', cmd_type='python')
s.run(check=True)
s.add('sagecal_annotate.py -s sagecal.model -c sagecal.clusters -o '+modelname+'_'+str(number_of_clusters)+'clusters.reg', log='annotate.log', cmd_type='python')
s.run(check=True)

logging.info('Preparing file list...')
all_mss=[]
for tc in sorted(glob.glob('split/tc*')):
    mss = sorted(glob.glob(tc+'/*MS'))
    mss = np.array(mss)[np.linspace(0, len(mss)-1, number_of_sbs).astype(int)].tolist()
    all_mss.append(mss)

for i, mss in enumerate(all_mss):
    logging.info('Cycle on hour '+str(i))

    # prepare sagecal file with list of mss
    check_rm('sagecal.mss')
    with open('sagecal.mss', 'w') as mss_file:
        for ms in mss:
            mss_file.write("%s\n" % ms)

    # add necessary columns
    logging.info('Add columns...')
    for ms in mss:
        s.add('addcol2ms.py -m '+ms+' -c CORRECTED_DATA,RESTORED_DATA', log=ms+'_addcol.log', cmd_type='python')
    s.run(check=True)

    # TODO: test smoothing
    #logging.info('BL-based smoothing...')
    #for ms in mss:
    #    s.add('BLavg.py -r -w -i DATA -o SMOOTHED_DATA '+ms, log=ms+'_smooth.log', cmd_type='python', processors='max')
    #s.run(check=True, max_threads=10)

    # -J 1 -> only phase
    # -r 10/50 -> increase regularization?
    # -y 14000 -> cut some of the longest BLs

    # sagecal calibration and empty field in CORRECTED_DATA
    logging.info('Run sagecal -> calibration...')
    s.add('mpirun -np '+str(number_of_sbs+1)+' ~sarod/bin/sagecal-mpi -n 1 -j 5 -f sagecal.mss -I DATA -O CORRECTED_DATA -s sagecal.model -p solutions-TC'+str(i)+'.txt -F 0 -c sagecal.clusters -A 10 -P 3 -Q 3 -r 10 -J 1 -t '+str(inttime), log='sagecal_cal-g'+str(i)+'.log', cmd_type='general', processors='max')
    s.run(check=True)

    # restored cluster in RESTORED_DATA
    logging.info('Run sagecal -> restoring...')
    with open('ignore.clusters', 'w') as ignore_file:
        for idc in xrange(number_of_clusters):
            if idc != restoreid: ignore_file.write("%i\n" % idc)

    for ms in mss:
        s.add('~sarod/bin/sagecal -n 4 -a 2 -d '+ms+' -I CORRECTED_DATA -O RESTORED_DATA -s sagecal.model -F 0 -c sagecal.clusters -t '+str(inttime)+' -k '+str(restoreid)+' -p '+ms+'.solutions -z ignore.clusters -J 1', log=ms+'_sagecal-restore.log', cmd_type='general', processors='max')
    s.run(check=True)

    # TODO: flagging of spikes

check_rm('img')
os.mkdir('img')

# imaging orig
logging.info('Imaging -> original...')
s.add('wsclean -reorder -name img/orig'+img_extname+' -size 5000 5000 -mem 100 -j 40 -baseline-averaging 2.0 -scale 5arcsec -weight briggs 0.0 -niter 10000 -no-update-model-required -mgain 0.7 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -datacolumn DATA '+' '.join([item for sublist in all_mss for item in sublist]), log='wsclean-orig.log', cmd_type='wsclean', processors='max')
s.run(check=True)

# imaging restored facet
logging.info('Imaging -> restored facet...')
s.add('wsclean -reorder -name img/restored'+img_extname+' -size 5000 5000 -mem 100 -j 40 -baseline-averaging 2.0 -scale 5arcsec -weight briggs 0.0 -niter 10000 -no-update-model-required -mgain 0.7 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -datacolumn RESTORED_DATA '+' '.join([item for sublist in all_mss for item in sublist]), log='wsclean-restored.log', cmd_type='wsclean', processors='max')
s.run(check=True)

# imaging empty field
logging.info('Imaging -> empty field...')
s.add('wsclean -reorder -name img/empty'+img_extname+' -size 5000 5000 -mem 100 -j 40 -baseline-averaging 2.0 -scale 5arcsec -weight briggs 0.0 -niter 1 -no-update-model-required -mgain 0.7 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -datacolumn CORRECTED_DATA '+' '.join([item for sublist in all_mss for item in sublist]), log='wsclean-empty.log', cmd_type='wsclean', processors='max')
s.run(check=True)
