#!/usr/bin/python
# perform DD-calibration with sagecal and imaging
# Input:
# all DIE-calibrated MSs splitted in SBs
# a sky model extrated during DIE calibration
# Output:
#

model = '~/scripts/autocal/LBAsurvey/toothbrush.LBA.skymodel'
modelname = 'LBA'
#model = '~/scripts/autocal/LBAsurvey/toothbrush.HBA150.skymodel' # add beam (-B 1)
#modelname = 'HBA'
number_of_sbs = 30
number_of_clusters = 10
restoreid = 2 # this is the cluster number that will be resotred
inttime = 120

###############################################################
import os, sys, glob
import numpy as np
from lib_pipeline import *

set_logger()
check_rm('logs')
s = Scheduler(dry=False)
os.makedirs('logs/mss')

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
for tc in xrange(10):
    mss =[]
    for c in np.linspace(0, 121, number_of_sbs).astype(int):
        ms = 'mss/TC0'+str(tc)+'-c'+str(c)+'.MS'
        if not os.path.exists(ms): continue
        mss.append(ms)
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

    # sagecal calibration and empty field in CORRECTED_DATA
    logging.info('Run sagecal -> calibration...')
    s.add('mpirun -np '+str(number_of_sbs+1)+' ~sarod/bin/sagecal-mpi -n 1 -j 5 -f sagecal.mss -I DATA -O CORRECTED_DATA -s sagecal.model -c sagecal.clusters -A 10 -P 3 -Q 3 -r 10 -F 0 -t '+str(inttime), log='sagecal_cal-g'+str(i)+'.log', cmd_type='general', processors='max')
    s.run(check=True)

    # restored cluster in RESTORED_DATA
    logging.info('Run sagecal -> restoring...')
    with open('ignore.clusters', 'w') as ignore_file:
        for idc in xrange(number_of_clusters):
            if idc != restoreid: ignore_file.write("%i\n" % idc)

    for ms in mss:
        s.add('~sarod/bin/sagecal -a 2 -s sagecal.model -c sagecal.clusters -F 0 -I CORRECTED_DATA -O RESTORED_DATA -t '+str(inttime)+' -n 4 -k '+str(restoreid)+' -d '+ms+' -p '+ms+'.solutions -z ignore.clusters', log=ms+'_sagecal-restore.log', cmd_type='general', processors='max')
    s.run(check=True)

check_rm('img')
os.mkdir('img')
img_extname = '_model'+modelname+'_nClusters'+str(number_of_clusters)+'_nSB'+str(number_of_sbs)+'_time'+str(inttime)

# imaging restored facet
logging.info('Imaging -> restored facet...')
s.add('wsclean -reorder -name img/restored'+img_extname+' -size 5000 5000 -mem 100 -j 40 -baseline-averaging 2.0 -scale 5arcsec -weight briggs 0.0 -niter 10000 -no-update-model-required -mgain 0.7 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -datacolumn RESTORED_DATA '+' '.join([item for sublist in all_mss for item in sublist]), log='wsclean-restored.log', cmd_type='wsclean', processors='max')
s.run(check=True)

# imaging empty field
logging.info('Imaging -> empty field...')
s.add('wsclean -reorder -name img/empty'+img_extname+' -size 5000 5000 -mem 100 -j 40 -baseline-averaging 2.0 -scale 5arcsec -weight briggs 0.0 -niter 1 -no-update-model-required -mgain 0.7 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -datacolumn CORRECTED_DATA '+' '.join([item for sublist in all_mss for item in sublist]), log='wsclean-empty.log', cmd_type='wsclean', processors='max')
s.run(check=True)

# imaging orig
logging.info('Imaging -> original...')
s.add('wsclean -reorder -name img/orig'+img_extname+' -size 5000 5000 -mem 100 -j 40 -baseline-averaging 2.0 -scale 5arcsec -weight briggs 0.0 -niter 10000 -no-update-model-required -mgain 0.7 -pol I -cleanborder 0 -joinchannels -fit-spectral-pol 2 -channelsout 10 -deconvolution-channels 5 -datacolumn DATA '+' '.join([item for sublist in all_mss for item in sublist]), log='wsclean-orig.log', cmd_type='wsclean', processors='max')
s.run(check=True)
