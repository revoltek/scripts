#!/usr/bin/env python

import os, sys, glob
from LiLF import lib_ms, lib_img, lib_util, lib_log

cycle='09'

# make dirs
#for obs in glob.glob('c*-o*'):
#    for tgt in glob.glob(obs+'/p*'):
#        if not os.path.isdir('c'+cycle+'-pointings/%s' % (os.path.basename(tgt))):
#            print tgt
#            os.system('mkdir c'+cycle+'-pointings/%s' % (os.path.basename(tgt)))

# run calibrator pipeline
#for cal in glob.glob('c*-o*/3c*'):
#    if not os.path.exists('/home/fdg/data/LBAsurvey/%s' % cal):
#        os.system('mkdir ~/data/LBAsurvey/%s' % cal)
#    os.chdir('/home/fdg/data/LBAsurvey/%s' % cal)
#    print ("-------- Working on: %s" % cal)
#    if not os.path.exists('data-bkp'):
#        os.system('mkdir data-bkp')
#        os.system('cp -r /lofar1/LBAsurvey/%s/*MS data-bkp' % cal)
#    os.system('rm *walker')
#    os.system('~/LiLF/pipelines/LOFAR_cal.py')

# Make links
#for cal in glob.glob('c*-o*/3c*'):
#    print ("-------- Working on: %s" % cal)
#    MSs = lib_ms.AllMSs( glob.glob(cal+'/data-bkp/*MS'), None, check_flags=False)
#    obsname = cal.split('/')[0]
#    calname = cal.split('/')[1]
#    for MS in MSs.getListObj():
#        obsid = MS.getObsID()
#        if not os.path.exists('%s/id%s_-_%s' % (obsname,obsid,calname)):
#            os.system('ln -s %s %s/id%s_-_%s' % (calname,obsname,obsid,calname))

# run timesplit pipeline
#for tgt in glob.glob('c*-o*/p*'):
#    os.chdir('/home/fdg/data/LBAsurvey/%s' % tgt)
#    print ("-------- Working on: %s" % tgt)
#    if not os.path.exists('data-bkp'):
#        os.system('mkdir data-bkp')
#        os.system('cp -r /lofar1/LBAsurvey/%s/*MS data-bkp' % tgt)
#    os.system('rm *walker')
#    os.system('~/LiLF/pipelines/LOFAR_timesplit.py')

# copy data
#for tgt in sorted(glob.glob('c*-pointings/p*')):
#    tgt_name = tgt.split('/')[1]
#    os.chdir('/home/fdg/data/LBAsurvey/%s' % tgt)
#    print ("-------- Working on: %s" % tgt)
#    if not os.path.exists('OLD-mss'):
#        os.system('mv mss OLD-mss')
#    os.system('mkdir mss')
#    for i, ms in enumerate(glob.glob('../../c*-o*/%s/mss/*' % tgt_name)):
#        os.system('cp -r %s mss/TC%02i.MS' % (ms,i))

# run self pipeline
#for tgt in sorted(glob.glob('c*-pointings/p*')):
#    os.chdir('/home/fdg/data/LBAsurvey/%s' % tgt)
#    print ("-------- Working on: %s" % tgt)
#    os.system('~/LiLF/pipelines/LOFAR_self.py')

# run dd pipeline
#for tgt in sorted(glob.glob('*-pointings/p*')):
#    os.chdir('/home/fdg/data/LBAsurvey/%s' % tgt)
#    print ("-------- Working on: %s" % tgt)
#    os.system('~/LiLF/pipelines/LOFAR_dd-serial.py')

# special run
for tgt in sorted(glob.glob('*-pointings/p*')):
    os.chdir('/home/fdg/data/LBAsurvey/%s' % tgt)
    print ("-------- Working on: %s" % tgt)
    mss = ','.join(glob.glob('mss-avg/*MS'))
    # DATA is with the peeled sources
    cmd = 'DDF.py --Debug-Pdb=never --Parallel-NCPU=64  --Cache-Dir . --Data-MS '+mss+' --Data-ColName DATA --Data-Sort 1 --Output-Mode Clean --Deconv-CycleFactor 0 --Deconv-MaxMinorIter 1000000 --Deconv-RMSFactor 2.0 --Deconv-FluxThreshold 0.0 --Deconv-Mode HMP --HMP-AllowResidIncrease 1.0 --Weight-Robust -0.5 --Image-NPix 8794 --CF-wmax 50000 --CF-Nw 100 --Beam-CenterNorm 1 --Beam-Smooth 1 --Beam-Model LOFAR --Beam-LOFARBeamMode A --Beam-NBand 1 --Beam-DtBeamMin 5 --Output-Also onNeds --Image-Cell 3.0 --Freq-NDegridBand 6 --Freq-NBand 6 --Mask-Auto 1 --Mask-SigTh 5.0 --GAClean-MinSizeInit 10 --GAClean-MaxMinorIterInitHMP 100000 --Facets-DiamMax 1.5 --Facets-DiamMin 0.1 --Weight-ColName WEIGHT_SPECTRUM --Output-Name img/wideDD-init --RIME-ForwardMode BDA-degrid --Output-RestoringBeam 15.0 --Deconv-MaxMajorIter 10 --Deconv-PeakFactor 0.005 --Mask-External ddcal/c01/images/wideDD-c01-mask.fits --Cache-Reset 0 >> logs/ddfacetM-init.log  2>&1'
    os.system(cmd)
