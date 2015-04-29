#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys, time, glob
from subprocess import Popen, PIPE
import numpy
import pyfits
import pyrap.tables as pt
import pyrap.images
import lofar.parmdb
import lofar.bdsm as bdsm
 
# import psutil
pi = numpy.pi

# to run casa when not logged in
# Step 0. Run in screen
# Step 1. Xvfb :5 &
# Step 2. setenv DISPLAY :5
# run the program in screen and to exit screen: screen -d

# parameters
msfile = 'tooth-concat.MS'
cluster = 'block6'
atrous_do = True
imsize = 4095
skymodel = 'toothbrush.GMRT150.skymodel'
nterms = 1
cellsizetime_a = 600
cellsizetime_p = 1
TEC = True

wplanes = 1
if imsize > 512:
    wplanes = 64
if imsize > 799:
    wplanes = 96
if imsize > 1023:
    wplanes = 128
if imsize > 1599:
    wplanes = 196
if imsize > 2049:
    wplanes = 256
if imsize > 3000:
    wplanes = 448
if imsize > 4095:
    wplanes = 512


def thread_cmd(cmd_list, ncpu=8):
    """
    Run a sequence of commands in parallell
    cmd_list: list of commands
    ncpu: max number of processes at the same time
    """
    import subprocess
    import time

    processes = set()

    for cmd in cmd_list:
        print "RUN: "+cmd
        processes.add(subprocess.Popen(cmd, shell=True))
        while len(processes) >= ncpu:
            # wait for the child process to finish
            time.sleep(5)
            # remove finished process
            processes.difference_update(set(p for p in processes if p.poll() is not None))

    #Check if all the child processes were closed
    for p in processes:
        if p.poll() is None:
            p.wait()

    return 0


def runpybdsm(image_name, threshpix=5, threshisl=3, atrous_do=False, do_mask=True):
    """
    Run source extractor and optionally
    Create a clean mask
    """

    # convert to boolean
    if atrous_do == True:
        print 'Changing island threshold to 4 because atrous_do=True'
        threshisl = 4.0

    soumodel = 'skymodels/'+image_name.split('.image')[0] + '.skymodel'
    os.system('rm ' + soumodel)

    # DO THE SOURCE DETECTION
    img = bdsm.process_image(
        image_name,
        mean_map='zero',
        rms_box=(70, 10),
        thresh_pix=int(threshpix),
        thresh_isl=int(threshisl),
        atrous_do=atrous_do,
        ini_method='curvature')

    img.export_image(img_type='sou_model', outfile=soumodel)

    # WRITE THE GAUSSIAN MODEL FITS
    if do_mask:
        mask_name = 'masks/'+image_name.split('.image')[0] + '.cleanmask'
        os.system('rm -r ' + mask_name)
        os.system('cp -r ' + image_name + ' ' + mask_name)

        gausmodel = 'skymodels/'+image_name.split('.image')[0] + '.gausmodel'
        os.system('rm ' + gausmodel)

        img.export_image(img_type='gaus_model', outfile=gausmodel)

        print 'Making mask:', mask_name
        img = pyrap.images.image(mask_name)
        pixels = numpy.copy(img.getdata())
        pixels_mask = 0. * numpy.copy(pixels)

        hdulist = pyfits.open(gausmodel)
        pixels_gs = hdulist[0].data

        gs_cut = 1e-3
        idx = numpy.where(pixels_gs > gs_cut)
        pixels_mask[idx] = 1.0

        img.putdata(pixels_mask)
        
        return soumodel, mask_name
    return soumodel


def make_image(
    msfile,
    cluster,
    callnumber,
    threshpix,
    threshisl,
    imsize,
    wplanes,
    cellsize = 4, # in arcsec
    nterms = 1,
    niter = 1000, # 7500 causes nasty clean artifacts in HBA
    uvmax = 10000, # in klambda
    robust = 0,
    casa = False, # use casa o awimager?
    atrous_do = True,
    do_mask = True,
    do_sourceext = True):
    """
    Create an image using casa_clean.py or awimager, the mask is created on the fly

    return:
    the image, the mask ('' if not), and the skymodel ('' if not)
    """

    mscale = 'False'
    if atrous_do:
        mscale = 'True'

    # do the masking
    skymodel = ''
    mask = ''
    if do_mask:

        imout = 'img/im' + callnumber + '_cluster' + cluster + '-nm'
        if casa: 
            os.system('casapy --nologger -c casa_clean.py ' + msfile + '  ' + imout + ' ' + 'None' + ' ' + '1mJy'
                  + ' ' + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale + ' ' + wplanes)
        else:
            cmd = 'awimager ms=' + msfile + ' image=' + imout + ' weight=briggs npix=8192 robust=' + str(robust) + ' cellsize=' + str(cellsize) + 'arcsec data=CORRECTED_DATA padding=1. niter=75000 stokes=I threshold=0. operation=mfclark oversample=8 timewindow=36000 wmax=7000 cyclefactor=1.75 gain=0.1 PBCut=0.005 UseLIG=1 UseWSplit=1 ApplyBeamCode=3 ApplyElement=0 TWElement=14 SpheSupport=15 MakeDirtyCorr=0 SingleGridMode=1 FindNWplanes=1 UVmin=0.08 UVmax=' + str(uvmax) + ' ChanBlockSize=1'
            os.system(cmd)

        if nterms > 1:
            _null, mask = runpybdsm(imout + '.image.tt0', threshpix, threshisl, atrous_do, do_mask=True)
        else:
            _null, mask = runpybdsm(imout + '.image', threshpix, threshisl, atrous_do, do_mask=True)

    # do the cleaning
    print 'IMAGING: ' + imout + ' ' + 'None' + ' ' + '1mJy' + ' ' + str(niter) + ' ' + str(nterms)
    imout = 'img/im' + callnumber + '_cluster' + cluster

    if casa: 
        os.system('casapy --nologger -c casa_clean.py ' + msfile + ' ' + imout + ' ' + mask + ' ' + '1mJy' + ' '
                  + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' ' + mscale + ' ' + wplanes)
    else:
        cmd = 'awimager ms=' + msfile + ' image=' + imout + ' weight=briggs npix=8192 robust=' + str(robust) + ' cellsize=' + str(cellsize) + 'arcsec data=CORRECTED_DATA padding=1. niter=75000 stokes=I threshold=0. operation=mfclark oversample=8 timewindow=36000 wmax=7000 cyclefactor=1.75 gain=0.1 PBCut=0.005 UseLIG=1 UseWSplit=1 ApplyBeamCode=3 ApplyElement=0 TWElement=14 SpheSupport=15 MakeDirtyCorr=0 SingleGridMode=1 FindNWplanes=1 UVmin=0.08 UVmax=' + str(uvmax) + ' ChanBlockSize=1'
        os.system(cmd)

    # convert to FITS and extract skymodel
    skymodel = ''
    if nterms > 1:
        os.system('image2fits in=' + imout + '.image.tt0' + ' ' + 'out=' + imout + '.fits')
        if sourceext:
            skymodel = runpybdsm(imout+'.image.tt0', threshpix, threshisl, atrous_do, do_mask=False)
    else:
        os.system('image2fits in=' + imout + '.image' + ' ' + 'out=' + imout + '.fits')
        if sourceext:
            skymodel = runpybdsm(imout+'.image', threshpix, threshisl, atrous_do, do_mask=False)

    return imout, mask, skymodel


def runbbs(msfile, skymodel, parset, parmdb='instrument', applycal=False):
    """
    Run BBS
    applycal: if true do not use the -f in calibrate command
    """
    
    log = 'log/'+msfile+'.log'

    if applycal:
        cmd = 'calibrate-stand-alone -t 8 -f --parmdb-name ' + parmdb + ' ' + msfile + ' ' + parset + ' ' + skymodel + '>>' + log
        #cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb + ' ' + msfile + ' ' + parset + ' ' + skymodel + '>>' + log
    else:
        cmd = 'calibrate-stand-alone -t 8 --parmdb-name ' + parmdb + ' ' + msfile + ' ' + parset + ' ' + skymodel + '>>' + log
        #cmd = 'calibrate-stand-alone -f --parmdb-name ' + parmdb + ' ' + msfile + ' ' + parset + ' ' + skymodel + '>>' + log

    os.system(cmd)
    return 0


###############################
#### APPLY CAL AMP SOLUTIONS
print "Applyting cal amp solutions"
#os.system('~/scripts/losoto/tools/H5parm_importer.py -v tooth_init.h5  globaldb-cal >> log/losoto_importer.log')
#os.system('~/scripts/losoto/losoto.py -v tooth_init.h5 parsets/losoto-reset.parset >> log/losoto.log')
#os.system('~/scripts/losoto/tools/H5parm_exporter.py -v -c tooth_init.h5 globaldb-cal >> log/losoto_exporter.log')

#for j, SBfile in enumerate(glob.glob('data/L*MS')):
#    os.system('rm -r '+SBfile+'/instrument')
#    os.system('cp -r globaldb-cal/sol000_instrument-'+str(j)+' '+SBfile+'/instrument')

#cmd_list = []
#for j, SBfile in enumerate(glob.glob('data/L*MS')):
#    cmd_list.append('calibrate-stand-alone '+SBfile+' parsets/bbs-init_corr.parset toothbrush.fakemodel.skymodel >> log/bbs-init_corr.log')
#thread_cmd(cmd_list)

#os.system('NDPPP parsets/NDPPP-concat.parset')

#sys.exit(0)

###############################
#### CALIB ON PH on GMRT model + IMAGING
print "Starting cycle 0"
skymodel = 'toothbrush.GMRT150.skymodel'
parset = 'parsets/bbs-global_ph.parset'

runbbs(msfile, skymodel, parset, 'instrument', False)
imout, mask, skymodel = make_image(msfile, cluster, callnumber='0', threshpix=10, threshisl=6, imsize=imsize, cellsize=4., wplanes=wplanes, nterms=nterms, atrous_do=atrous_do)

sys.exit(0)

###############################
### CALIB ON AMP AND PH
print "Starting cycle 1"

# create skymodel for BBS
#os.system('casapy2bbs.py -m ' + mask + ' '
#           + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout + '.skymodel')

# FFT
#os.system('casapy --nologger -c ft.py '
#               + msinputlist + ' ' + imout + '.model' + ' '
#              + str(nterms) + ' ' + str(wplanes))

# phase+TEC at 10 s
runbbs(msfile, skymodel, 'parsets/bbs-global_fast.parset', 'instrument-fast', False)
os.system('~/scripts/losoto/tools/H5parm_importer.py -v -s fast -i instrument-fast '+msfile+'1.h5 '+msfile+' >> log/losoto_importer.log')
# DEBUG imaging
__null, __null, __null = make_image(msfile, cluster, callnumber='1high', threshpix=10, threshisl=6, imsize=imsize, cellsize=4.,  wplanes=wplanes, nterms=nterms, atrous_do=atrous_do)
# amp (gain) on 10 min
runbbs(msfile, skymodel, 'parsets/bbs-global_slow.parset', 'instrument-slow', False)
os.system('~/scripts/losoto/tools/H5parm_importer.py -v -s slow -i instrument-slow '+msfile+'1.h5 '+msfile+' >> log/losoto_importer.log')
os.system('~/scripts/losoto/losoto.py -v '+msfile+'1.h5 parsets/losoto-normalize.parset >> log/losoto.log')
# losoto to combine
os.system('~/scripts/losoto/losoto.py -v '+msfile+'1.h5 parsets/losoto-merge.parset >> log/losoto.log')
os.system('~/scripts/losoto/tools/H5parm_exporter.py -v -s comb -c '+msfile+'1.h5 '+msfile+' >> log/losoto_exporter.log')
# applycal
runbbs(msfile, skymodel, 'parsets/bbs-global_corr.parset', 'instrument', True)
# imaging
imout, mask, skymodel_high = make_image(msfile, cluster, callnumber='1high', threshpix=10, threshisl=6, imsize=imsize, cellsize=4., wplanes=wplanes, nterms=nterms, atrous_do=atrous_do)
# pybdsm + sub
runbbs(msfile, skymodel_high, 'parsets/bbs-global_sub_corr.parset', 'instrument', True)
# imaging large scales/field (res ~ 1')
imout, mask, skymodel_low = make_image(msfile, cluster, callnumber='1low', threshpix=10, threshisl=6, imsize=imsize, cellsize=20, uvmax=20, robust=1, wplanes=wplanes, nterms=nterms, atrous_do=atrous_do)
# merge
os.system('cat '+skymodel_high+' > skymodels/'+msfile+'1comb.skymodel')
os.system('cat '+skymodel_low+' | grep -v # >> skymodels/'+msfile+'1comb.skymodel')
skymodel = 'skymodels/'+msfile+'1comb.skymodel'

sys.exit(1)


























####################
### CALIBRATE WITH BBS PHASE ONLY 2 ###
# create skymodel for BBS

os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')
if FFT:
    os.system('casapy --nologger -c ft.py '
               + msinputlist + ' ' + imout + '.model' + ' '
              + str(nterms) + ' ' + str(wplanes))

# phase only calibrate

skymodel = imout + '.skymodel'
if len(mslist) == 34:
    parset = create_scalarphase_parset(cellsizetime_p, TEC, '11,11,12',
            FFT)
if len(mslist) == 16:
    parset = create_scalarphase_parset(cellsizetime_p, TEC, '8,8', FFT)
else:
    parset = create_scalarphase_parset(cellsizetime_p, TEC,
            str(len(mslist)), FFT)
runbbs(  # NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)
    mslist,
    skymodel,
    parset,
    'instrument',
    False,
    TEC,
    )

######################################
### MAKE IMAGE 2 ###

imout, mask = make_image(datadir,cluster,'2',10,10,nterms,atrous_do,imsize,wplane)

####################
### CALIBRATE WITH BBS PHASE+AMP 1 ###

os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')
if FFT:
    os.system('casapy --nologger -c ft.py '
               + msinputlist + ' ' + imout + '.model' + ' '
              + str(nterms) + ' ' + str(wplanes))

skymodel = imout + '.skymodel'
if len(mslist) == 34:
    parset = create_scalarphase_parset_p(cellsizetime_p, TEC, '11,11,12'
            , FFT)
if len(mslist) == 16:
    parset = create_scalarphase_parset_p(cellsizetime_p, TEC, '8,8',
            FFT)
else:
    parset = create_scalarphase_parset_p(cellsizetime_p, TEC,
            str(len(mslist)), FFT)

# solve +apply phases

runbbs(
    mslist,
    skymodel,
    parset,
    'instrument_phase0',
    False,
    TEC,
    )

# solve amps

parmdb = 'instrument_amps0'
parset = create_amponly_parset(cellsizetime_a, FFT)
runbbs(
    mslist,
    skymodel,
    parset,
    parmdb,
    False,
    False,
    )

for ms in mslist:

  # remove outliers from the solutions

    if phasors:
        os.system('python /home/weeren/scripts/rx42_hba/smoothcal_rx42.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')
    else:
        os.system('python /home/weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

if smooth:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb + '_smoothed',
        True,
        False,
        )
else:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb,
        True,
        False,
        )

### MAKE IMAGE 3 ###

imout, mask = make_image(datadir,cluster,'3',10,10,nterms,atrous_do,imsize,wplane)

#### CALIBRATE  BBS PHASE+AMP 2 ###
# make model

os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')
if FFT:
    os.system('casapy --nologger -c ft.py '
               + msinputlist + ' ' + imout + '.model' + ' '
              + str(nterms) + ' ' + str(wplanes))

# parmdb keep from previous step

skymodel = imout + '.skymodel'

# pre-apply amp solutions

if TEC:
    if smooth:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p_TEC.parset'
                ,
            parmdb + '_smoothed',
            True,
            False,
            )
    else:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p_TEC.parset'
                ,
            parmdb,
            True,
            False,
            )
else:
    if smooth:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset'
                ,
            parmdb + '_smoothed',
            True,
            False,
            )
    else:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset'
                ,
            parmdb,
            True,
            False,
            )
if smooth:
    pre_apply_parmdb = parmdb + '_smoothed'  # need in merged if this is last calibration round
else:
    pre_apply_parmdb = parmdb  # need in merged if this is last calibration round

# phase only cal

skymodel = imout + '.skymodel'
if len(mslist) == 34:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
            '11,11,12', FFT)
if len(mslist) == 16:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC, '8,8',
            FFT)
else:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
            str(len(mslist)), FFT)
runbbs(
    mslist,
    skymodel,
    parset,
    'instrument_phase1',
    False,
    TEC,
    )

# solve amps

parmdb = 'instrument_amps1'
parset = create_amponly_parset(cellsizetime_a, FFT)
runbbs(
    mslist,
    skymodel,
    parset,
    parmdb,
    False,
    False,
    )

for ms in mslist:

  # remove outliers from the solutions

    if phasors:
        os.system('python /home/weeren/scripts/rx42_hba/smoothcal_rx42.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')
    else:
        os.system('python /home/weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

if smooth:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb + '_smoothed',
        True,
        False,
        )
else:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb,
        True,
        False,
        )

### MAKE IMAGE 4 ###

imout, mask = make_image(datadir,cluster,'4',10,10,nterms,atrous_do,imsize,wplane)

### CREATE FINAL MODEL ###

skymodelf = 'im_cluster' + cluster + '.final.skymodel'
os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + skymodelf)
if FFT:
    os.system('casapy --nologger -c ft.py '
               + msinputlist + ' ' + imout + '.model' + ' '
              + str(nterms) + ' ' + str(wplanes))

### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
### INCLUDES SPLINE INTERPOLARION OF AMPS ###

if merge_parmdb:

    if phasors:
        dummyparset = \
            '/home/weeren/scripts/rx42_hba/scalarphase+amp.parset'
    else:
        if TEC:
            dummyparset = \
                '/home/weeren/scripts/rx42_hba/scalarphase+ap+TEC.parset'
        else:
            dummyparset = \
                '/home/weeren/scripts/rx42_hba/scalarphase+ap.parset'

    dummyparmdb = 'instrument_template'

   # runbbs(mslist, skymodel,dummyparset, dummyparmdb, True, False)
   # TO SPEED THINGS UP, hard coded for RX42

    for ms in mslist:
        os.system('rm -rf ' + ms + '/' + dummyparmdb)
        os.system('cp -r /home/weeren/scripts/rx42_hba/' + dummyparmdb
                  + ' ' + ms + '/')

    if smooth:
        parmdb_a = 'instrument_amps1_smoothed'  # last/best ampplitude(+phase) parmdb
    else:
        parmdb_a = 'instrument_amps1'  # last/best ampplitude(+phase) parmdb
    parmdb_p = 'instrument_phase1'  # last/best CommonScalarPhase parmdb
    parmdbout = 'instrument_merged'

    for ms in mslist:
        create_merged_parmdb(
            ms,
            ms + '/' + pre_apply_parmdb,
            ms + '/' + parmdb_a,
            ms + '/' + parmdb_p,
            ms + '/' + dummyparmdb,
            ms + '/' + parmdbout,
            cellsizetime_a,
            cellsizetime_p,
            )

##################################################################
##################################################################
##################################################################
##################################################################

sys.exit()

### CALIBRATE BBS PHASE+AMP 3 ###
# create skymodel

os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')
if FFT:
    os.system('casapy --nologger -c ft.py '
               + msinputlist + ' ' + imout + '.model' + ' '
              + str(nterms) + ' ' + str(wplanes))

# parmdb keep from previous step

skymodel = imout + '.skymodel'

# pre-apply amp solutions

if TEC:
    if smooth:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p_TEC.parset'
                ,
            parmdb + '_smoothed',
            True,
            False,
            )
    else:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p_TEC.parset'
                ,
            parmdb,
            True,
            False,
            )
else:
    if smooth:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset'
                ,
            parmdb + '_smoothed',
            True,
            False,
            )
    else:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset'
                ,
            parmdb,
            True,
            False,
            )

if smooth:
    pre_apply_parmdb = parmdb + '_smoothed'  # need in merged if this is last calibration round
else:
    pre_apply_parmdb = parmdb  # need in merged if this is last calibration round

# solve +apply phases

if len(mslist) == 34:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
            '11,11,12', FFT)
if len(mslist) == 16:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC, '8,8',
            FFT)
else:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
            str(len(mslist)), FFT)
runbbs(
    mslist,
    skymodel,
    parset,
    'instrument_phase2',
    False,
    TEC,
    )

# solve amps

parmdb = 'instrument_amps2'
parset = create_amponly_parset(cellsizetime_a, FFT)
runbbs(
    mslist,
    skymodel,
    parset,
    parmdb,
    False,
    False,
    )

for ms in mslist:

  # remove outliers from the solutions

    if phasors:
        os.system('python /home/weeren/scripts/rx42_hba/smoothcal_rx42.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')
    else:
        os.system('python /home/weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

if smooth:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb + '_smoothed',
        True,
        False,
        )
else:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb,
        True,
        False,
        )

### MAKE IMAGE 5 ###

imout, mask = make_image(datadir,cluster,'5',8,8,nterms,atrous_do,imsize,wplane)

### CREATE FINAL MODEL ###

skymodelf = 'im_cluster' + cluster + '.final.skymodel'
os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + skymodelf)
if FFT:
    os.system('casapy --nologger -c ft.py '
               + msinputlist + ' ' + imout + '.model' + ' '
              + str(nterms) + ' ' + str(wplanes))

### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
### INCLUDES SPLINE INTERPOLARION OF AMPS ###

if merge_parmdb:

    if phasors:
        dummyparset = \
            '/home/weeren/scripts/rx42_hba/scalarphase+amp.parset'
    else:
        if TEC:
            if FFT:
                dummyparset = \
                    '/home/weeren/scripts/rx42_hba/scalarphase+ap+TEC+FFT.parset'
            else:
                dummyparset = \
                    '/home/weeren/scripts/rx42_hba/scalarphase+ap+TEC.parset'
        else:
            dummyparset = \
                '/home/weeren/scripts/rx42_hba/scalarphase+ap.parset'

    dummyparmdb = 'instrument_template'

    runbbs(
        mslist,
        skymodel,
        dummyparset,
        dummyparmdb,
        True,
        False,
        )
    if h:
        parmdb_a = 'instrument_amps2_smoothed'  # last/best ampplitude(+phase) parmdb
    else:
        parmdb_a = 'instrument_amps2'  # last/best ampplitude(+phase) parmdb
    parmdb_p = 'instrument_phase2'  # last/best CommonScalarPhase parmdb
    parmdbout = 'instrument_merged'

    for ms in mslist:
        create_merged_parmdb(
            ms,
            ms + '/' + pre_apply_parmdb,
            ms + '/' + parmdb_a,
            ms + '/' + parmdb_p,
            ms + '/' + dummyparmdb,
            ms + '/' + parmdbout,
            cellsizetime_a,
            cellsizetime_p,
            )

##################################################################
##################################################################
##################################################################
##################################################################

sys.exit()

### CALIBRATE BBS PHASE+AMP 4 ###

os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')
if FFT:
    os.system('casapy --nologger -c ft.py '
               + msinputlist + ' ' + imout + '.model' + ' '
              + str(nterms) + ' ' + str(wplanes))

# parmdb keep from previous step

skymodel = imout + '.skymodel'
if smooth:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset',
        parmdb + '_smoothed',
        True,
        False,
        )
else:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset',
        parmdb,
        True,
        False,
        )

if smooth:
    pre_apply_parmdb = parmdb + '_smoothed'  # need in merged if this is last calibration round
else:
    pre_apply_parmdb = parmdb  # need in merged if this is last calibration round

# solve +apply phases

if len(mslist) == 34:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
            '11,11,12', FFT)
if len(mslist) == 16:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC, '8,8',
            FFT)
else:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
            str(len(mslist)), FFT)
runbbs(
    mslist,
    skymodel,
    parset,
    'instrument_phase3',
    False,
    TEC,
    )

# solve amps

parmdb = 'instrument_amps3'
parset = create_amponly_parset(cellsizetime_a, FFT)
runbbs(
    mslist,
    skymodel,
    parset,
    parmdb,
    False,
    False,
    )

for ms in mslist:

  # remove outliers from the solutions

    if phasors:
        os.system('python smoothcal_rx42.py ' + ms + ' ' + ms + '/'
                  + parmdb + ' ' + ms + '/' + parmdb + '_smoothed')
    else:
        os.system('python /home/weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

if smooth:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb + '_smoothed',
        True,
        False,
        )
else:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb,
        True,
        False,
        )

### MAKE IMAGE 6 ###

imout, mask = make_image(datadir,cluster,'6',8,6,nterms,atrous_do,imsize,wplane)

### CREATE FINAL MODEL ###

skymodelf = 'im_cluster' + cluster + '.final.skymodel'
os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + skymodelf)
if FFT:
    os.system('casapy --nologger -c ft.py '
               + msinputlist + ' ' + imout + '.model' + ' '
              + str(nterms) + ' ' + str(wplanes))

### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
### INCLUDES SPLINE INTERPOLARION OF AMPS ###

if merge_parmdb:

    if phasors:
        dummyparset = \
            '/home/weeren/scripts/rx42_hba/scalarphase+amp.parset'
    else:
        if TEC:
            dummyparset = \
                '/home/weeren/scripts/rx42_hba/scalarphase+ap+TEC.parset'
        else:
            dummyparset = \
                '/home/weeren/scripts/rx42_hba/scalarphase+ap.parset'

    dummyparmdb = 'instrument_template'

    runbbs(
        mslist,
        skymodel,
        dummyparset,
        dummyparmdb,
        True,
        False,
        )
    if smooth:
        parmdb_a = 'instrument_amps3_smoothed'  # last/best ampplitude(+phase) parmdb
    else:
        parmdb_a = 'instrument_amps3'  # last/best ampplitude(+phase) parmdb
    parmdb_p = 'instrument_phase3'  # last/best CommonScalarPhase parmdb
    parmdbout = 'instrument_merged'

    for ms in mslist:
        create_merged_parmdb(
            ms,
            ms + '/' + pre_apply_parmdb,
            ms + '/' + parmdb_a,
            ms + '/' + parmdb_p,
            ms + '/' + dummyparmdb,
            ms + '/' + parmdbout,
            cellsizetime_a,
            cellsizetime_p,
            )

sys.exit()

### CALIBRATE BBS PHASE+AMP 5 ###

os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')
if FFT:
    os.system('casapy --nologger -c ft.py '
               + msinputlist + ' ' + imout + '.model' + ' '
              + str(nterms) + ' ' + str(wplanes))

# parmdb keep from previous step

skymodel = imout + '.skymodel'

# pre-apply amp solutions

if TEC:
    if smooth:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p_TEC.parset'
                ,
            parmdb + '_smoothed',
            True,
            False,
            )
    else:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p_TEC.parset'
                ,
            parmdb,
            True,
            False,
            )
else:
    if smooth:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset'
                ,
            parmdb + '_smoothed',
            True,
            False,
            )
    else:
        runbbs(
            mslist,
            skymodel,
            '/home/weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset'
                ,
            parmdb,
            True,
            False,
            )
if smooth:
    pre_apply_parmdb = parmdb + '_smoothed'  # need in merged if this is last calibration round
else:
    pre_apply_parmdb = parmdb  # need in merged if this is last calibration round

skymodel = imout + '.skymodel'

# solve +apply phases

if len(mslist) == 34:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
            '11,11,12', FFT)
if len(mslist) == 16:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC, '8,8',
            FFT)
else:
    parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
            str(len(mslist)), FFT)
runbbs(
    mslist,
    skymodel,
    parset,
    'instrument_phase4',
    False,
    TEC,
    )

# solve amps

parmdb = 'instrument_amps4'
parset = create_amponly_parset(cellsizetime_a, FFT)
runbbs(
    mslist,
    skymodel,
    parset,
    parmdb,
    False,
    False,
    )

for ms in mslist:

  # remove outliers from the solutions

    if phasors:
        os.system('python /home/weeren/scripts/rx42_hba/smoothcal_rx42.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')
    else:
        os.system('python /home/weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

if smooth:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb + '_smoothed',
        True,
        False,
        )
else:
    runbbs(
        mslist,
        skymodel,
        '/home/weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
        parmdb,
        True,
        False,
        )

### MAKE IMAGE 7 ###

imout, mask = make_image(datadir,cluster,'7',8,6,nterms,atrous_do,imsize,wplane)

### CREATE FINAL MODEL ###

skymodelf = 'im_cluster' + cluster + '.final.skymodel'
os.system('/home/weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
           + '-t ' + str(nterms) + ' ' + imout + '.model ' + skymodelf)

### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
### INCLUDES SPLINE INTERPOLARION OF AMPS ###

if merge_parmdb:

    if phasors:
        dummyparset = \
            '/home/weeren/scripts/rx42_hba/scalarphase+amp.parset'
    else:
        if TEC:
            dummyparset = \
                '/home/weeren/scripts/rx42_hba/scalarphase+ap+TEC.parset'
        else:
            dummyparset = \
                '/home/weeren/scripts/rx42_hba/scalarphase+ap.parset'

    dummyparmdb = 'instrument_template'

    runbbs(
        mslist,
        skymodel,
        dummyparset,
        dummyparmdb,
        True,
        False,
        )
    if smooth:
        parmdb_a = 'instrument_amps3_smoothed'  # last/best ampplitude(+phase) parmdb
    else:
        parmdb_a = 'instrument_amps3'  # last/best ampplitude(+phase) parmdb
    parmdb_p = 'instrument_phase3'  # last/best CommonScalarPhase parmdb
    parmdbout = 'instrument_merged'

    for ms in mslist:
        create_merged_parmdb(
            ms,
            ms + '/' + pre_apply_parmdb,
            ms + '/' + parmdb_a,
            ms + '/' + parmdb_p,
            ms + '/' + dummyparmdb,
            ms + '/' + parmdbout,
            cellsizetime_a,
            cellsizetime_p,
            )

