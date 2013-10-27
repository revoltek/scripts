#!/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('GTK')
import numpy
import os
import sys
import lofar.parmdb
import lofar.expion.parmdbmain
from scipy import interpolate
import time
from subprocess import Popen, PIPE
import pyrap.tables as pt
import uuid

# import psutil

pi = numpy.pi


# to run casa when not logged in
# Step 0. Run in screen
# Step 1. Xvfb :5 &
# Step 2. setenv DISPLAY :5
# run the program in screen and to exit screen: screen -d

def create_merged_parmdb_spline(
    parmdb_a,
    parmdb_p,
    parmdb_t,
    parmdb_out,
    cellsizetime_a,
    cellsizetime_b,
    ):

    pdb_a = lofar.parmdb.parmdb(parmdb_a)
    pdb_p = lofar.parmdb.parmdb(parmdb_p)
    pdb_t = lofar.parmdb.parmdb(parmdb_t)

    parms_a = pdb_a.getValuesGrid('*')
    parms_p = pdb_p.getValuesGrid('*')
    parms_t = pdb_t.getValuesGrid('*')

    os.system('rm -rf ' + parmdb_out)

    keynames_p = parms_p.keys()
    keynames_a = parms_a.keys()

    for key in keynames_p:

   # copy over the scalar phases

        scalarp = numpy.copy(parms_p[key]['values'][:, 0])
        parms_t[key]['values'][:, 0] = scalarp

    for key in keynames_a:
        print key
        tmp1 = numpy.copy(parms_t[key]['values'][:, 0])
        tmp2 = numpy.copy(parms_a[key]['values'][:, 0])

        time_t = 1.0 * numpy.arange(len(tmp1))
        time_a = numpy.float(cellsizetime_a / cellsizetime_p) \
            * numpy.arange(len(tmp2)) + 0.5 * cellsizetime_a \
            / cellsizetime_p

        el = (numpy.max(time_t) - numpy.max(time_a[:-1])) / 2 \
            + numpy.max(time_a[:-1])
        time_a[len(time_a) - 1] = el

   # to avoid edge effects in spline add values edges to be the last values from BBS

        time_a = numpy.copy(numpy.append(0, time_a))
        time_a = numpy.append(time_a, numpy.max(time_t))

        val_start = tmp2[0]
        val_end = tmp2[len(tmp2) - 1]

        tmp2 = numpy.copy(numpy.append(val_start, tmp2))
        tmp2 = numpy.append(tmp2, val_end)

        x = numpy.copy(time_a)
        y = numpy.copy(tmp2)

   # SPLINE INTERPOL of AMPS

        tck = interpolate.splrep(x, y, s=0)
        xnew = time_t
        ynew = interpolate.splev(xnew, tck, der=0)

        parms_t[key]['values'][:, 0] = ynew

    lofar.expion.parmdbmain.store_parms(parmdb_out, parms_t,
            create_new=True)
    return


def create_merged_parmdb(
    ms,
    pre_apply_parmdb,
    parmdb_a,
    parmdb_p,
    parmdb_t,
    parmdb_out,
    cellsizetime_a,
    cellsizetime_b,
    ):

    pdb_pre = lofar.parmdb.parmdb(pre_apply_parmdb)
    pdb_a = lofar.parmdb.parmdb(parmdb_a)
    pdb_p = lofar.parmdb.parmdb(parmdb_p)
    pdb_t = lofar.parmdb.parmdb(parmdb_t)

    parms_pre = pdb_pre.getValuesGrid('*')
    parms_a = pdb_a.getValuesGrid('*')
    parms_p = pdb_p.getValuesGrid('*')
    parms_t = pdb_t.getValuesGrid('*')

    os.system('rm -rf ' + parmdb_out)

    keynames_p = parms_p.keys()
    keynames_a = parms_a.keys()

    length = numpy.int(cellsizetime_b)
    for key in keynames_p:

   # copy over the CommonScalar phases

   # just to get axis length

        tmp1 = numpy.copy(parms_t[key]['values'][:, 0])
        for idx in range(len(tmp1)):
            el = numpy.float(idx) / numpy.float(length)
            el = numpy.int(numpy.floor(el))

            parms_t[key]['values'][idx, 0] = \
                numpy.copy(parms_p[key]['values'][el, 0])

    pol_list = ['0:0', '1:1']
    gain = 'Gain'
    anttab = pt.table(ms + '/ANTENNA')
    antenna_list = anttab.getcol('NAME')
    anttab.close()

 # length = numpy.int(cellsizetime_a/cellsizetime_b)

    length = numpy.int(cellsizetime_a)

    for pol in pol_list:
        for antenna in antenna_list:
            real1 = numpy.copy(parms_pre[gain + ':' + pol + ':Real:'
                               + antenna]['values'][:, 0])
            real2 = numpy.copy(parms_a[gain + ':' + pol + ':Real:'
                               + antenna]['values'][:, 0])
            imag1 = numpy.copy(parms_pre[gain + ':' + pol + ':Imag:'
                               + antenna]['values'][:, 0])
            imag2 = numpy.copy(parms_a[gain + ':' + pol + ':Imag:'
                               + antenna]['values'][:, 0])

     # just to get axis length

            tmp1 = numpy.copy(parms_t[gain + ':' + pol + ':Imag:'
                              + antenna]['values'][:, 0])

            G1 = real1 + 1j * imag1
            G2 = real2 + 1j * imag2
            Gnew = G1 * G2

            for idx in range(len(tmp1)):
                el = numpy.float(idx) / numpy.float(length)
                el = numpy.int(numpy.floor(el))

                parms_t[gain + ':' + pol + ':Imag:' + antenna]['values'
                        ][idx, 0] = numpy.copy(numpy.imag(Gnew[el]))
                parms_t[gain + ':' + pol + ':Real:' + antenna]['values'
                        ][idx, 0] = numpy.copy(numpy.real(Gnew[el]))

 # for key in keynames_a:
 #  print key
 #  tmp1 = numpy.copy(parms_t[key]['values'][:,0])
 #  tmp2 = numpy.copy(parms_a[key]['values'][:,0])
 #  length = numpy.int(cellsizetime_a/cellsizetime_b)

 #  for idx in range(len(tmp1)):
 #
 #    el =  numpy.float(idx)/numpy.float(length)
 #    el =  numpy.int(numpy.floor(el))
 #    #print idx, el
 #    parms_t[key]['values'][idx,0] = tmp2[el]

    lofar.expion.parmdbmain.store_parms(parmdb_out, parms_t,
            create_new=True)
    return


def make_image(
    mslist,
    cluster,
    callnumber,
    threshpix,
    threshisl,
    nterms,
    atrous_do,
    imsize,
    ):

    do_mask = True
    niter = 1000  # 7500 causes nasty clean artifacts
    mscale = 'False'

    average = True  # average the data a lot to speed up the imaging,

 # ONLY average for small FOV, otherwise timesmearing is a problem
 # average data

    if average:
        for ms in mslist:
            ndppp_parset = ms + '_NDPPP.parset'
            os.system('rm -f ' + ndppp_parset)
            output = ms + '.tmpavg'
            f = open(ndppp_parset, 'w')
            f.write('msin = %s\n' % ms)
            f.write('msin.datacolumn = CORRECTED_DATA\n')
            f.write('msout = %s\n' % output)
            f.write('steps=[avg]\n')
            f.write('rficonsole.type=aoflagger\n')
            f.write('avg.type = squash\n')
            f.write('avg.freqstep = 1\n')
            f.write('avg.timestep = 12\n')  # is the default
            f.close()
            cmd = 'ps -u weeren | grep NDPPP | wc -l'
            output = numpy.int(Popen(cmd, shell=True,
                               stdout=PIPE).communicate()[0])

            while output > 2:
                time.sleep(10)
                output = numpy.int(Popen(cmd, shell=True,
                                   stdout=PIPE).communicate()[0])

                pid = Popen('pidof NDPPP', shell=True,
                            stdout=PIPE).communicate()[0]

      # print pid

                pid_list = pid.split(' ')
                for pidnr in pid_list:

        # print 'PID', pidnr

                    os.system('kill -s CONT ' + str(pidnr))

            os.system('NDPPP ' + ndppp_parset + '&')

    output = numpy.int(Popen(cmd, shell=True,
                       stdout=PIPE).communicate()[0])
    while output > 0:
        time.sleep(10)
        output = numpy.int(Popen(cmd, shell=True,
                           stdout=PIPE).communicate()[0])
        pid = Popen('pidof NDPPP', shell=True,
                    stdout=PIPE).communicate()[0]

      # print pid

        pid_list = pid.split(' ')
        for pidnr in pid_list:

        # print 'PID', numpy.int(pidnr)

            os.system('kill -s CONT ' + str(pidnr))

    ms = ''
    for m in mslist:
        if average:
            ms = ms + ' ' + m + '.tmpavg'
        else:
            ms = ms + ' ' + m

    imout = 'im' + callnumber + '_cluster' + cluster + 'nm'
    print ms + ' ' + imout + ' ' + 'None' + ' ' + '1mJy' + ' ' \
        + str(niter) + ' ' + str(nterms)

    if do_mask:
        os.system('casapy --nologger -c ~weeren/scripts/rx42_hba/casapy_cleanv3.py '
                   + ms + ' ' + imout + ' ' + 'None' + ' ' + '1mJy'
                  + ' ' + str(niter) + ' ' + str(nterms) + ' '
                  + str(imsize) + ' ' + mscale)

   # make mask

        if nterms > 1:
            os.system('python ~weeren/scripts/rx42_hba/makecleanmask.py --threshpix '
                       + str(threshpix) + ' --threshisl '
                      + str(threshisl) + ' --atrous_do '
                      + str(atrous_do) + ' ' + imout + '.image.tt0')
        else:
            os.system('python ~weeren/scripts/rx42_hba/makecleanmask.py --threshpix '
                       + str(threshpix) + ' --threshisl '
                      + str(threshisl) + ' --atrous_do '
                      + str(atrous_do) + ' ' + imout + '.image')

   # clean image with manual mask
        mask = imout + '.cleanmask'

    imout = 'im' + callnumber + '_cluster' + cluster

    os.system('casapy --nologger -c ~weeren/scripts/rx42_hba/casapy_cleanv3.py '
               + ms + ' ' + imout + ' ' + mask + ' ' + '1mJy' + ' '
              + str(niter) + ' ' + str(nterms) + ' ' + str(imsize) + ' '
               + mscale)

 # convert to FITS
    if nterms > 1:
        os.system('image2fits in=' + imout + '.image.tt0' + ' ' + 'out='
                   + imout + '.fits')
    else:
        os.system('image2fits in=' + imout + '.image' + ' ' + 'out='
                  + imout + '.fits')

    if average:
        print 'rm -rf ' + ms
        os.system('rm -rf ' + ms)

    return (imout, mask)


def runbbs(
    mslist,
    skymodel,
    parset,
    parmdb,
    applycal,
    TEC,
    ):

 # NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)
 # NOTE WORK FROM DATA  if TEC

    if TEC:
        if len(mslist) == 10000:  # 34 does not work now!!!!!!

   # this is a very special case for a full run, manually here to run with 3 solver controls

            vdslist = ''
            os.system('makevds /home/weeren/cep2.clusterdesc ' + ms)
            vdslist = vdslist + ms + '.vds '
        else:

            vdslist = ''
            for ms in mslist:
                os.system('makevds /home/weeren/cep2.clusterdesc ' + ms)
                vdslist = vdslist + ms + '.vds '
            gds = 'tec.gds'
            print vdslist
            os.system('combinevds ' + gds + ' ' + vdslist)

            key = str(uuid.uuid4())
            cmd1 = 'calibrate -f --key ' + key \
                + ' --cluster-desc /home/weeren/cep2.clusterdesc --instrument-name ' \
                + parmdb + ' '
            cmd2 = '--db ldb001 --db-user postgres ' + gds + ' ' \
                + parset + ' ' + skymodel + ' ' + '. > ' + gds \
                + '.bbslog'
            cmd = cmd1 + cmd2
            print cmd
            os.system(cmd)
    else:

        for ms in mslist:
            log = ms + '.bbslog'
            if applycal:
                cmd = 'calibrate-stand-alone --parmdb-name ' + parmdb \
                    + ' ' + ms + ' ' + parset + ' ' + skymodel + '>' \
                    + log + '&'
            else:
                cmd = 'calibrate-stand-alone -f --parmdb-name ' \
                    + parmdb + ' ' + ms + ' ' + parset + ' ' + skymodel \
                    + '>' + log + '&'
            print cmd
            os.system(cmd)
        time.sleep(10)

        done = 0
        while done < len(mslist):
            done = 0
            for ms in mslist:
                cmd = \
                    "grep 'INFO - bbs-reducer terminated successfully.' " \
                    + ms + '.bbslog'
                output = Popen(cmd, shell=True,
                               stdout=PIPE).communicate()[0]
                if output[0:4] == 'INFO':
                    done = done + 1
                    print ms, 'is done'
            time.sleep(5)

    return


def create_scalarphase_parset(timestep, TEC, groups):
    bbs_parset = 'scalarphase.parset'
    os.system('rm -f ' + bbs_parset)
    f = open(bbs_parset, 'w')

    if TEC:
        f.write('Strategy.InputColumn = DATA\n')
        f.write('Strategy.UseSolver   = T\n')
    else:

        f.write('Strategy.InputColumn = MODEL_DATA\n')
        f.write('Strategy.UseSolver   = F\n')

    f.write('Strategy.ChunkSize   = 500\n')
    f.write('Strategy.Steps       = [solve,correct]\n')

    f.write('Step.solve.Model.Sources                = []\n')
    f.write('Step.solve.Model.Cache.Enable           = T\n')
    f.write('Step.solve.Model.Phasors.Enable         = F\n')
    f.write('Step.solve.Model.DirectionalGain.Enable = F\n')
    f.write('Step.solve.Model.Gain.Enable            = F\n')

    if TEC:
        f.write('Step.solve.Model.TEC.Enable            = T\n')
        f.write('Step.solve.Solve.CalibrationGroups     = [%s]\n'
                % groups)
        f.write('Step.solve.Model.Clock.Enable          = T\n')

    f.write('Step.solve.Model.Rotation.Enable        = F\n')
    f.write('Step.solve.Model.CommonScalarPhase.Enable = T \n')
    f.write('Step.solve.Operation                    = SOLVE\n')

    if TEC:

    # f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*"]\n')

        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*","Clock:*"]\n'
                )
    else:
        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*"]\n'
                )

    f.write('Step.solve.Solve.CellSize.Freq          = 0\n')
    f.write('Step.solve.Solve.CellSize.Time          = %s\n'
            % str(timestep))
    f.write('Step.solve.Solve.CellChunkSize          = 500\n')
    f.write('Step.solve.Solve.PropagateSolutions     = T\n')
    f.write('Step.solve.Solve.Options.MaxIter        = 100\n')
    f.write('Step.solve.Solve.Options.LMFactor       = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs    = F\n')
    f.write('Step.solve.Solve.Options.UseSVD         = T\n')
    f.write('Step.solve.Model.Beam.Enable            = F\n')
    f.write('Step.solve.Solve.UVRange		   = [80]\n')
    f.write('Step.solve.Solve.Mode                   = COMPLEX \n')

    f.write('Step.correct.Model.Sources                 = []\n')
    f.write('Step.correct.Model.CommonScalarPhase.Enable= T\n')
    f.write('Step.correct.Model.Cache.Enable            = T\n')

    if TEC:
        f.write('Step.correct.Model.TEC.Enable  = T\n')
        f.write('Step.correct.Model.Clock.Enable  = T\n')

    f.write('Step.correct.Model.DirectionalGain.Enable  = F\n')
    f.write('Step.correct.Model.Gain.Enable             = F\n')
    f.write('Step.correct.Model.Phasors.Enable          = F\n')
    f.write('Step.correct.Operation                     = CORRECT\n')
    f.write('Step.correct.Output.Column                 = CORRECTED_DATA\n'
            )
    f.write('Step.correct.Model.Beam.Enable             = F\n')
    f.write('Step.correct.Output.WriteCovariance        = F\n')

    f.close()
    return bbs_parset


def create_scalarphase_parset_p(timestep, TEC, groups):
    bbs_parset = 'scalarphase_p.parset'
    os.system('rm -f ' + bbs_parset)
    f = open(bbs_parset, 'w')
    if TEC:
        f.write('Strategy.InputColumn = DATA\n')
        f.write('Strategy.UseSolver   = T\n')
    else:
        f.write('Strategy.InputColumn = MODEL_DATA\n')
        f.write('Strategy.UseSolver   = F\n')

    f.write('Strategy.ChunkSize   = 500\n')
    f.write('Strategy.Steps       = [solve,correct]\n')

    f.write('Step.solve.Model.Sources                = []\n')
    f.write('Step.solve.Model.Cache.Enable           = T\n')
    f.write('Step.solve.Model.Phasors.Enable         = F\n')
    f.write('Step.solve.Model.DirectionalGain.Enable = F\n')
    f.write('Step.solve.Model.Gain.Enable            = F\n')

    if TEC:
        f.write('Step.solve.Model.TEC.Enable            = T\n')
        f.write('Step.solve.Model.Clock.Enable            = T\n')
        f.write('Step.solve.Solve.CalibrationGroups     = [%s]\n'
                % groups)

    f.write('Step.solve.Model.Rotation.Enable        = F\n')
    f.write('Step.solve.Model.CommonScalarPhase.Enable = T \n')
    f.write('Step.solve.Operation                    = SOLVE\n')

    if TEC:

    # f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*"]\n')

        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*","Clock:*"]\n'
                )
    else:
        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*"]\n'
                )

    f.write('Step.solve.Solve.CellSize.Freq          = 0\n')
    f.write('Step.solve.Solve.CellSize.Time          = %s\n'
            % str(timestep))
    f.write('Step.solve.Solve.CellChunkSize          = 500\n')
    f.write('Step.solve.Solve.PropagateSolutions     = T\n')
    f.write('Step.solve.Solve.Options.MaxIter        = 100\n')
    f.write('Step.solve.Solve.Options.LMFactor       = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs    = F\n')
    f.write('Step.solve.Solve.Options.UseSVD         = T\n')
    f.write('Step.solve.Model.Beam.Enable            = F\n')
    f.write('Step.solve.Solve.UVRange		   = [80]\n')
    f.write('Step.solve.Solve.Mode                   = COMPLEX \n')

    f.write('Step.correct.Model.Sources                 = []\n')
    f.write('Step.correct.Model.CommonScalarPhase.Enable= T\n')
    f.write('Step.correct.Model.Cache.Enable            = T\n')
    f.write('Step.correct.Model.DirectionalGain.Enable  = F\n')

    if TEC:
        f.write('Step.correct.Model.TEC.Enable  = T\n')
        f.write('Step.correct.Model.Clock.Enable  = T\n')

    f.write('Step.correct.Model.Gain.Enable             = F\n')
    f.write('Step.correct.Model.Phasors.Enable          = F\n')
    f.write('Step.correct.Operation                     = CORRECT\n')
    f.write('Step.correct.Output.Column                 = CORRECTED_DATA_PHASE\n'
            )
    f.write('Step.correct.Model.Beam.Enable             = F\n')
    f.write('Step.correct.Output.WriteCovariance        = F\n')

    f.close()
    return bbs_parset


def create_amponly_parset(timestep):
    bbs_parset = 'amplitudeonly.parset'
    os.system('rm -f ' + bbs_parset)
    f = open(bbs_parset, 'w')

    f.write('Strategy.InputColumn = CORRECTED_DATA_PHASE\n')
    f.write('Strategy.ChunkSize   = 480\n')
    f.write('Strategy.UseSolver   = F\n')
    f.write('Strategy.Steps       = [solve]\n')

    f.write('Step.solve.Model.Sources                  = []\n')
    f.write('Step.solve.Model.Cache.Enable             = T\n')
    f.write('Step.solve.Model.Phasors.Enable           = F\n')  # # T
    f.write('Step.solve.Model.DirectionalGain.Enable   = F\n')
    f.write('Step.solve.Model.Gain.Enable              = T\n')
    f.write('Step.solve.Model.Rotation.Enable          = F\n')
    f.write('Step.solve.Model.CommonScalarPhase.Enable = F\n')
    f.write('Step.solve.Operation                      = SOLVE\n')
    f.write('Step.solve.Solve.Parms                    = ["Gain:0:0:*","Gain:1:1:*"]\n'
            )  # #  ["Gain:0:0:Ampl:*","Gain:1:1:Ampl:*"]
    f.write('Step.solve.Solve.CellSize.Freq            = 0\n')
    f.write('Step.solve.Solve.CellSize.Time            = %s\n'
            % str(timestep))
    f.write('Step.solve.Solve.CellChunkSize            = 480\n')
    f.write('Step.solve.Solve.PropagateSolutions       = T\n')
    f.write('Step.solve.Solve.Options.MaxIter          = 100\n')
    f.write('Step.solve.Solve.Options.LMFactor         = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs      = F\n')
    f.write('Step.solve.Solve.Options.UseSVD           = T\n')
    f.write('Step.solve.Model.Beam.Enable              = F\n')
    f.write('Step.solve.Solve.UVRange		     = [80]\n')  # # [600]
    f.write('Step.solve.Solve.Mode                     = COMPLEX\n')

    f.close()
    return bbs_parset


def create_scalarphase_parset_p2(timestep, TEC, groups):
    bbs_parset = 'scalarphase_p2.parset'
    os.system('rm -f ' + bbs_parset)
    f = open(bbs_parset, 'w')
    f.write('Strategy.InputColumn = CORRECTED_DATA_AMP\n')
    f.write('Strategy.ChunkSize   = 500\n')

    if TEC:
        f.write('Strategy.UseSolver   = T\n')
    else:
        f.write('Strategy.UseSolver   = F\n')

    f.write('Strategy.Steps       = [solve,correct]\n')

    f.write('Step.solve.Model.Sources                = []\n')
    f.write('Step.solve.Model.Cache.Enable           = T\n')
    f.write('Step.solve.Model.Phasors.Enable         = F\n')
    f.write('Step.solve.Model.DirectionalGain.Enable = F\n')
    f.write('Step.solve.Model.Gain.Enable            = F\n')
    if TEC:
        f.write('Step.solve.Model.TEC.Enable            = T\n')
        f.write('Step.solve.Model.Clock.Enable          = T\n')
        f.write('Step.solve.Solve.CalibrationGroups     = [%s]\n'
                % groups)

    f.write('Step.solve.Model.Rotation.Enable        = F\n')
    f.write('Step.solve.Model.CommonScalarPhase.Enable = T \n')
    f.write('Step.solve.Operation                    = SOLVE\n')

    if TEC:
#        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*"]\n'
                )
        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*","TEC:*","Clock:*"]\n'
                )
    else:
        f.write('Step.solve.Solve.Parms                  = ["CommonScalarPhase:*"]\n'
                )

    f.write('Step.solve.Solve.CellSize.Freq          = 0\n')
    f.write('Step.solve.Solve.CellSize.Time          = %s\n'
            % str(timestep))
    f.write('Step.solve.Solve.CellChunkSize          = 500\n')
    f.write('Step.solve.Solve.PropagateSolutions     = T\n')
    f.write('Step.solve.Solve.Options.MaxIter        = 100\n')
    f.write('Step.solve.Solve.Options.LMFactor       = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs    = F\n')
    f.write('Step.solve.Solve.Options.UseSVD         = T\n')
    f.write('Step.solve.Model.Beam.Enable            = F\n')
    f.write('Step.solve.Solve.UVRange		   = [80]\n')
    f.write('Step.solve.Solve.Mode                   = COMPLEX \n')

    f.write('Step.correct.Model.Sources                 = []\n')
    f.write('Step.correct.Model.CommonScalarPhase.Enable= T\n')
    f.write('Step.correct.Model.Cache.Enable            = T\n')
    f.write('Step.correct.Model.DirectionalGain.Enable  = F\n')
    f.write('Step.correct.Model.Gain.Enable             = F\n')

    if TEC:
        f.write('Step.correct.Model.TEC.Enable  = T\n')
        f.write('Step.correct.Model.Clock.Enable  = T\n')

    f.write('Step.correct.Model.Phasors.Enable          = F\n')
    f.write('Step.correct.Operation                     = CORRECT\n')
    f.write('Step.correct.Output.Column                 = CORRECTED_DATA_PHASE\n'
            )
    f.write('Step.correct.Model.Beam.Enable             = F\n')
    f.write('Step.correct.Output.WriteCovariance        = F\n')

    f.close()
    return bbs_parset


el = len(sys.argv)

mslist = sys.argv[1:el - 7]
cluster = str(sys.argv[el - 7])
atrous_do = str(sys.argv[el - 6])
imsize = numpy.int(sys.argv[el - 5])

nterms = numpy.int(sys.argv[el - 4])  # only 1 to 3 is supported !!
cellsizetime_a = numpy.int(sys.argv[el - 3])
cellsizetime_p = numpy.int(sys.argv[el - 2])
TECi = str(sys.argv[el - 1])

TEC = False
if TECi == 'True':
    TEC = True

print 'mslist', mslist
print 'source', cluster
print 'atrous_do', atrous_do
print 'imsize', imsize

# sys.exit()

merge_parmdb = True
phasors = False  # if true only solve for amps on long timescales

# cellsizetime_a          = 120
# cellsizetime_p          = 1
# nterms                  = 1  # only 1 to 3 is supported !!

# if merge_parmdb:
#   pre_apply_parmdb = 'instrument_amps0' + '_smoothed' # need in merged if this is last calibration round
#   if phasors:
#      dummyparset = '~weeren/scripts/rx42_hba/scalarphase+amp.parset'
#   else:
#      dummyparset = '~weeren/scripts/rx42_hba/scalarphase+ap.parset'
#   dummyparmdb = 'instrument_template'
#
#   #runbbs(mslist, skymodel,dummyparset, dummyparmdb, True, False)
#
#   parmdb_a    = 'instrument_amps1_smoothed'  # last/best ampplitude(+phase) parmdb
#   parmdb_p    = 'instrument_phase1'          # last/best CommonScalarPhase parmdb
#   parmdbout   = 'instrument_merged'
#
#   for ms in mslist:
#     create_merged_parmdb(ms, ms+'/'+pre_apply_parmdb, ms+'/'+parmdb_a, ms+'/'+parmdb_p, ms+'/'+dummyparmdb,ms+'/'+parmdbout,cellsizetime_a,cellsizetime_p)
#
# sys.exit()

# runbbs(mslist, 'bla.skymodel', 'bla.parset', 'parmdbbla', False, True)
# sys.exit()

#### MAKE IMAGE 0 ###

(imout, mask) = make_image(
    mslist,
    cluster,
    '0',
    10,
    6,
    nterms,
    atrous_do,
    imsize,
    )

#####################

### CALIBRATE WITH BBS PHASE ONLY 1 ###
# create skymodel for BBS

os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')

# phase only calibrate

skymodel = imout + '.skymodel'

# if len(mslist) == 34:
#  parset   = create_scalarphase_parset(cellsizetime_p, TEC, '11,11,12')
# if len(mslist) == 16:
#  parset   = create_scalarphase_parset(cellsizetime_p, TEC, '8,8')
# else:

parset = create_scalarphase_parset(cellsizetime_p, TEC,
                                   str(len(mslist)))
runbbs(
    mslist,
    skymodel,
    parset,
    'instrument',
    False,
    TEC,
    )

# NOTE WORK FROM MODEL_DATA (contains correct phase data from 10SB calibration)
######################################

### MAKE IMAGE 1 ###

(imout, mask) = make_image(
    mslist,
    cluster,
    '1',
    10,
    10,
    nterms,
    atrous_do,
    imsize,
    )

####################

### CALIBRATE WITH BBS PHASE ONLY 2 ###
# create skymodel for BBS

os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')

# phase only calibrate

skymodel = imout + '.skymodel'

# if len(mslist) == 34:
#  parset   = create_scalarphase_parset(cellsizetime_p, TEC, '11,11,12')
# if len(mslist) == 16:
#  parset   = create_scalarphase_parset(cellsizetime_p, TEC, '8,8')
# else:

parset = create_scalarphase_parset(cellsizetime_p, TEC,
                                   str(len(mslist)))
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

(imout, mask) = make_image(
    mslist,
    cluster,
    '2',
    10,
    10,
    nterms,
    atrous_do,
    imsize,
    )

####################

### CALIBRATE WITH BBS PHASE+AMP 1 ###

os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')

skymodel = imout + '.skymodel'

# if len(mslist) == 34:
#  parset   = create_scalarphase_parset_p(cellsizetime_p, TEC, '11,11,12')
# if len(mslist) == 16:
#  parset   = create_scalarphase_parset_p(cellsizetime_p, TEC, '8,8')
# else:

parset = create_scalarphase_parset_p(cellsizetime_p, TEC,
        str(len(mslist)))

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
parset = create_amponly_parset(cellsizetime_a)
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
        os.system('python ~weeren/scripts/rx42_hba/smoothcal_rx42.py '
                  + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')
    else:
        os.system('python ~weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

runbbs(
    mslist,
    skymodel,
    '~weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
    parmdb + '_smoothed',
    True,
    False,
    )

### MAKE IMAGE 3 ###

(imout, mask) = make_image(
    mslist,
    cluster,
    '3',
    10,
    10,
    nterms,
    atrous_do,
    imsize,
    )

#### CALIBRATE  BBS PHASE+AMP 2 ###
# make model

os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')

# parmdb keep from previous step

skymodel = imout + '.skymodel'

# pre-apply amp solutions

if TEC:
    runbbs(
        mslist,
        skymodel,
        '~weeren/scripts/rx42_hba/apply_amplitudeonly_p_TEC.parset',
        parmdb + '_smoothed',
        True,
        False,
        )
else:
    runbbs(
        mslist,
        skymodel,
        '~weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset',
        parmdb + '_smoothed',
        True,
        False,
        )

pre_apply_parmdb = parmdb + '_smoothed'  # need in merged if this is last calibration round

# phase only cal

skymodel = imout + '.skymodel'

# if len(mslist) == 34:
#  parset   = create_scalarphase_parset_p2(cellsizetime_p, TEC, '11,11,12')
# if len(mslist) == 16:
#  parset   = create_scalarphase_parset_p2(cellsizetime_p, TEC, '8,8')
# else:

parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
        str(len(mslist)))
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
parset = create_amponly_parset(cellsizetime_a)
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
        os.system('python ~weeren/scripts/rx42_hba/smoothcal_rx42.py '
                  + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')
    else:
        os.system('python ~weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

runbbs(
    mslist,
    skymodel,
    '~weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
    parmdb + '_smoothed',
    True,
    False,
    )

### MAKE IMAGE 4 ###

(imout, mask) = make_image(
    mslist,
    cluster,
    '4',
    10,
    10,
    nterms,
    atrous_do,
    imsize,
    )

### CREATE FINAL MODEL ###

skymodelf = 'im_cluster' + cluster + '.final.skymodel'
os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + skymodelf)

### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
### INCLUDES SPLINE INTERPOLARION OF AMPS ###

if merge_parmdb:

    if phasors:
        dummyparset = '~weeren/scripts/rx42_hba/scalarphase+amp.parset'
    else:
        if TEC:
            dummyparset = \
                '~weeren/scripts/rx42_hba/scalarphase+ap+TEC+clock.parset'
        else:
            dummyparset = \
                '~weeren/scripts/rx42_hba/scalarphase+ap.parset'

    dummyparmdb = 'instrument_template'

    runbbs(
        mslist,
        skymodel,
        dummyparset,
        dummyparmdb,
        True,
        False,
        )

    parmdb_a = 'instrument_amps1_smoothed'  # last/best ampplitude(+phase) parmdb
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

os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')

# parmdb keep from previous step

skymodel = imout + '.skymodel'

# pre-apply amp solutions

if TEC:
    runbbs(
        mslist,
        skymodel,
        '~weeren/scripts/rx42_hba/apply_amplitudeonly_p_TEC.parset',
        parmdb + '_smoothed',
        True,
        False,
        )
else:
    runbbs(
        mslist,
        skymodel,
        '~weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset',
        parmdb + '_smoothed',
        True,
        False,
        )

pre_apply_parmdb = parmdb + '_smoothed'  # need in merged if this is last calibration round

# solve +apply phases
# if len(mslist) == 34:
#  parset   = create_scalarphase_parset_p2(cellsizetime_p, TEC, '11,11,12')
# if len(mslist) == 16:
#  parset   = create_scalarphase_parset_p2(cellsizetime_p, TEC, '8,8')
# else:

parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
        str(len(mslist)))
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
parset = create_amponly_parset(cellsizetime_a)
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
        os.system('python ~weeren/scripts/rx42_hba/smoothcal_rx42.py '
                  + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')
    else:
        os.system('python ~weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

runbbs(
    mslist,
    skymodel,
    '~weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
    parmdb + '_smoothed',
    True,
    False,
    )

### MAKE IMAGE 5 ###

(imout, mask) = make_image(
    mslist,
    cluster,
    '5',
    8,
    8,
    nterms,
    atrous_do,
    imsize,
    )

### CREATE FINAL MODEL ###

skymodelf = 'im_cluster' + cluster + '.final.skymodel'
os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + skymodelf)

### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
### INCLUDES SPLINE INTERPOLARION OF AMPS ###

if merge_parmdb:

    if phasors:
        dummyparset = '~weeren/scripts/rx42_hba/scalarphase+amp.parset'
    else:
        if TEC:
            dummyparset = \
                '~weeren/scripts/rx42_hba/scalarphase+ap+TEC+clock.parset'
        else:
            dummyparset = \
                '~weeren/scripts/rx42_hba/scalarphase+ap.parset'

    dummyparmdb = 'instrument_template'

    runbbs(
        mslist,
        skymodel,
        dummyparset,
        dummyparmdb,
        True,
        False,
        )

    parmdb_a = 'instrument_amps2_smoothed'  # last/best ampplitude(+phase) parmdb
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

os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')

# parmdb keep from previous step

skymodel = imout + '.skymodel'
runbbs(
    mslist,
    skymodel,
    '~weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset',
    parmdb + '_smoothed',
    True,
    False,
    )
pre_apply_parmdb = parmdb + '_smoothed'  # need in merged if this is last calibration round

# solve +apply phases
# if len(mslist) == 34:
#  parset   = create_scalarphase_parset_p2(cellsizetime_p, TEC, '11,11,12')
# if len(mslist) == 16:
#  parset   = create_scalarphase_parset_p2(cellsizetime_p, TEC, '8,8')
# else:

parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
        str(len(mslist)))
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
parset = create_amponly_parset(cellsizetime_a)
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
        os.system('python ~weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

runbbs(
    mslist,
    skymodel,
    '~weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
    parmdb + '_smoothed',
    True,
    False,
    )

### MAKE IMAGE 6 ###

(imout, mask) = make_image(
    mslist,
    cluster,
    '6',
    8,
    6,
    nterms,
    atrous_do,
    imsize,
    )

### CREATE FINAL MODEL ###

skymodelf = 'im_cluster' + cluster + '.final.skymodel'
os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + skymodelf)

### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
### INCLUDES SPLINE INTERPOLARION OF AMPS ###

if merge_parmdb:

    if phasors:
        dummyparset = '~weeren/scripts/rx42_hba/scalarphase+amp.parset'
    else:
        if TEC:
            dummyparset = \
                '~weeren/scripts/rx42_hba/scalarphase+ap+TEC.parset'
        else:
            dummyparset = \
                '~weeren/scripts/rx42_hba/scalarphase+ap.parset'

    dummyparmdb = 'instrument_template'

    runbbs(
        mslist,
        skymodel,
        dummyparset,
        dummyparmdb,
        True,
        False,
        )

    parmdb_a = 'instrument_amps3_smoothed'  # last/best ampplitude(+phase) parmdb
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

os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + imout
          + '.skymodel')

# parmdb keep from previous step

skymodel = imout + '.skymodel'

# pre-apply amp solutions

if TEC:
    runbbs(
        mslist,
        skymodel,
        '~weeren/scripts/rx42_hba/apply_amplitudeonly_p_TEC.parset',
        parmdb + '_smoothed',
        True,
        False,
        )
else:
    runbbs(
        mslist,
        skymodel,
        '~weeren/scripts/rx42_hba/apply_amplitudeonly_p.parset',
        parmdb + '_smoothed',
        True,
        False,
        )

pre_apply_parmdb = parmdb + '_smoothed'  # need in merged if this is last calibration round

skymodel = imout + '.skymodel'

# solve +apply phases

# if len(mslist) == 34:
#  parset   = create_scalarphase_parset_p2(cellsizetime_p, TEC, '11,11,12')
# if len(mslist) == 16:
#  parset   = create_scalarphase_parset_p2(cellsizetime_p, TEC, '8,8')
# else:

parset = create_scalarphase_parset_p2(cellsizetime_p, TEC,
        str(len(mslist)))
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
parset = create_amponly_parset(cellsizetime_a)
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
        os.system('python ~weeren/scripts/rx42_hba/smoothcal_rx42.py '
                  + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')
    else:
        os.system('python ~weeren/scripts/rx42_hba/smoothcal_rx42_nophasors.py '
                   + ms + ' ' + ms + '/' + parmdb + ' ' + ms + '/'
                  + parmdb + '_smoothed')

# apply amps

runbbs(
    mslist,
    skymodel,
    '~weeren/scripts/rx42_hba/apply_amplitudeonly.parset',
    parmdb + '_smoothed',
    True,
    False,
    )

### MAKE IMAGE 7 ###

(imout, mask) = make_image(
    mslist,
    cluster,
    '7',
    8,
    6,
    nterms,
    atrous_do,
    imsize,
    )

### CREATE FINAL MODEL ###

skymodelf = 'im_cluster' + cluster + '.final.skymodel'
os.system('~weeren/scripts/rx42_hba/casapy2bbs.py -m ' + mask + ' '
          + '-t ' + str(nterms) + ' ' + imout + '.model ' + skymodelf)

### CREATED MERGED PARMDB SCALARPHASE+AMPS ###
### INCLUDES SPLINE INTERPOLARION OF AMPS ###

if merge_parmdb:

    if phasors:
        dummyparset = '~weeren/scripts/rx42_hba/scalarphase+amp.parset'
    else:
        if TEC:
            dummyparset = \
                '~weeren/scripts/rx42_hba/scalarphase+ap+TEC.parset'
        else:
            dummyparset = \
                '~weeren/scripts/rx42_hba/scalarphase+ap.parset'

    dummyparmdb = 'instrument_template'

    runbbs(
        mslist,
        skymodel,
        dummyparset,
        dummyparmdb,
        True,
        False,
        )

    parmdb_a = 'instrument_amps3_smoothed'  # last/best ampplitude(+phase) parmdb
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

