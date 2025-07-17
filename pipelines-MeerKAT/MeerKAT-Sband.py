#!/usr/bin/python3

import casatasks as casa
import numpy as np
from casacore import tables as t
import os
import sys
import logging

os.system('irainit load casapy release-5.5.0-149.el7')

# Older CASA
from recipes.almapolhelpers import *
# Current CASA
#from almatasks import *

invis   = 'RawData/m87sband-flipped.MS'
calms   = 'MS_Files/m87sband-cal.MS'
tgtms   = 'MS_Files/m87sband-tgt.MS'
tgtavgms   = 'MS_Files/m87sband-tgt-avg.MS'
ref_ant = 'm002'

# Name your gain tables
tab = {'K_tab' : 'delay.cal',
       'B_tab' : 'bandpass.cal',
       'Gp_tab': 'gain_p.cal',
       'Ga_tab': 'gain_a.cal',
       'Tsec_tab' : 'T_sec.cal',
       'Gpsec_tab' : 'gain_p_sec.cal',
       # Pol cal tables #
       'Kcross_tab': 'kcross.cal',
       'Xf_tab': 'Xf.cal',
       'Df_tab': 'Df.cal'}
inpath = 'CASA_Tables/' # use ./ for $PWD

FluxCal = 'J1939-6342' # one of the BandPassCals
BandPassCal = 'J1939-6342,J0408-6545'
PolCal = 'J1331+3030'
PhaseTargetDic = {'J1150-0023':'M87'} # PhaseCal <--> Target pairs

# Fix some variables
for name in tab:
    tab[name] = os.path.join(inpath, tab[name])
PhaseCal = ','.join(PhaseTargetDic.keys())
Targets = ','.join(PhaseTargetDic.values())
CalibFields = ','.join([BandPassCal, PolCal, PhaseCal])

tricolour_command = f'singularity run --bind $(pwd) -B {os.getcwd()} ~/storage/tricolour.simg tricolour'

#############################
### Logs Setting up and functions
log_file = os.path.join(invis + '.log')
casa_log = os.path.join(invis + '_casa.log')

log_level = logging.DEBUG
logging.basicConfig(filename=log_file,
        format='%(asctime)s %(name)s:%(funcName)s\t%(message)s',
        datefmt='%Y-%m-%d %H:%M', filemode='w',
        level=log_level)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
logger = logging.getLogger(__name__)

#### NEEDED FOR J0408-6545 from https://skaafrica.atlassian.net/wiki/spaces/ESDKB/pages/1481408634/Flux+and+bandpass+calibration
def casa_flux_model(lnunu0, iref, *args):
    """
    Compute model:
    iref * 10**lnunu0 ** (args[0] + args[1] * lnunu0 + args[1] * lnunu0 ** 2 + args[0] * lnunu0 ** 3)
    """
    exponent = np.sum([arg * (lnunu0 ** (power))
                       for power, arg in enumerate(args)], axis=0)
    return iref * (10**lnunu0) **(exponent)

def fit_flux_model(nu, s, nu0, sigma, sref, order=5):
    """
    Fit a flux model of given order from :
    S = fluxdensity *(freq/reffreq)**(spix[0]+spix[1]*log(freq/reffreq)+..)
    Very rarely, the requested fit fails, in which case fall
    back to a lower order, iterating until zeroth order. If all
    else fails return the weighted mean of the components.
    Finally convert the fitted parameters to a
    katpoint FluxDensityModel:
    log10(S) = a + b*log10(nu) + c*log10(nu)**2 + ...
    Parameters
    ----------
    nu : np.ndarray
        Frequencies to fit in Hz
    s : np.ndarray
        Flux densities to fit in Jy
    nu0 : float
        Reference frequency in Hz
    sigma : np.ndarray
        Errors of s
    sref : float
        Initial guess for the value of s at nu0
    order : int (optional)
        The desired order of the fitted flux model (1: SI, 2: SI + Curvature ...)
    """
    from scipy.optimize import curve_fit

    init = [sref, -0.7] + [0] * (order - 1)
    lnunu0 = np.log10(nu/nu0)
    for fitorder in range(order, -1, -1):
        try:
            popt, _ = curve_fit(casa_flux_model, lnunu0, s, p0=init[:fitorder + 1], sigma=sigma)
        except RuntimeError:
            log.warn("Fitting flux model of order %d to CC failed. Trying lower order fit." %
                     (fitorder,))
        else:
            coeffs = np.pad(popt, ((0, order - fitorder),), "constant")
            return [nu0] +  coeffs.tolist()
    # Give up and return the weighted mean
    coeffs = [np.average(s, weights=1./(sigma**2))] + [0] * order
    return [nu0]+  coeffs.tolist()

def convert_flux_model(nu=np.linspace(0.9,2,200)*1e9 , a=1,b=0,c=0,d=0,Reffreq= 1.0e9) :
    """
    Convert a flux model from the form:
    log10(S) = a + b*log10(nu) + c*log10(nu)**2 + ...
    to an ASA style flux model in the form:
    S = fluxdensity *(freq/reffreq)**(spix[0]+spix[1]*log(freq/reffreq)+..)
    Parameters
    ----------
    nu : np.ndarray
        Frequencies to fit in Hz
    a,b,c,d : float
        parameters of a log flux model.
    Reffreq : float
        Reference frequency in Hz
    returns :
    reffreq,fluxdensity,spix[0],spix[1],spix[2]
    """
    MHz = 1e6
    S = 10**(a + b*np.log10(nu/MHz) +c*np.log10(nu/MHz)**2 + d*np.log10(nu/MHz)**3)
    return fit_flux_model(nu, S , Reffreq,np.ones_like(nu),sref=1 ,order=3)

##############################
# Change RECEPTOR_ANGLE : DEFAULT IS -90DEG but should be fixed with the initial swap
casa.tb.open(invis+'/FEED', nomodify=False)
feed_angle = casa.tb.getcol('RECEPTOR_ANGLE')
new_feed_angle = np.zeros(feed_angle.shape)
casa.tb.putcol('RECEPTOR_ANGLE', new_feed_angle)
casa.tb.close()

###########################
# Split the calibrators
casa.split(vis = invis, outputvis = calms, field = f"{BandPassCal},{PolCal},{PhaseCal}", datacolumn = 'data', spw = '0:210~3841')

# Standard flagging for shadowing, zero-clip, and auto-correlation
casa.flagdata(vis=calms, flagbackup=False, mode='shadow')
casa.flagdata(vis=calms, flagbackup=False, mode='manual', autocorr=True)
casa.flagdata(vis=calms, flagbackup=False, mode='clip', clipzeros=True, clipminmax=[0.0, 100.0])

# Set flux density scale
for cal in set(FluxCal.split(',')+BandPassCal.split(',')+PolCal.split(',')):

    if cal == 'J1939-6342':
        casa.setjy(vis = calms, field = cal, standard = 'Stevens-Reynolds 2016', usescratch = True)
    elif cal == 'J0408-6545':
        a=-0.9790; b=3.3662; c=-1.1216; d=0.0861
        reffreq,fluxdensity,spix0,spix1,spix2 =  convert_flux_model(np.linspace(0.9,2,200)*1e9,a,b,c,d)
        casa.setjy(vis = calms, field = cal, usescratch = True, \
            spix = [spix0, spix1, spix2, 0], fluxdensity = fluxdensity, reffreq = '%f Hz'%(reffreq), standard = 'manual', )
    elif cal == 'J1331+3030':
            reffreq = '1.284GHz'
            # Stokes flux density in Jy
            I =  15.74331687; Q = 0.8628247336; U = 1.248991241; V = 0
            # Spectral Index in 2nd order
            alpha =   [-0.4881438633844231, -0.17025582235426978]
            casa.setjy(vis=calms, field=cal, standard="manual", fluxdensity=[I,Q,U,V], spix=alpha, reffreq=reffreq,
                polindex=[0.0964244115642966, 0.018277345381024372, -0.07332409550519345, 0.3253188415787851, -0.37228554528542385],
                polangle=[0.48312994184873537, 0.12063061722082152, -0.04180094935296229, 0.12832951565823608],
                rotmeas=0.12,usescratch=True)
    else:
        print("Unknown calibrator ", cal)
        sys.exit()

################################################################
# backup flags
casa.flagmanager(vis = calms, mode = 'save', versionname = 'PreCal')

# Run TFCrop on the DATA column
casa.flagdata(vis=calms, mode='tfcrop', field=CalibFields,
        ntime='scan', timecutoff=5.0, freqcutoff=5.0, timefit='line',
        freqfit='line', extendflags=False, timedevscale=5., freqdevscale=5.,
        extendpols=True, growaround=False, action='apply', flagbackup=False,
        overwrite=True, writeflags=True, datacolumn='DATA')
casa.flagdata(vis=calms, mode='extend', field=CalibFields,
        datacolumn='data', clipzeros=True, ntime='scan', extendflags=False,
        extendpols=True, growtime=80., growfreq=80., growaround=False,
        flagneartime=False, flagnearfreq=False, action='apply',
        flagbackup=False, overwrite=True, writeflags=True)

### Basic calibration
# Delay calibration
casa.gaincal(vis=calms, field=BandPassCal, caltable=tab['K_tab'], gaintype='K', refant=ref_ant)
# Gani calibration
casa.gaincal(vis=calms, field=BandPassCal, caltable=tab['Gp_tab'], gaintype='G', calmode='p', gaintable=[tab['K_tab']], refant=ref_ant)
casa.gaincal(vis=calms, field=BandPassCal, caltable=tab['Ga_tab'], gaintype='G', calmode='a', gaintable=[tab['K_tab'],tab['Gp_tab']], refant=ref_ant)
# one can now combine the scans and use different B as diagnostics
casa.bandpass(vis=calms, field=BandPassCal, caltable=tab['B_tab'], bandtype='B', gaintable=[tab['K_tab'],tab['Gp_tab'],tab['Ga_tab']], combine='scan', refant=ref_ant)

# Restore original falgs
casa.flagmanager(vis=calms, mode='restore', versionname='PreCal')

# applycal
casa.applycal(vis=calms,field=CalibFields, gaintable=[tab['K_tab'],tab['Gp_tab'],tab['Ga_tab'],tab['B_tab']])

# Flag with tricolour


### Redo basic calibtration
# Delay calibration
casa.gaincal(vis=calms, field=BandPassCal, caltable=tab['K_tab'], gaintype='K', refant=ref_ant)
# PLOT: plotms(vis=tab['K_tab']) 
# Gani calibration
casa.gaincal(vis=calms, field=BandPassCal, caltable=tab['Gp_tab'], gaintype='G', calmode='p', gaintable=[tab['K_tab']], refant=ref_ant)
casa.gaincal(vis=calms, field=BandPassCal, caltable=tab['Ga_tab'], gaintype='G', calmode='a', gaintable=[tab['K_tab'],tab['Gp_tab']], refant=ref_ant)
# one can now combine the scans and use different B as diagnostics
casa.bandpass(vis=calms, field=BandPassCal, caltable=tab['B_tab'], bandtype='B', gaintable=[tab['K_tab'],tab['Gp_tab'],tab['Ga_tab']], combine='scan', refant=ref_ant)
# Leackage
casa.polcal(vis=calms,
   caltable=tab['Df_tab'],field=FluxCal, poltype='Df', solint='inf', refant=ref_ant,
   gaintable=[tab['K_tab'], tab['Gp_tab'], tab['Ga_tab'], tab['B_tab']])

############################################################################
# Bootrap secondary calibrator
casa.gaincal(vis=calms, caltable='sec.Gp', field='J1733-1304', gaintype='G', calmode='p', solint='inf', combine='',refant=ref_ant, gaintable=['bp.K', 'bp.Ga', 'bp.B', 'bp.Df'])
casa.gaincal(vis=calms, caltable='sec.T', field='J1733-1304', gaintype='T', calmode='a', solnorm=True, solint='inf', combine='',refant=ref_ant, gaintable=['bp.K', 'sec.Gp', 'bp.Ga', 'bp.B', 'bp.Df'])
#image the secondary and selfcal
casa.tclean(vis=calms,field=PhaseCal,cell='0.5arcsec',imsize=512,niter=1000,imagename=PhaseCal+'-selfcal',weighting='briggs',robust=-0.2,
       datacolumn= 'corrected',deconvolver= 'mtmfs',nterms=2,specmode='mfs',interactive=False)
gaincal(vis=calms, caltable='sec.Gp', field='J1733-1304', gaintype='G', calmode='p', solint='inf', combine='',refant=ref_ant, gaintable=['bp.K', 'bp.Ga', 'bp.B', 'bp.Df'])
# QUESTION: can I remove solnorm here?
gaincal(vis=calms, caltable='sec.T', field='J1733-1304', gaintype='T', calmode='a', solnorm=True, solint='inf', combine='',refant=ref_ant, gaintable=['bp.K', 'sec.Gp', 'bp.Ga', 'bp.B', 'bp.Df'])

# Solve for polarization angle
# KCROSS
casa.gaincal(vis=calms,
    caltable=tab['Kcross_tab'], field=PolCal, gaintype='KCROSS', refant=ref_ant, solint='inf', parang=True,
    gaintable=[tab['K_tab'], tab['Gp_tab'], tab['Ga_tab'], tab['B_tab'], tab['Df_tab']])
# Xf that is constant within a scan
casa.polcal(vis=msFile,
   caltable=tab['Xf_tab'], field=PolCal, poltype='Xf',
   uvrange=myuvrange, solint='inf', refant=ref_ant,
   combine='scan', preavg=-1.,
   gaintable=[tab['B_tab'], tab['Ga_tab'], tab['K_tab'], tab['Kcross_tab'], tab['Df_tab'], tab['Tsec_tab'], tab['Gpsec_tab']])

# Final applycal to all sources
casa.applycal(vis=calms, field='*', gaintable=['bp.K','sec.Gp','bp.Ga','bp.B','sec.T','bp.Df','pol.Xf'], parang=True)

# test images
tclean(vis=calms,field=xcal,cell='0.5arcsec',imsize=512,niter=1000,imagename=xcal+'-selfcal',weighting='briggs',robust=-0.2,datacolumn= 'corrected',deconvolver= 'mtmfs',\
       nterms=2,specmode='mfs',interactive=False)

###############################################################################
# Target

applycal(vis  = calms, parang = True, calwt = False, field = '',\
    gaintable = [ktab, btab, kxtab, ptab_xyf, dgen, ftab],\
    gainfield = ['', '', '', '', '', ''],\
    interp    = ['nearest,linear','nearest,linearflag','nearest,linear','nearest,linearflag','linear,linearflag','linear,linear'])


wsclean -verbose -log-time -no-update-model-required -j 64 \
    -field 2 -weight briggs -0.5 -size 2048 2048 -scale 1.0asec -channels-out 48 \
    -no-mf-weighting -weighting-rank-filter 3 -auto-mask 5 -auto-threshold 1.0 \
    -taper-gaussian 6 -pol IQUV -data-column CORRECTED_DATA -niter 10000 -gain 0.05 \
    -mgain 0.9 -join-polarizations -join-channels -squared-channel-joining -gridder wgridder \
    -padding 1.3 -local-rms -name img/3c286_wsclean_FullStokes MS_Files/m87sband-cal.MS

# IMAGE 3C138... OTHER PYTHON FILE

### WORK ON TARGET ...
# Split target ...
split(vis = invis, outputvis = tgtms, field = "M87", datacolumn = 'data', spw = '0:210~3841')

# Change RECEPTOR_ANGLE : DEFAULT IS -90DEG but should be fixed with the initial swap
casa.tb.open(calms+'/FEED', nomodify=False)
feed_angle = casa.tb.getcol('RECEPTOR_ANGLE')
new_feed_angle = np.zeros(feed_angle.shape)
casa.tb.putcol('RECEPTOR_ANGLE', new_feed_angle)
casa.tb.close()

applycal(vis  = tgtms, parang = True, calwt = False, field = '',\
    gaintable = [ktab, btab, kxtab, ptab_xyf, dgen, ftab],\
    gainfield = ['', '', '', '', '', gcal],\
    interp    = ['nearest,linear','nearest,linearflag','nearest,linear','nearest,linearflag','linear,linearflag','linear,linear'])

flagmanager(vis = tgtms, mode = 'save', versionname = 'ApplyCal')

aoflagger-setup
aoflagger -v -j 32 -strategy meerkat_custom20230417.lua -column CORRECTED_DATA MS_Files/xxx.MS

split(vis = tgtms, outputvis = tgtavgms, datacolumn = 'corrected', width = 8)
flagdata(vis=tgtavgms, mode='manual', spw = '0:66~72,0:106~114,0:202~210,0:278~285')

# Casa
for i in range(30):
    tgtavgms   = 'MS_Files/m87sband-tgt-avg.MS'
    gaincal(vis=tgtavgms, caltable='selfcal%02i.G' %i, solint='8s', refant='m002', parang=False)
    gaincal(vis=tgtavgms, caltable='selfcal%02i.K' %i, solint='8s', refant='m002', gaintype='K', gaintable='selfcal%02i.G' %i, parang=False)
    bandpass(vis=tgtavgms, caltable='selfcal%02i.B' %i, combine='', solint='300s', gaintable=['selfcal%02i.G' %i, 'selfcal%02i.K' %i], refant='m002', parang=False)
    applycal(vis=tgtavgms, gaintable=['selfcal%02i.B' %i, 'selfcal%02i.G' %i, 'selfcal%02i.K' %i], interp=['linear,linearflag','linear', 'linear,linearflag'], parang=False)
    os.system('singularity exec ~/storage/pill.simg wsclean -name img/m87-test%02i -reorder -parallel-reordering 5 -parallel-gridding 12 -j 64 -mem 100 -update-model-required -weight briggs 0.0 -size 2500 2500 -scale 0.7arcsec -channels-out 454 -deconvolution-channels 8 -pol XX,YY -data-column CORRECTED_DATA -niter 10000000 -auto-threshold 2 -gain 0.1 -mgain 0.5 -join-channels -multiscale -fit-spectral-pol 3 -multiscale -no-mf-weighting -fits-mask m87-07asec-2500.fits MS_Files/m87sband-tgt-avg.MS/' % i)


# DP3
import os, glob
mss = sorted(glob.glob('MS_Files/m87sband-tgt-full-scan*.MS'))
for i in range(30):
    print(f'Cycle: {i}')
    for ms in mss:
        print(f'Working on {ms}...')
        # solve
        os.system(f'DP3 DP3-sol.parset msin={ms} msout=. sol.h5parm={ms}/ph-{i}.h5 sol.mode=diagonalphase sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=10e6 >> DP3.log')
        os.system(f'DP3 DP3-sol.parset msin={ms} msout=. sol.h5parm={ms}/amp-{i}.h5 sol.mode=diagonalamplitude sol.solint=20 sol.nchan=1 sol.smoothnessconstraint=30e6 >> DP3.log')
        # correct
        os.system(f'DP3 DP3-cor.parset msin={ms} msout=. cor1.parmdb={ms}/ph-{i}.h5 cor2.parmdb={ms}/amp-{i}.h5 >> DP3.log')
    # clean
    print(f'Cleaning...')
    imgname = "img/m87-dp3%02i" % i
    os.system("rm -r wsclean_concat.MS")
    os.system(f"taql select from {mss} giving wsclean_concat.MS as plain")
    #os.system(f'wsclean -name {imgname} -reorder -parallel-reordering 5 -parallel-gridding 12 -j 64 -mem 100 -no-update-model-required -weight briggs 0.0 -size 2500 2500 -scale 0.7arcsec -channels-out 454 -deconvolution-channels 8 -pol XX,YY,XY,YX -data-column CORRECTED_DATA -niter 10000000 -auto-threshold 2 -gain 0.1 -mgain 0.5 -join-channels -multiscale -fit-spectral-pol 3 -multiscale -no-mf-weighting -fits-mask m87-07asec-2500.fits wsclean_concat.MS >> wsclean.log')
    os.system(f'wsclean -name {imgname} -reorder -parallel-reordering 5 -parallel-gridding 12 -j 64 -mem 100 -no-update-model-required -weight briggs 0.0 -size 2500 2500 -scale 0.7arcsec -channels-out 454 -deconvolution-channels 8 -pol IQUV -data-column CORRECTED_DATA -niter 10000000 -auto-threshold 2 -gain 0.1 -mgain 0.5 -join-channels -multiscale -fit-spectral-pol 3 -multiscale -no-mf-weighting -join-polarizations -fits-mask m87-07asec-2500.fits wsclean_concat.MS >> wsclean.log')
    for ms in mss:
        print(f'Predict on {ms}...')
        # predict
        os.system(f'wsclean  -predict -padding 1.8 -j 64 -name {imgname} -channels-out 461 {ms} >> wsclean.log')
