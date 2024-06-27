##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### Full-polarisation calibration script for MeerKAT S-band data
##### Original version C.J.Riseley (June 2023) inspired by techniques developed and discussed
##### in Riseley et al. (2015) and Riseley (2016)
#####

casapy-setup release-5.5.0-149.el7

# Older CASA
from recipes.almapolhelpers import *
# Current CASA
#from almatasks import *
invis   = 'RawData/m87sband.MS'
calms   = 'MS_Files/m87sband-cal.MS'
tgtms   = 'MS_Files/m87sband-tgt.MS'
tgtavgms   = 'MS_Files/m87sband-tgt-avg.MS'
ref_ant = 'm002'

# Name your gain tables
gtab     = "CASA_Tables/calib.gcal0"
ktab     = "CASA_Tables/calib.kcal"
gtab_all = "CASA_Tables/calib.gcal_0"
gtab_pol = "CASA_Tables/calib.gcal_1"
btab     = "CASA_Tables/calib.bandpass"
ftab     = "CASA_Tables/calib.fluxscale"
kxtab    = "CASA_Tables/calib.xdelcal"
ptab_amb = "CASA_Tables/calib.xy0amb"
ptab_xyf = "CASA_Tables/calib.xyf"
dtab     = "CASA_Tables/calib.leak0"
dgen     = "CASA_Tables/calib.leakgen"

fcal  = 'J1939-6342'
bpcal = 'J1939-6342,J0408-6545'
gcal  = 'J1150-0023' #'J0217+017'
xcal  = 'J1331+3030' #'J0521+1638'

# SPLIT OUT CALIBRATORS ...
split(vis = invis, outputvis = calms, field = "{0},{1},{2}".format(bpcal, gcal, xcal),\
    datacolumn = 'data', spw = '0:200~3838')

# Shadow flagging:
flagdata(vis = calms, mode = 'shadow')

# SPW flagging : in S1 sub-band, only edge channels need flagging... learn lessons from L-band channel edge flagging:
#flagdata(vis = calms, mode = 'manual', spw = '0:0~256,0:3838~4095', flagbackup = True)
## Dual RFI spike near 2 GHz // channels 180-190 approx // seems most visible in cross-hands, but still visible in parallel hands.
## Now taken care of by split command

# RFI flagging on parallel hands with TFCROP
flagdata(vis=calms, mode="tfcrop", autocorr=False, inpfile="", reason="any",\
        tbuff=0.0, spw="", field="*", clipoutside=True, correlation='ABS_XX,YY',\
        channelavg=False, chanbin=1, timeavg=False, timebin="0s", clipzeros=True,\
        ntime="scan", combinescans=False, timecutoff=4.0, freqcutoff=3.0,\
        timefit="line", freqfit="poly", maxnpieces=7, flagdimension="freqtime",\
        usewindowstats="std", halfwin=1, extendflags=True, winsize=3, timedev="",\
        freqdev="", timedevscale=5.0, freqdevscale=5.0, spectralmax=1000000.0, spectralmin=0.0,\
        antint_ref_antenna="", minchanfrac=0.6, verbose=False, extendpols=True, growtime=50.0,\
        growfreq=50.0, growaround=False, flagneartime=False, flagnearfreq=False, minrel=0.0,\
        maxrel=1.0, minabs=0, maxabs=-1, spwchan=False, spwcorr=False,\
        basecnt=False, fieldcnt=False, name="Summary", action="apply", display="",\
        flagbackup=True, savepars=False, cmdreason="", outfile="", overwrite=True,\
        writeflags=True)

# RFI flagging on cross hands with RFLAG
flagdata(vis=calms, mode="rflag", spw="", field="*", correlation='XY,YX',\
        timedevscale=4.0, freqdevscale=3.0,\
        extendflags=True, winsize=3, timedev="", freqdev="",\
        antint_ref_antenna="", minchanfrac=0.6, verbose=False,\
        name="Summary", action="apply", display="", flagbackup=True,\
        overwrite=True, writeflags=True)

# GAIN CALIB RFI flagging on parallel hands with RFLAG
flagdata(vis=calms, mode="rflag", spw="", field=gcal, correlation='XX,YY',\
        timedevscale=4.0, freqdevscale=3.0,\
        extendflags=True, winsize=3, timedev="", freqdev="",\
        antint_ref_antenna="", minchanfrac=0.6, verbose=False,\
        name="Summary", action="apply", display="", flagbackup=True,\
        overwrite=True, writeflags=True)

#### MANUAL NOTES ON RFI / DATA QA : CBID 1688876155
## J0217+017
## SCAN 6 SOME RESIDUAL RFI AROUND CHAN 750 + LOW AMPLITUDES AROUND CHAN 2450
## Dual RFI spike near 2 GHz // channels 180-190 approx // seems most visible in cross-hands, but still visible in parallel hands.
## NOW TAKEN CARE OF BY ABOVE CALLS TO RFLAG

# Change RECEPTOR_ANGLE : DEFAULT IS -90DEG ... MAY LEAVE AS IT IS?
#tb.open(calms+'/FEED', nomodify=False)
#feed_angle = tb.getcol('RECEPTOR_ANGLE')
#new_feed_angle = np.zeros(feed_angle.shape)
#tb.putcol('RECEPTOR_ANGLE', new_feed_angle)
#tb.close()

# Clear previous calibration and restore flags
#clearcal(calms)
#flagmanager(vis = calms, mode = 'restore', versionname = 'PreCal')
flagmanager(vis = calms, mode = 'save', versionname = 'PreCal')

# Set flux density scale
default("setjy")

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
    from scipy.optimize import curve_fit
    from scipy.special import binom
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

# J0408-6545
a=-0.9790
b=3.3662
c=-1.1216
d=0.0861
################################################################

for cal in set(fcal.split(',')+bpcal.split(',')):

    if cal == 'J1939-6342':
        setjy(vis = calms, field = "{0}".format(fcal), spw = "", selectdata = False, timerange = "", scan = "", \
            standard = 'Stevens-Reynolds 2016', scalebychan = True, useephemdir = False, usescratch = True)
    elif cal == 'J0408-6545'
        reffreq,fluxdensity,spix0,spix1,spix2 =  convert_flux_model(np.linspace(0.9,2,200)*1e9,a,b,c,d)
        setjy(vis = calms, field = bpcal.split(',')[1], usescratch = True, scalebychan = True,\
            spix = [spix0, spix1, spix2, 0], fluxdensity = fluxdensity, reffreq = '%f Hz'%(reffreq), standard = 'manual', )
    else:
        print("Unknown calibrator ", cal)
        sys.exit()

# Initial phase calibration
gaincal(vis = calms, caltable = gtab, selectdata = True,\
    solint = "inf", field = bpcal, refant = ref_ant,\
    solnorm = False, gaintype = "G", calmode = "p",
    minsnr = 3.0, parang = False)

# Delay calibration
gaincal(vis = calms, caltable = ktab, selectdata = True,\
    solint = "inf", field = bpcal, combine = "",\
    refant = ref_ant, solnorm = False, gaintype = "K",\
    gaintable = [gtab], gainfield = [''], interp = ['linear'],\
    parang = False)

# Bandpass calibration
bandpass(vis = calms, caltable = btab, selectdata = True,\
    solint = "inf", field = bpcal, combine = "",\
    refant = ref_ant, solnorm = False, bandtype = "B",\
    gaintable = [ktab], gainfield = [''],\
    interp = ['linear'], parang = False)

# Gain calibration: bandpass calibrators
#gaincal(vis = calms, caltable = gtab_all, selectdata = True,\
#    solint = "20s", field = bpcal, combine = "",\
#    refant = ref_ant, gaintype = "G", calmode = "ap",\
#    gaintable = [ktab, btab], gainfield = ['', ''],\
#    interp = ['nearest', 'linear,linear'],\

gaincal(vis = calms, caltable = gtab_all, selectdata = True,\
    solint = "20s", field = bpcal.split(',')[0], combine = "",\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linearflag'],\
    parang = False)

gaincal(vis = calms, caltable = gtab_all, selectdata = True,\
    solint = "20s", field = bpcal.split(',')[1], combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linearflag'],\
    parang = False)

# Gain calibration: gain calibrator
gaincal(vis = calms, caltable = gtab_all, selectdata = True,\
    solint = "20s", field = gcal, combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linearflag'],\
    parang = False)

# Gain calibration: polarisation calibrator
gaincal(vis = calms, caltable = gtab_all, selectdata = True,\
    solint = "20s", field = xcal, combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linearflag'],\
    parang = False)

### Polarisation calibration
qu = qufromgain(gtab_all, paoffset=-90)

#Fields: 4
#ID   Code Name
#0    T    J0408-6545   
#1    T    J1150-0023   
#2    T    J1331+3030
#3    T    J1939-6342   

"""
Latitude =  -30.7120094176
Found as many as 4 fields.
Can't discern an ALMA bandname from: none
Found as many as 1 spws.
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
Fld= 0 Spw= 0 Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
(B=none, PA offset=-90.0deg) Gx/Gy= 0.926416138209 Q= 0.0395901298801 U= -0.0029847551957 P= 0.0397024828884 X= -2.15572631543
For field id =  0  there are  1 good spws.
Spw mean: Fld= 0 Q= 0.0395901298801 U= -0.0029847551957 (rms= 0.0 0.0 ) P= 0.0397024828884 X= -2.15572631543
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
Fld= 1 Spw= 0 Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
(B=none, PA offset=-90.0deg) Gx/Gy= 1.00263564675 Q= 0.000391777676417 U= -0.0417012260384 P= 0.0417030663484 X= -44.7308646638
For field id =  1  there are  1 good spws.
Spw mean: Fld= 1 Q= 0.000391777676417 U= -0.0417012260384 (rms= 0.0 0.0 ) P= 0.0417030663484 X= -44.7308646638
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
Fld= 2 Spw= 0 Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to -1.0
(B=none, PA offset=-90.0deg) Gx/Gy= 1.02356942122 Q= 0.0610014399346 U= 0.0789237847053 P= 0.0997503857952 X= 26.1494984528
For field id =  2  there are  1 good spws.
Spw mean: Fld= 2 Q= 0.0610014399346 U= 0.0789237847053 (rms= 0.0 0.0 ) P= 0.0997503857952 X= 26.1494984528
Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
Fld= 3 Spw= 0 Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
(B=none, PA offset=-90.0deg) Gx/Gy= 1.00051184435 Q= -0.000301214682363 U= -7.41574149283e-05 P= 0.00031020897321 X= -83.0845749062
For field id =  3  there are  1 good spws.
Spw mean: Fld= 3 Q= -0.000301214682363 U= -7.41574149283e-05 (rms= 0.0 0.0 ) P= 0.00031020897321 X= -83.0845749062
"""

# Cross-hand delay calibration
gaincal(vis = calms, caltable = kxtab, selectdata = True,\
    solint = "inf", field = xcal, combine = "scan",\
    refant = ref_ant, gaintype = "KCROSS",\
    smodel=[1.0, 0.0, 1.0, 0],\
    gaintable = [ktab, btab, gtab_all], gainfield = ['', '', xcal],\
    interp = ['nearest', 'linear,linearflag', 'linear'],\
    parang = False)

# -0.0078nsec

# XY phase calibration
gaincal(vis = calms, caltable = ptab_amb, selectdata = True,\
    solint = "inf", field = xcal, combine = "scan",\
    refant = ref_ant, gaintype = "XYf+QU", calmode = "ap",\
    smodel = [1.0, 0.0, 1.0, 0], preavg = 200.0, \
    gaintable = [ktab, btab, gtab_all, kxtab], gainfield = ['', '', xcal, xcal],\
    interp = ['nearest', 'linear,linearflag', 'linear', 'linear'],\
    parang = False)

"""
Spw = 0 (ich=1819/3639): 
X-Y phase = 15.9441255136 deg.
Fractional Poln: Q = 0.0486539267004, U = 0.0800458192825; P = 0.0936725037009, X = 29.3538407165deg.
Net (over baselines) instrumental polarization: -0.00793615696959
"""

S = xyamb(xytab = ptab_amb, qu = qu[2], xyout = ptab_xyf)

"""
Expected QU =  (0.061001439934595385, 0.07892378470526927)
Spw = 0: Found QU = [ 0.04865393  0.08004582]
      ...KEEPING X-Y phase 13.7561486581 deg
Ambiguity resolved (spw mean): Q= 0.0486539267004 U= 0.0800458192825 (rms= 0.0 0.0 ) P= 0.0936725027315 X= 29.3538398993
Returning the following Stokes vector: [1.0, 0.048653926700353622, 0.080045819282531738, 0.0]
"""

plotcal(caltable = ptab_amb, xaxis = 'time', yaxis = 'amp', poln = '/', subplot = 111, field = '2',\
    showgui = False, figfile = 'gxgy_vs_parang.png')

plotcal(caltable = ptab_amb, xaxis = 'chan', yaxis = 'phase', subplot = 111, field = '2',\
    showgui = False, figfile = 'xy_phase.png')

# RE-DERIVE GAIN TABLES
gaincal(vis = calms, caltable = gtab_pol, selectdata = True,\
    solint = "20s", field = bpcal.split(',')[0], combine = "",\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linearflag'],\
    parang = True)

gaincal(vis = calms, caltable = gtab_pol, selectdata = True,\
    solint = "20s", field = bpcal.split(',')[1], combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linearflag'],\
    parang = True)

# Gain calibration: gain calibrator
gaincal(vis = calms, caltable = gtab_pol, selectdata = True,\
    solint = "20s", field = gcal, combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linearflag'],\
    parang = True)

# Gain calibration: polarisation calibrator
gaincal(vis = calms, caltable = gtab_pol, selectdata = True,\
    solint = "20s", field = xcal, combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap", smodel = S,\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linearflag'],\
    parang = True)

qu = qufromgain(gtab_pol, paoffset = -90)
# now the P for the pol cal should be ~0

"""
Fld= 2 Spw= 0 Can't discern an ALMA bandname from: none
Unresolved bandname: default band position angle set to 0.0
(B=none, PA offset=-90.0deg) Gx/Gy= 0.999482530243 Q= 0.00479459789233 U= 0.00731266016854 P= 0.00874432202573 X= 28.374430667
For field id =  2  there are  1 good spws.
Spw mean: Fld= 2 Q= 0.00479459789233 U= 0.00731266016854 (rms= 0.0 0.0 ) P= 0.00874432202573 X= 28.374430667
"""

# Leakage calibration
polcal(vis = calms, caltable = dtab, selectdata = True,\
    solint = 'inf', field = xcal, combine = 'scan',\
    refant = '', poltype = 'Dflls', smodel = S, preavg = 200,\
    gaintable = [ktab, btab, gtab_pol, kxtab, ptab_xyf],\
    gainfield = ['', '', xcal, xcal, xcal],\
    interp = ['nearest', 'linear,linearflag', 'nearest', 'nearest', 'nearest,linearflag'])

# Generalise leakage solution to cross-hands
Dgen(dtab = dtab, dout = dgen)

# Flux scale correction
fluxscale(vis = calms, caltable = gtab_pol, reference = [fcal],\
    transfer = [bpcal.split(',')[1], gcal, xcal], fluxtable = ftab,\
    append = False, incremental = False)

### APPLY TO CALIBRATORS AND IMAGE 3C138 ETC ...
flagmanager(vis = calms, mode = 'save', versionname = 'AfterCal')
#flagmanager(vis = calms, mode = 'restore', versionname = 'AfterCal')

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
split(vis = invis, outputvis = tgtms, field = "M87", datacolumn = 'data', spw = '0:200~3838')

applycal(vis  = tgtms, parang = True, calwt = False, field = '',\
    gaintable = [ktab, btab, kxtab, ptab_xyf, dgen, ftab],\
    gainfield = ['', '', '', '', '', gcal],\
    interp    = ['nearest,linear','nearest,linearflag','nearest,linear','nearest,linearflag','linear,linearflag','linear,linear'])

flagmanager(vis = tgtms, mode = 'save', versionname = 'ApplyCal')

aoflagger-setup
aoflagger -v -j 32 -strategy meerkat_custom20230417.lua -column CORRECTED_DATA MS_Files/xxx.MS

split(vis = tgtms, outputvis = tgtavgms, datacolumn = 'corrected', width = 8)
flagdata(vis=tgtavgms, mode='manual', spw = '0:66~72,0:106~114,0:202~210,0:278~285')

for i in range(30):
    tgtavgms   = 'MS_Files/m87sband-tgt-avg.MS'
    gaincal(vis=tgtavgms, caltable='selfcal%02i.G' %i, solint='8s', refant='m002', parang=False)
    gaincal(vis=tgtavgms, caltable='selfcal%02i.K' %i, solint='8s', refant='m002', gaintype='K', gaintable='selfcal%02i.G' %i, parang=False)
    bandpass(vis=tgtavgms, caltable='selfcal%02i.B' %i, combine='', solint='300s', gaintable=['selfcal%02i.G' %i, 'selfcal%02i.K' %i], refant='m002', parang=False)
    applycal(vis=tgtavgms, gaintable=['selfcal%02i.B' %i, 'selfcal%02i.G' %i, 'selfcal%02i.K' %i], interp=['linear,linearflag','linear', 'linear,linearflag'], parang=False)
    os.system('singularity exec ~/storage/pill.simg wsclean -name img/m87-test%02i -reorder -parallel-reordering 5 -parallel-gridding 12 -j 64 -mem 100 -update-model-required -weight briggs 0.0 -size 2500 2500 -scale 0.7arcsec -channels-out 454 -deconvolution-channels 8 -pol I -data-column CORRECTED_DATA -niter 10000000 -auto-threshold 2 -gain 0.1 -mgain 0.5 -join-channels -multiscale -fit-spectral-pol 5 -multiscale -no-mf-weighting -fits-mask m87-07asec-2500.fits MS_Files/m87sband-tgt-avg.MS/' % i)
