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
ref_ant = 'm052'

# Name your gain tables
gtab     = "CASA_Tables/calib.gcal0"
ktab     = "CASA_Tables/calib.kcal"
gtab_all = "CASA_Tables/calib.gcal_0"
gtab_pol = "CASA_Tables/calib.gcal_1"
btab     = "CASA_Tables/calib.bandpass"
ftab     = "CASA_Tables/calib.fluxscale"
kxtab    = "CASA_Tables/calib.xdelcal"
ptab_amb = "CASA_Tables/calib.xy0amb"
dtab     = "CASA_Tables/calib.leak0"
dgen     = "CASA_Tables/calib.leakgen"

fcal  = 'J1939-6342'
bpcal = 'J1939-6342,J0408-6545'
gcal  = 'J1150-0023' #'J0217+017'
xcal  = 'J1331+3030' #'J0521+1638'

#Fields: 4
#ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#0    T    J1939-6342          19:39:25.030000 -63.42.45.60000 J2000   0         117040
#1    T    J0217+017           02:17:48.950000 +01.44.49.70000 J2000   1         260260
#2    T    J0521+1638          05:21:09.890000 +16.38.22.10000 J2000   2         160160
#3    T    J0408-6545          04:08:20.380000 -65.45.09.10000 J2000   3         113960

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
clearcal(calms)
flagmanager(vis = calms, mode = 'restore', versionname = 'PreApplyCal')

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
    interp = ['nearest', 'linear,linear'],\
    parang = False)

gaincal(vis = calms, caltable = gtab_all, selectdata = True,\
    solint = "20s", field = bpcal.split(',')[1], combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linear'],\
    parang = False)

# Gain calibration: gain calibrator
gaincal(vis = calms, caltable = gtab_all, selectdata = True,\
    solint = "20s", field = gcal, combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linear'],\
    parang = False)

# Gain calibration: polarisation calibrator
gaincal(vis = calms, caltable = gtab_all, selectdata = True,\
    solint = "20s", field = xcal, combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linear'],\
    parang = False)

### Polarisation calibration
from almatasks import *
qu = qufromgain(gtab_all, paoffset=-90)

#Fields: 4
#ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#0    T    J1939-6342          19:39:25.030000 -63.42.45.60000 J2000   0         117040
#1    T    J0217+017           02:17:48.950000 +01.44.49.70000 J2000   1         260260
#2    T    J0521+1638          05:21:09.890000 +16.38.22.10000 J2000   2         160160
#3    T    J0408-6545          04:08:20.380000 -65.45.09.10000 J2000   3         113960


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### CBID 1688875155
# NB: default band position angle will be offset by -90deg.
# Latitude =  -30.7123828219
# Found as many as 4 fields.

# Fld= 0 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 1.02054763496 Q= -0.0100031147749 U= -0.00150779746044 P= 0.0101161137984 X= -85.7140882916
# For field id =  0  there are  1 good spws.
# Spw mean: Fld= 0 Q= -0.0100031147749 U= -0.00150779746044 (rms= 0.0 0.0 ) P= 0.0101161137984 X= -85.7140882916

# Fld= 1 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 1.0008047196 Q= -0.0151563213578 U= -0.00681893051843 P= 0.0166196236575 X= -77.8883501736
# For field id =  1  there are  1 good spws.
# Spw mean: Fld= 1 Q= -0.0151563213578 U= -0.00681893051843 (rms= 0.0 0.0 ) P= 0.0166196236575 X= -77.8883501736

# Fld= 2 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 1.00690416179 Q= 0.0820857901697 U= -0.0379964624936 P= 0.0904533476982 X= -12.4194074097
# For field id =  2  there are  1 good spws.
# Spw mean: Fld= 2 Q= 0.0820857901697 U= -0.0379964624936 (rms= 0.0 0.0 ) P= 0.0904533476982 X= -12.4194074097

# Fld= 3 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 0.974287701203 Q= -0.0125278883621 U= 0.00672615503328 P= 0.0142193230621 X= 75.8844455453
# For field id =  3  there are  1 good spws.
# Spw mean: Fld= 3 Q= -0.0125278883621 U= 0.00672615503328 (rms= 0.0 0.0 ) P= 0.0142193230621 X= 75.8844455453

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### CBID 1688961378
# NB: default band position angle will be offset by -90deg.
# Latitude =  -30.7124439492
# Found as many as 4 fields.

# Fld= 0 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 1.05722265854 Q= -0.0269933847427 U= -0.00207566159036 P= 0.0270730713239 X= -87.8014413766
# For field id =  0  there are  1 good spws.
# Spw mean: Fld= 0 Q= -0.0269933847427 U= -0.00207566159036 (rms= 0.0 0.0 ) P= 0.0270730713239 X= -87.8014413766

# Fld= 1 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 0.999414278364 Q= -0.0195204453057 U= -0.005181511351 P= 0.0201964314871 X= -82.5671038283
# For field id =  1  there are  1 good spws.
# Spw mean: Fld= 1 Q= -0.0195204453057 U= -0.005181511351 (rms= 0.0 0.0 ) P= 0.0201964314871 X= -82.5671038283

# Fld= 2 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 1.00002942461 Q= 0.0790099558732 U= -0.0411348850778 P= 0.089076663035 X= -13.751399364
# For field id =  2  there are  1 good spws.
# Spw mean: Fld= 2 Q= 0.0790099558732 U= -0.0411348850778 (rms= 0.0 0.0 ) P= 0.089076663035 X= -13.751399364

# Fld= 3 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 0.988333185697 Q= -0.00577434036083 U= 0.0015650451292 P= 0.00598267271871 X= 82.4175985735
# For field id =  3  there are  1 good spws.
# Spw mean: Fld= 3 Q= -0.00577434036083 U= 0.0015650451292 (rms= 0.0 0.0 ) P= 0.00598267271871 X= 82.4175985735

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

# Cross-hand delay calibration
gaincal(vis = calms, caltable = kxtab, selectdata = True,\
    solint = "inf", field = xcal, combine = "scan",\
    refant = ref_ant, gaintype = "KCROSS",\
    smodel=[1, 0.0, 1.0, 0],\
    gaintable = [ktab, btab, gtab_all], gainfield = ['', '', xcal],\
    interp = ['nearest', 'linear,linear', 'linear'],\
    parang = False)

# XY phase calibration
gaincal(vis = calms, caltable = ptab_amb, selectdata = True,\
    solint = "inf", field = xcal, combine = "scan",\
    refant = ref_ant, gaintype = "XYf+QU", calmode = "ap",\
    smodel = [1, 0.0, 1.0, 0], preavg = 200.0, \
    gaintable = [ktab, btab, gtab_all, kxtab], gainfield = ['', '', xcal, xcal],\
    interp = ['nearest', 'linear,linear', 'linear', 'linear'],\
    parang = False)

#### CBID 1688875155
# Spw = 0 (ich=1819/3639):
#  X-Y phase = 1.14334378498 deg.
#  Fractional Poln: Q = 0.0860520303249, U = -0.0418415255845; P = 0.0956852402694, X = -12.9653252811deg.
#  Net (over baselines) instrumental polarization: 0.000243483938023

#### CBID 1688961378
# Spw = 0 (ich=1819/3639):
#  X-Y phase = 1.12583587571 deg.
#  Fractional Poln: Q = 0.0848468840122, U = -0.0403640158474; P = 0.0939587550388, X = -12.7208647728deg.
#  Net (over baselines) instrumental polarization: 0.001725094285

ptab_xyf = "CASA_Tables/calib.xyf"
S = xyamb(xytab = ptab_amb, qu = qu[2], xyout = ptab_xyf)

#### CBID 1688875155
# Expected QU =  (0.082085790169724032, -0.037996462493649449)
# Spw = 0: Found QU = [ 0.08605203 -0.04184153]
#       ...KEEPING X-Y phase 2.50516551205 deg
# Ambiguity resolved (spw mean): Q= 0.0860520303249 U= -0.0418415255845 (rms= 0.0 0.0 ) P= 0.0956852401694 X= -12.9653258675
# Returning the following Stokes vector: [1.0, 0.086052030324935913, -0.041841525584459305, 0.0]

#### CBID 1688961378
# Expected QU =  (0.079009955873199977, -0.04113488507779748)
# Spw = 0: Found QU = [ 0.08484688 -0.04036402]
#       ...KEEPING X-Y phase 2.57798789117 deg
# Ambiguity resolved (spw mean): Q= 0.0848468840122 U= -0.0403640158474 (rms= 0.0 0.0 ) P= 0.0939587542591 X= -12.7208643829
# Returning the following Stokes vector: [1.0, 0.08484688401222229, -0.040364015847444534, 0.0]


plotcal(caltable = ptab_xy0, xaxis = 'time', yaxis = 'amp', poln = '/', subplot = 111, field = '2',\
    showgui = False, figfile = 'gxgy_vs_parang.png')

plotcal(caltable = ptab_xy0, xaxis = 'chan', yaxis = 'phase', subplot = 111, field = '2',\
    showgui = False, figfile = 'xy_phase.png')

# RE-DERIVE GAIN TABLES
gaincal(vis = calms, caltable = gtab_pol, selectdata = True,\
    solint = "20s", field = bpcal.split(',')[0], combine = "",\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linear'],\
    parang = True)

gaincal(vis = calms, caltable = gtab_pol, selectdata = True,\
    solint = "20s", field = bpcal.split(',')[1], combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linear'],\
    parang = True)

# Gain calibration: gain calibrator
gaincal(vis = calms, caltable = gtab_pol, selectdata = True,\
    solint = "20s", field = gcal, combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap",\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linear'],\
    parang = True)

# Gain calibration: polarisation calibrator
gaincal(vis = calms, caltable = gtab_pol, selectdata = True,\
    solint = "20s", field = xcal, combine = "", append = True,\
    refant = ref_ant, gaintype = "G", calmode = "ap", smodel = S,\
    gaintable = [ktab, btab], gainfield = ['', ''],\
    interp = ['nearest', 'linear,linear'],\
    parang = True)

qu = qufromgain(gtab_pol, paoffset = -90)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### CBID 1688875155
# NB: default band position angle will be offset by -90deg.
# Latitude =  -30.7123828219
# Found as many as 4 fields.

# Fld= 0 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 1.02073413666 Q= -0.0100925079755 U= -0.00151334670885 P= 0.0102053385783 X= -85.7360871
# For field id =  0  there are  1 good spws.
# Spw mean: Fld= 0 Q= -0.0100925079755 U= -0.00151334670885 (rms= 0.0 0.0 ) P= 0.0102053385783 X= -85.7360871

# Fld= 1 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 1.00080507502 Q= -0.0151565788922 U= -0.00681936562871 P= 0.0166200370425 X= -77.8878483446
# For field id =  1  there are  1 good spws.
# Spw mean: Fld= 1 Q= -0.0151565788922 U= -0.00681936562871 (rms= 0.0 0.0 ) P= 0.0166200370425 X= -77.8878483446

# Fld= 2 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 0.993709620832 Q= -0.00274628934727 U= -0.00147377772418 P= 0.00311674926152 X= -75.8900338841
# For field id =  2  there are  1 good spws.
# Spw mean: Fld= 2 Q= -0.00274628934727 U= -0.00147377772418 (rms= 0.0 0.0 ) P= 0.00311674926152 X= -75.8900338841

# Fld= 3 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 0.974118030194 Q= -0.0126171593133 U= 0.0067325336383 P= 0.0143010390646 X= 75.9577787294
# For field id =  3  there are  1 good spws.
# Spw mean: Fld= 3 Q= -0.0126171593133 U= 0.0067325336383 (rms= 0.0 0.0 ) P= 0.0143010390646 X= 75.9577787294

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### CBID 1688961378
# NB: default band position angle will be offset by -90deg.
# Latitude =  -30.7124439492
# Found as many as 4 fields.

# Fld= 0 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 1.0572526425 Q= -0.0270067452389 U= -0.00207672972317 P= 0.0270864743875 X= -87.8013983782
# For field id =  0  there are  1 good spws.
# Spw mean: Fld= 0 Q= -0.0270067452389 U= -0.00207672972317 (rms= 0.0 0.0 ) P= 0.0270864743875 X= -87.8013983782

# Fld= 1 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 0.999414485693 Q= -0.0195205520088 U= -0.0051818811676 P= 0.0201966295001 X= -82.5666356503
# For field id =  1  there are  1 good spws.
# Spw mean: Fld= 1 Q= -0.0195205520088 U= -0.0051818811676 (rms= 0.0 0.0 ) P= 0.0201966295001 X= -82.5666356503

# Fld= 2 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 0.985585719856 Q= -0.00594534031595 U= -0.00652437934417 P= 0.00882692456629 X= -66.1706727008
# For field id =  2  there are  1 good spws.
# Spw mean: Fld= 2 Q= -0.00594534031595 U= -0.00652437934417 (rms= 0.0 0.0 ) P= 0.00882692456629 X= -66.1706727008

# Fld= 3 Spw= 0 Can't discern an ALMA bandname from: none
# Unresolved bandname: default band position angle set to 0.0
# (B=none, PA offset=-90.0deg) Gx/Gy= 0.988230972488 Q= -0.00582672369282 U= 0.00156633718414 P= 0.00603358278031 X= 82.4767416744
# For field id =  3  there are  1 good spws.
# Spw mean: Fld= 3 Q= -0.00582672369282 U= 0.00156633718414 (rms= 0.0 0.0 ) P= 0.00603358278031 X= 82.4767416744

# Leakage calibration
polcal(vis = calms, caltable = dtab, selectdata = True,\
    solint = 'inf', field = xcal, combine = 'scan',\
    refant = '', poltype = 'Dflls', smodel = S, preavg = 200,\
    gaintable = [ktab, btab, gtab_pol, kxtab, ptab_xyf],\
    gainfield = ['', '', xcal, xcal, xcal],\
    interp = ['nearest', 'linear,linear', 'nearest', 'nearest', 'nearest'])

# Generalise leakage solution to cross-hands
Dgen(dtab = dtab, dout = dgen)

# Flux scale correction
fluxscale(vis = calms, caltable = gtab_pol, reference = [fcal],\
    transfer = [bpcal.split(',')[1], gcal, xcal], fluxtable = ftab,\
    append = False, incremental = False)

### APPLY TO CALIBRATORS AND IMAGE 3C138 ETC ...
flagmanager(vis = calms, mode = 'save', versionname = 'PreApplyCal')
#flagmanager(vis = calms, mode = 'restore', versionname = 'PreApplyCal')

applycal(vis  = calms, parang = True, calwt = False, field = '',\
    gaintable = [ktab, btab, kxtab, ptab_xyf, dgen, ftab],\
    gainfield = ['', '', '', '', '', ''],\
    interp    = ['nearest,linear','nearest,linear','nearest,linear','nearest,linear','linear,linear','linear,linear'])


# IMAGE 3C138... OTHER PYTHON FILE
#Fields: 4
#ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#0    T    J1939-6342          19:39:25.030000 -63.42.45.60000 J2000   0         117040
#1    T    J0217+017           02:17:48.950000 +01.44.49.70000 J2000   1         260260
#2    T    J0521+1638          05:21:09.890000 +16.38.22.10000 J2000   2         160160
#3    T    J0408-6545          04:08:20.380000 -65.45.09.10000 J2000   3         113960



### WORK ON TARGET ...
# Split target ...
rawms  = '../RawData/1688961378.MS'
targms = 'MS_Files/1688961378_HCG15.MS'
split(vis = rawms, outputvis = targms, field = "HCG15", datacolumn = 'data', spw = '0:200~3838')


applycal(vis  = targms, parang = True, calwt = False, field = '',\
    gaintable = [ktab, btab, kxtab, ptab_xyf, dgen, ftab],\
    gainfield = ['', '', '', '', '', 'J0217+017'],\
    interp    = ['nearest,linear','nearest,linear','nearest,linear','nearest,linear','linear,linear','linear,linear'])

flagmanager(vis = 'MS_Files/1688961378_HCG15.MS', mode = 'save', versionname = 'ApplyCal')

aoflagger-setup
aoflagger -v -j 32 -strategy meerkat_custom20230417.lua -column CORRECTED_DATA MS_Files/1688961378_HCG15.MS

split(vis = 'MS_Files/1688961378_HCG15.MS', outputvis = 'MS_Files/1688961378_HCG15_corr-avg.MS', datacolumn = 'corrected', width = 8)

