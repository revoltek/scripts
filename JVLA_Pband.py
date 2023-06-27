# calibration P band virgo data in casa 6.5.5
from casatasks.private import tec_maps

msbkp = 'bkp/13B-091.Cconf.ms'
ms=ms.replace('bkp/','')
bpcal = '1331+305=3C286'
bpcalscan = 3
plcal = 'OQ208'

# other variables
caldir = ms+'-cal'
imgdir = ms+'-img'

# hanningsmooth
hanningsmooth(vis=bkpms,outputvis=ms,datacolumn='data')

# generate calibration tables
gencal(vis=ms,caltable=caldir+'/antpos.cal',caltype='antpos')
tec_maps.create(vis=ms,doplot=True,imname=caldir+'/iono')
gencal(vis=ms,caltable=caldir+'/tec.cal',caltype='tecim',infile=caldir+'/iono.IGS_TEC.im')
gencal(vis=ms,caltype='rq',caltable=caldir+'/rq.cal')
plotms(vis=ms,xaxis='freq',yaxis='amp',antenna='ea01',correlation='XX,YY', field=bpcal,\
           plotrange=[0.2,0.5,0.0,100.0], coloraxis='spw',xlabel='Frequency',ylabel='Amplitude',iteraxis='baseline',\
           plotfile=caldir+'/amp-init.png')

flagdata(vis=ms,mode='manual',scan='1')
plotms(vis=ms,xaxis='freq',yaxis='amp',antenna='ea01',correlation='XX,YY', field=bpcal,\
           plotrange=[0.2,0.5,0.0,100.0], coloraxis='spw',xlabel='Frequency',ylabel='Amplitude',iteraxis='baseline',\
           plotfile=caldir+'/amp-HS.png')

# setjy
clearcal(vis=ms)
# Fix polarisation model
execfile('/hshome/baq1889/scripts/JVLA_fixflux.py')
I, alpha, _ = get_I('3c286.dat')
PF = get_PF('3c286.dat')
PA = get_PA('3c286.dat')

setjy(vis=ms, field=bpcal, scalebychan=True, standard="manual", fluxdensity=[I,0,0,0], spix=alpha, reffreq='3.0GHz', polindex=PF, polangle=PA, interpolation="nearest", usescratch=True)
#setjy(vis=ms, standard='Scaife-Heald 2012', field=bpcal, usescratch=True)

# flag 1
def printflag(summary_old, summary_new):
    axis = 'scan'
    for value, stats in summary_new[ axis ].iteritems():
        old_stats = summary_old[ axis ][ value ]
        print ('%s %s: %5.1f percent flagged additionally' % ( axis, value, 100. * ( stats[ 'flagged' ] - old_stats[ 'flagged' ] ) / stats[ 'total' ] ))

summary_1 = flagdata(vis=ms, mode='summary')
flagdata(vis=ms, field='*', mode='tfcrop', datacolumn='data', timecutoff=4., freqcutoff=3., maxnpieces=5,\
  action='apply', display='report', flagbackup=True, combinescans=True, ntime='3600s', correlation='ABS_XY,ABS_YX')
flagdata(vis=ms, field='*', mode='tfcrop', datacolumn='data', timecutoff=3., freqcutoff=3., maxnpieces=2,\
  action='apply', display='report', flagbackup=False, combinescans=True, ntime='3600s', correlation='ABS_XX,ABS_YY')
flagdata(vis=ms, mode='extend')
summary_2 = flagdata(vis=ms, mode='summary')
printflag(summary_1, summary_2)

# quick BP claib
gaincal(vis=ms, caltable=caldir+'/pband.G0', gaintype='G', calmode='p', solint='int', field=bpcal,refant='ea09',\
  gaintable=[caldir+'/antpos.cal',caldir+'/rq.cal',caldir+'/tec.cal'])
gaincal(vis=ms, caltable=caldir+'/pband.K0', gaintype='K', solint='inf', field=bpcal,refant='ea09',\
  gaintable=[caldir+'/antpos.cal',caldir+'/rq.cal',caldir+'/tec.cal',caldir+'/pband.G0'])
bandpass(vis=ms, caltable=caldir+'/pband.B0', solint='inf', field=bpcal,refant='ea09', minsnr=2.0,\
  gaintable=[caldir+'/antpos.cal',caldir+'/rq.cal',caldir+'/tec.cal',caldir+'/pband.G0',caldir+'/pband.K0'])
applycal(vis=ms, field=bpcal, applymode='calflagstrict',\
  gaintable=[caldir+'/antpos.cal',caldir+'/rq.cal',caldir+'/tec.cal',caldir+'/pband.G0',caldir+'/pband.K0',caldir+'/pband.B0'])

# flag 2
flagdata(vis=ms, field=bpcal, mode='rflag', datacolumn='corrected', timedevscale=4., freqdevscale=3.,\
  action='apply',  flagbackup=True, combinescans=True, ntime='3600s', maxnpieces=5)
flagdata(vis=ms, field=bpcal, mode='rflag', datacolumn='corrected', timedevscale=4., freqdevscale=3.,\
  action='apply',  flagbackup=True, combinescans=True, ntime='3600s', maxnpieces=5)
summary_3 = flagdata(vis=ms , mode='summary')
printflag(summary_2, summary_3)

# delay
gaincal(vis=ms, caltable=caldir+'/pband.G1pre', field=bpcal, solint='int', refant='ea09', spw='5:60~63', gaintype='G',calmode='p', minsnr=5, \
    gaintable=[caldir+'/antpos.cal',caldir+'/rq.cal',caldir+'/tec.cal'],parang=True)
plotms(vis=caldir+'/pband.G1pre', coloraxis='antenna1', xaxis='time', yaxis='phase', plotfile=caldir+'/pband.G1pre.png', overwrite=True)
gaincal(vis=ms, caltable=caldir+'/pband.K1', field=bpcal, solint='inf', refant='ea09', gaintype='K',\
  gaintable=[caldir+'/antpos.cal',caldir+'/rq.cal',caldir+'/tec.cal',caldir+'/pband.G1pre'],spwmap=[[],[],[],[5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]],parang=True)
plotms(vis=caldir+'/pband.K1', coloraxis='antenna1', xaxis='spw',  plotfile=caldir+'/pband.K1.png', overwrite=True)

# bp
bandpass( vis=ms, caltable=caldir+'/pband.B1', field=bpcal, solint='inf', refant='ea09', minsnr=3.0,\
                gaintable=[caldir+'/antpos.cal',caldir+'/tec.cal',caldir+'/rq.cal',caldir+'/pband.K1'], \
                interp=['','','','nearest,nearestflag','nearest,nearest'], parang=True)
plotms(vis=caldir+'/pband.B1', coloraxis='antenna1', yaxis='amp', iteraxis='spw', plotfile=caldir+'/pband.B1a.png', overwrite=True)
plotms(vis=caldir+'/pband.B1', coloraxis='antenna1', yaxis='phase', iteraxis='spw', plotfile=caldir+'/pband.B1p.png', overwrite=True)

# gain
gaincal( vis=ms, caltable=caldir+'/pband.G1', field=bpcal, solint = 'int',\
             refant = 'ea09', minsnr = 3.0, gaintype = 'G', calmode = 'ap',\
             gaintable=[caldir+'/antpos.cal',caldir+'/tec.cal',caldir+'/rq.cal',caldir+'/pband.K1', caldir+'/pband.B1'], parang=True,\
             interp = ['','','','nearest,nearestflag', 'nearest,nearestflag' ])
plotms(vis=caldir+'/pband.G1', coloraxis='antenna1', yaxis='amp', iteraxis='spw', plotfile=caldir+'/pband.G1a.png', overwrite=True)
plotms(vis=caldir+'/pband.G1', coloraxis='antenna1', yaxis='phase', iteraxis='spw', plotfile=caldir+'/pband.G1p.png', overwrite=True)

smoothcal(vis=ms, tablein=caldir+'/pband.G1', caltable=caldir+'/pband.Gs1', smoothtype='median', smoothtime = 60.*60.)
plotms(vis=caldir+'/pband.Gs1', coloraxis='antenna1', yaxis='amp', plotfile=caldir+'/pband.G1a.png', overwrite=True)
plotms(vis=caldir+'/pband.Gs1', coloraxis='antenna1', yaxis='phase', plotfile=caldir+'/pband.G1p.png', overwrite=True)

# Polcal
# https://library.nrao.edu/public/memos/evla/EVLAM_201.pdf (page 23)
# GAINCAL‘KCROSS’smodel=provided→POLCAL‘Xf’smodel=provided→ POLCAL‘Dflls’smodel=providedrefant=none(→iterate)

# X-Y cross delay calib
gaincal(vis=ms, caltable=caldir+'/pband.Kc1', field=bpcal, scan='', solint='inf', refant='ea09', minsnr=3.0, gaintype='KCROSS',\
  gaintable=[ caldir+'/antpos.cal',caldir+'/tec.cal',caldir+'/rq.cal',caldir+'/pband.K1', caldir+'/pband.B1', caldir+'/pband.G1'], parang=True,\
  interp=['','','','nearest,nearestflag', 'nearest,nearestflag', 'nearest,nearestflag'] )
plotms(vis=caldir+'/pband.Kc1', coloraxis='spw', plotfile=caldir+'/pband.Kc1.png', overwrite=True)

# Leakage
polcal(vis=ms, caltable=caldir+'/pband.Dflls1', field=bpcal, solint='inf', refant='ea09', minsnr=3.0, poltype = 'Dflls',\
  gaintable=[ caldir+'/antpos.cal',caldir+'/tec.cal',caldir+'/rq.cal',caldir+'/pband.K1', caldir+'/pband.B1', caldir+'/pband.G1', caldir+'/pband.Kc1'], \
  interp=['','','','nearest,nearestflag', 'nearest,nearestflag', 'nearest,nearestflag', 'nearest,nearestflag'] )
plotms(vis=caldir+'/pband.Dflls1', coloraxis='corr', xaxis='Frequency', yaxis='Amp', iteraxis='antenna', plotfile=caldir+'/pband.Dflls1.png', overwrite=True)
plotms(vis=caldir+'/pband.Dflls1', coloraxis='corr', xaxis='Frequency', yaxis='Phase', iteraxis='antenna', plotfile=caldir+'/pband.Dflls1.png', overwrite=True)

# angle
polcal(vis=ms, caltable=caldir+'/pband.Xf1', field=bpcal, solint='inf,2MHz', refant='ea09', minsnr=3.0, poltype = 'Xf',\
  gaintable=[ caldir+'/antpos.cal',caldir+'/tec.cal',caldir+'/rq.cal',caldir+'/pband.K1', caldir+'/pband.B1', caldir+'/pband.G1', caldir+'/pband.Kc1', caldir+'/pband.Dflls1'], \
  interp=['','','','nearest,nearestflag', 'nearest,nearestflag', 'nearest,nearestflag', 'nearest,nearestflag', 'nearest,nearestflag'] )
plotms(vis=caldir+'/pband.Xf1', coloraxis='spw', xaxis='Frequency', yaxis='Phase', plotfile=caldir+'/pband.Xf1.png', overwrite=True)


###
# apply
applycal(vis=ms, applymode='calflagstrict', flagbackup=True, \
  gaintable=[ caldir+'/antpos.cal',caldir+'/tec.cal',caldir+'/rq.cal',caldir+'/pband.K1', caldir+'/pband.B1', caldir+'/pband.Gs1', caldir+'/pband.Kc1'], parang=True, \
  interp=['','','','nearest,nearestflag', 'nearest,nearestflag', 'nearest,nearestflag', 'nearest,nearestflag'] )
