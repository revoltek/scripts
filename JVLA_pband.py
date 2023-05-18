# calibration P band virgo data in casa 6.5.5
from casatasks.private import tec_maps

ms = 'bkp/13B-091.Cconf.ms'
bpcal = '1331+305=3C286'

# other variables
caldir = ms+'-cal'

# generate calibration tables
gencal(vis=ms,caltable=caldir+'/antpos.cal',caltype='antpos')
tec_maps.create(vis=ms,doplot=True,imname=caldir+'/iono')
gencal(vis=ms,caltable=caldir+'/tec.cal',caltype='tecim',infile=caldir+'/iono.IGS_TEC.im')
gencal(vis=ms,caltype='rq',caltable=caldir+'/rq.cal')
plotms(vis=ms,xaxis='freq',yaxis='amp',antenna='ea01',correlation='XX,YY', field=bpcal,\
           plotrange=[0.2,0.5,0.0,100.0], coloraxis='spw',xlabel='Frequency',ylabel='Amplitude',iteraxis='baseline',\
           plotfile=caldir+'/amp-init.png')

# hanningsmooth
hanningsmooth(vis=ms,outputvis=ms.replace('bkp/',''),datacolumn='data')
ms=ms.replace('bkp/','')
flagdata(vis=ms,mode='manual',scan='1')
plotms(vis=ms,xaxis='freq',yaxis='amp',antenna='ea01',correlation='XX,YY', field=bpcal,\
           plotrange=[0.2,0.5,0.0,100.0], coloraxis='spw',xlabel='Frequency',ylabel='Amplitude',iteraxis='baseline',\
           plotfile=caldir+'/amp-HS.png')

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

# setjy
clearcal(vis=ms)
setjy(vis=ms, standard='Scaife-Heald 2012', field=bpcal)
# delay
gaincal(vis=ms, caltable=caldir+'/pband.K1', field=bpcal, solint='inf', refant='ea09', gaintype='K',\
  gaintable=[caldir+'/antpos.cal',caldir+'/rq.cal',caldir+'/tec.cal'],parang=True)
plotcal(caltable=caldir+'/pband.K1', xaxis='antenna', yaxis='delay', markersize=2.0,\
            plotrange=[0,30,-50.,50.],figfile=caldir+'/delay.png')
# bp
bandpass( vis=ms, caltable=caldir+'/pband.B1', field=bpcal, solint='inf', refant='ea09', minsnr=3.0,\
                gaintable=[caldir+'/antpos.cal',caldir+'/tec.cal',caldir+'/rq.cal',caldir+'/pband.K1'], parang=True,\
                interp=['','','','nearest,nearestflag'])
plotcal(caltable=caldir+'/pband.B1', xaxis='freq', yaxis='amp', iteration='antenna', markersize=2.0)
plotcal(caltable=caldir+'/pband.B1', xaxis='freq', yaxis='phase', iteration='antenna', plotrange=[ 200.,500.,-180.,180. ], markersize = 2.0)

# gain
gaincal( vis=ms, caltable=caldir+'/pband.G1', field=bpcal, solint = 'int',\
             refant = 'ea09', minsnr = 3.0, gaintype = 'G', calmode = 'ap',\
             gaintable=[caldir+'/antpos.cal',caldir+'/tec.cal',caldir+'/rq.cal',caldir+'/pband.K1', caldir+'/pband.B1'], parang=True,\
             interp = ['','','','nearest,nearestflag', 'nearest,nearestflag' ])


