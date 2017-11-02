import os, glob
import numpy as np
import subprocess
pi = np.pi

# fix for rficonsole LD_LIBRARY_PATH (cannot handle libgfortran provided by casapy)
#for the AOFlagger in the Reinout environment
#my_env = os.environ.copy()
#my_env["LD_LIBRARY_PATH"]="/lib:/usr/lib:/home/rvweeren/software/lofar/opt/LofIm/lib64:/home/rvweeren/software/root/lib:/home/rvweeren/software/heasoft-6.16/x86_64-unknown-linux-gnu-libc2.19-0/lib"
#GDG: no my_env in the subprocess for AOFlagger

# parang
# setjy modimage
# interp nearest
# timeragne
# check that modimage is right band (twice in scripts)
# no tfcrop flaggin, problem for calibrator scans with more than 50% flagged due to slews for example

# Nov 2015
# fixed tfcrop problem by unflag and flag action with tb.open
# linearflag in applycal last step!!
# changed to Taylor Butler 2013
# fixed aoflagger call with LD_LIBRARY_PATH and new rfis

basevis = 'C1.MS' # something.MS
spwlist = ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']

targetname = "XDCPJ0044-20"
phasecalname = "J2357-1125"
bandpasscalname = "0542+498=3C147"
polangcalname = "0521+166=3C138"

allspw_dataset   = 'C1_allspw.MS'
allspw_target    = 'C1_target.MS'

path_plot='plots/'
path_cal='cals/'
num_ant=27
num_plots=(num_ant/4)
ref_ant = 'ea20'

# listobs
if os.path.exists(basevis.replace('MS', 'listobs')): os.system('rm '+basevis.replace('MS', 'listobs'))
listobs( vis=basevis, listfile=basevis.replace('MS', 'listobs') )

for spw_id in spwlist:

  path_plot_spw = path_plot+'/spw%02i/' % int(spw_id)
  if not os.path.exists(path_plot_spw):
    os.makedirs(path_plot_spw)

  uvlimit = "<1000klambda"
  uvlimitp= "<1000klambda"
  
  hsmooth_file     = basevis+'_spw%02i.ms' % int(spw_id)
  outname_phasecal = basevis+'_spw%02i.phasecal.ms' % int(spw_id)
  outname_bandpass = basevis+'_spw%02i.bandpass.ms' % int(spw_id)
  outname_polang    = basevis+'_spw%02i.polang.ms' % int(spw_id)
  outname_target     = basevis+'_spw%02i.target.ms' % int(spw_id) 
  
  calgceff       = path_cal+'cal.gceff.spw.' + str(spw_id)
  calantpos      = path_cal+'cal.antpos.spw.'+ str(spw_id)
  calgcal1       = path_cal+'cal.gcal1.spw.' + str(spw_id) # P
  calbp1         = path_cal+'cal.bcal1.spw.' + str(spw_id)
  calgcal2       = path_cal+'cal.gcal2.spw.' + str(spw_id) # P
  calk1          = path_cal+'cal.kcal1.spw.' + str(spw_id)
  calbp2         = path_cal+'cal.bcal2.spw.' + str(spw_id)
  calgcal3       = path_cal+'cal.gcal3.spw.' + str(spw_id) # AP
  calkross1      = path_cal+'cal.kross1.spw.'+ str(spw_id)
  caldf1         = path_cal+'cal.dcal1.spw.' + str(spw_id)
  calxf1         = path_cal+'cal.xcal1.spw.' + str(spw_id)
  calgcal_loc    = path_cal+'cal.gcalloc.spw.' + str(spw_id) 

  os.system('rm -rf ' + hsmooth_file) 
  os.system('rm -rf ' + outname_phasecal)
  os.system('rm -rf ' + outname_bandpass)
  os.system('rm -rf ' + outname_polang)
  os.system('rm -rf ' + outname_target)
  os.system('rm -rf ' + calgceff)
  os.system('rm -rf '+calbp1+' '+calbp2+' '+calgcal1+' '+calgcal2+' '+caldf1+' '+\
             calxf1+' '+calantpos+' '+calk1+' '+calkross1+' '+calgcal3)
  os.system('rm -rf ' + allspw_dataset)
  os.system('rm -rf ' + allspw_target)
  
  #apply the hanning smooth because of the RFI
  default("hanningsmooth")
  hanningsmooth(vis=basevis, datacolumn="data", outputvis=hsmooth_file, spw=spw_id)

  # remove setup scan
  default("flagdata")
  flagdata(vis=hsmooth_file, flagbackup=False, mode='manual', scan='1')

  # FLAG SHADOW
  default("flagdata")
  flagdata(vis=hsmooth_file,mode="shadow",flagbackup=False)

  # first half of spw3 is shit
  #if spw_id == '04':    
  #  flagdata(vis=hsmooth_file, mode="manual", spw="0:0~32", flagbackup=False)

  # TFCROP FLAG
  #flagdata(vis=hsmooth_file,mode="tfcrop",autocorr=False,inpfile="",reason="any",spw="",field="",antenna="",   
  #         uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",             	   
  #         clipminmax=[],datacolumn="DATA",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0,
  #         quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  
  #         ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        
  #         maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         
  #         freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, 
  #         growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      
  #         maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",
  #         flagbackup=False,savepars=False,cmdreason="",outfile="")
  
  # GENCAL GAIN CURVE 
  gencal(vis=hsmooth_file,caltable=calgceff,caltype="gceff",spw="",antenna="",pol="",parameter=[])
  
  # ANTPOS
  gencal(vis=hsmooth_file,caltable=calantpos,caltype='antpos',spw="",antenna="",pol="",parameter=[])
  
  if os.path.isdir(calantpos):  
     print 'Antenna position corrections found'
     applycal(vis=hsmooth_file, gaintable=[calgceff,calantpos], calwt=False, flagbackup=False)	
  else:  
     print 'Antenna position corrections NOT found'
     applycal(vis=hsmooth_file, gaintable=calgceff, calwt=False, flagbackup=False)	

  # BANDPASS
  default("split")
  split(vis=hsmooth_file,outputvis=outname_bandpass,datacolumn="corrected",field=bandpasscalname)
        
  # POLANG
  default("split")
  split(vis=hsmooth_file,outputvis=outname_polang,datacolumn="corrected",field=polangcalname)
 
  #Phasecal, + bandpass for fluxscale
  default("split")
  split(vis=hsmooth_file,outputvis=outname_phasecal,datacolumn="corrected",field=bandpasscalname+','+phasecalname)
    
  # TARGET
  default("split")
  split(vis=hsmooth_file,outputvis=outname_target,datacolumn="corrected",field=targetname) 
  
  os.system('rm -rf ' + hsmooth_file)   
  
  setjy(vis=outname_bandpass, field=bandpasscalname, spw="", usescratch=True)
  setjy(vis=outname_phasecal, field=bandpasscalname, spw="", usescratch=True)
  setjy(vis=outname_polang, field=polangcalname, spw="", usescratch=True)
 
  # Fix polarisation
  tb.open(outname_polang, nomodify=False)
  data = tb.getcol('MODEL_DATA')
  if '3C286' in polangcalname:
    RLphase = 33.0*2.0
    polfrac = 0.095 # at 1.45GHz
  elif '3C138' in polangcalname:
    RLphase = -11.*2.0
    polfrac = 0.075

  StokesI = (0.5*(np.copy(data[0,:,:]) + np.copy(data[3,:,:])))
  Q = StokesI*polfrac*np.cos(RLphase*pi/180.0)
  U = StokesI*polfrac*np.sin(RLphase*pi/180.0)
  RL = (Q + (complex(0,1)*U))
  LR = (Q - (complex(0,1)*U))
  
  data[1,:,:] = RL
  data[2,:,:] = LR
  tb.putcol('MODEL_DATA', data)
  tb.flush()
  tb.close()

  ##################
  # INITPH 1
  gaincal(vis=outname_bandpass, caltable=calgcal1, spw="0:25~35", calmode='p', solint="int", refant=ref_ant, parang=True)
  
  figure_name='initph'
  for i in range(num_plots+1):
	  ant_plot=str(i*4)+'~'+str(i*4+3)
	  plotcal(caltable=calgcal1, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna',
				plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  # BANDPASS 1
  bandpass(vis=outname_bandpass, caltable=calbp1, solint="inf", combine="scan,field", refant=ref_ant, gaintable=[calgcal1], parang=True)
	   
  figure_name='initbpassA'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=calbp1, xaxis='freq', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna', \
				plotrange=[], figfile=path_plot_spw+figure_name+str(i+1)+'.png')
    
  figure_name='initbpassP'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=calbp1, xaxis='freq', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
				plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')
              
  # apply above solution to the individual sources for flagging purposes (rficonsole)
  applycal(vis=outname_bandpass, gaintable=[calgcal1,calbp1], interp=['linear','linear,linear'], flagbackup=False)
  applycal(vis=outname_polang, gaintable=[calgcal1,calbp1], interp=['linear','linear,linear'], flagbackup=False)

  #################
  # INITPH 2
  gaincal(vis=outname_bandpass, caltable=calgcal2, spw="0:25~35", calmode='p', solint="int", refant=ref_ant, parang=True)
	  
  figure_name='initph2'
  for i in range(num_plots+1):
	  ant_plot=str(i*4)+'~'+str(i*4+3)
	  plotcal(caltable=calgcal2, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
				plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  # DELAY
  gaincal(vis=outname_bandpass, caltable=calk1, spw="0:2~61", selectdata=True, uvrange=uvlimit,\
          solint="inf", combine="scan,field", refant=ref_ant, gaintype="K", gaintable=[calgcal2], parang=True)
         
  figure_name='delay'
  plotcal(caltable=calk1, xaxis='antenna', yaxis='delay', figfile=path_plot_spw+figure_name+'.png')

  # BANDPASS 2 
  gt = [calgcal2,calk1]
  bandpass(vis=outname_bandpass, caltable=calbp2, selectdata=True, uvrange=uvlimit, solint="inf", combine="scan,field", \
           refant=ref_ant, gaintable=gt, parang=True)

  figure_name='bpassA'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=calbp2, xaxis='freq', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna', 
				plotrange=[], figfile=path_plot_spw+figure_name+str(i+1)+'.png')
  figure_name='bpassP'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=calbp2, xaxis='freq', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna',
				plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  # GAIN CAL
  gt = [calk1,calbp2]
  gaincal(vis=outname_bandpass, caltable=calgcal3, calmode='ap', solint="inf", refant=ref_ant, gaintable=gt, interp=['linear','nearest'], parang=True)

  figure_name='finalamp'
  plotcal(caltable=calgcal3, xaxis='antenna', yaxis='amp', figfile=path_plot_spw+figure_name+'.png')
  figure_name='finalph'
  plotcal(caltable=calgcal3, xaxis='antenna', yaxis='phase', plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+'.png')

  ##################
  # KCROSS: residual differences on the reference antenna
  gt = [calgcal3,calbp2,calk1]
  gaincal(vis=outname_polang, caltable=calkross1, spw="0:2~61", selectdata=True, uvrange=uvlimitp, \
          solint="inf", combine="scan", refant=ref_ant, gaintype="KCROSS",append=False, gaintable=gt, interp=['linear','linear','linear,linear'], parang=True)
          
  figure_name='KCROSS'
  plotcal(caltable=calkross1, xaxis='antenna', yaxis='delay',figfile=path_plot_spw+figure_name+'.png')
  
  # Df on unpol cal (if using phasecal solve with Df+QU)
  gt = [calgcal3, calbp2, calk1, calkross1]
  polcal(vis=outname_bandpass, caltable=caldf1, selectdata=True, uvrange=uvlimit, solint="inf", combine="scan", \
         refant=ref_ant, poltype="Df", gaintable=gt, interp=['linear','linear','linear','linear'], spwmap=[])
         
  figure_name='Df_real'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=caldf1, xaxis='freq', yaxis='real', subplot=221, \
              antenna=ant_plot, iteration='antenna', figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='Df_imag'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=caldf1, xaxis='freq', yaxis='imag', subplot=221, \
              antenna=ant_plot, iteration='antenna', figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  # Xf
  gt = [calgcal3, calbp2, calk1, calkross1, caldf1]
  polcal(vis=outname_bandpass, caltable=calxf1, selectdata=True, uvrange=uvlimit, solint="inf", combine="scan", \
         refant=ref_ant, poltype="Xf", gaintable=gt, interp=['linear','linear','linear','linear','nearest'])
      
  figure_name='Xf.png'
  plotcal(caltable=calxf1, xaxis='chan', yaxis='phase', figfile=path_plot_spw+figure_name)

  ###################
  # Local cal (amp for elevation dependent problems) 
  gt = [calbp2,calk1,calkross1,caldf1,calxf1]
  gaincal(vis=outname_phasecal, caltable=calgcal_loc, spw="0:2~61", \
          selectdata=True, uvrange=uvlimit, solint="10min", refant=ref_ant, \
          gaintype="G", calmode="ap", gaintable=gt, interp=['linear,linear','nearest','nearest','nearest,linear'], parang=True)
          
  #figure_name='loc_amp-noscale'
  #for i in range(num_plots+1):
  #    ant_plot=str(i*4)+'~'+str(i*4+3)
  #    plotcal(caltable=calgcal_loc, xaxis='time', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna', \
  #            figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  # FLUXSCALE (fitorder not "used" as we have only 1 spw)
  fluxscale(vis=outname_phasecal, caltable=calgcal_loc, fluxtable=calgcal_loc, reference=[bandpasscalname], transfer=phasecalname)
            
  figure_name='loc_amp'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=calgcal_loc, xaxis='time', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna', \
              figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='loc_phase'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=calgcal_loc, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', plotrange=[-1,-1,-180,180], \
              figfile=path_plot_spw+figure_name+str(i+1)+'.png')

#  # APPLY
#  gt = [calgcal_loc,calbp2,calk1,calkross1,caldf1,calxf1]
#  applycal(vis=outname_target, gaintable=gt, interp=['linear','nearest','nearest','nearest','nearest','nearest'], parang=True, calwt=False, flagbackup=False)
#
#  # TFCROP FLAG
#  flagdata(vis=outname_target,mode="tfcrop",autocorr=False,inpfile="",reason="any",spw="",field="",antenna="",   \
#           uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",           \
#  	   clipminmax=[],datacolumn="corrected",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0,\
#  	   quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  \
#  	   ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        \
#  	   maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         \
#  	   freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, \
#  	   growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      \
#  	   maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",\
#  	   flagbackup=False,savepars=False,cmdreason="",outfile="")
#  	   
#  #RFLAG FLAG   
#  flagdata(vis=outname_target,mode="rflag",autocorr=False,inpfile="",reason="any",spw="",field="",antenna="",   \
#	uvrange="",timerange="",correlation="",scan="",intent="",array="",observation="",feed="",  \
#	clipminmax=[],datacolumn="corrected",clipoutside=True,channelavg=False,clipzeros=False,quackinterval=1.0,\
#  	quackmode="beg",quackincrement=False,tolerance=0.0,addantenna="",lowerlimit=0.0,upperlimit=90.0,  \
#  	ntime="scan",combinescans=False,timecutoff=4.0,freqcutoff=3.0,timefit="line",freqfit="poly",        \
#  	maxnpieces=7,flagdimension="freqtime",usewindowstats="none",halfwin=1,winsize=3,timedev="",         \
#  	freqdev="",timedevscale=5.0,freqdevscale=5.0,spectralmax=1000000.0,spectralmin=0.0,extendpols=True, \
#  	growtime=50.0,growfreq=50.0,growaround=False,flagneartime=False,flagnearfreq=False,minrel=0.0,      \
#  	maxrel=1.0,minabs=0,maxabs=-1,spwchan=False,spwcorr=False,basecnt=False,action="apply",display="",\
#  	flagbackup=False,savepars=False,cmdreason="",outfile="")
  	
  
  # Apply again, but now with the calibration flagging on
  gt = [calgcal_loc,calbp2,calk1,calkross1,caldf1,calxf1]
  applycal(vis=outname_target, gaintable=gt, interp=['linear','linear,linearflag','nearest','nearest','nearest,linearflag','nearest,linearflag'], parang=True, calwt=False, flagbackup=False)

#put together all the spw in one dataset (that is still divided into 16 spw)
all_spwlist = sorted(glob.glob('spw??.target.ms'))
concat(vis=all_spwlist, concatvis=allspw_dataset)
os.system('rm -rf spw??.target.ms')
	
#split and avareging; we ignore the initial and final channel because they are bad (total of 16 channels)
split(vis=allspw_dataset, outputvis=allspw_target, datacolumn='corrected', field=targetname, spw='*:7~54', width=4, timebin='10s')
os.system('rm -rf ' + allspw_dataset)

#run AOFlagger to better remove the RFI
subprocess.call('aoflagger -strategy ../flagonavg_rlandlr.rfis -column DATA '+allspw_target, shell=True)

#further flagging after rl and lr correlation inspection -- on spw
#plotms(vis=allspw_target, xaxis='uvdist',yaxis='amp',correlation='rl,lr',coloraxis='antenna1')

#badspw='1,2,3,8,9'    
#flagdata(vis=allspw_target,mode="manual",spw=badspw,flagbackup=False)
