# Run with: nohup casa -c ~/scripts/VLApipeline/JVLA_Lband-cal.py

import os, glob
import numpy as np
import subprocess
pi = np.pi


# TODO: document all entries and move use-related input to the top

basevis = 'B2.MS-orig' # something.MS
basevisnoMS = basevis.replace('.MS','')
no_spws=16 #Number of spws
spwlist = [str(s) for s in range(no_spws)]

uvlimit = "<100000klambda"
uvlimitp= "<100000klambda"

loggername=basevisnoMS+'_cal.log'
casalog.setlogfile(loggername) #Set casa log file as loggername 

flaglogfile=basevisnoMS+'_cal_flag_stat.log' #Set name of file to hold only flag_stat and applycal output
next_start_line=0

targetname = "XDCPJ0044-20"
phasecalname = "J2357-1125"
bandpasscalname = "0542+498=3C147"
polangcalname = "0521+166=3C138"
instpolname = "0542+498=3C147"
basevisnoMS = basevis.replace('.MS','')
# TODO: makes this automatic
calibrator_field_ids='2,3' #Field IDs of calibrator sources for pre-calibration flagging

path_plot = basevisnoMS+'_plots/'
path_data = basevisnoMS+'_data/'
num_ant = 27
num_plots = (num_ant/4)
ref_ant = 'ea20'

# TODO: I updated the path so that it's absolute
flagging_strategy='~/scripts/VLApipeline/Optimised_JVLA_Lband-flag_with_original.rfis'

def flag_stat(ms):
  '''Prints summary of flagging per antenna, array, correlation, field, observation, scan, spw and total flagging.'''
  default('flagdata')
  t = flagdata(vis=ms, mode='summary', field='', scan='', spwchan=False, spwcorr=False, basecnt=False, action='calculate', flagbackup=False, savepars=False)
  log = 'Flag statistics:'
  log += '\nAntenna, '
  for k in sorted(t['antenna']):
      log += k +': %d.2%% - ' % (100.*t['antenna'][k]['flagged']/t['antenna'][k]['total'])
  log += '\nCorrelation, '
  for k, v in t['correlation'].items():
      log += k +': %d.2%% - ' % (100.*v['flagged']/v['total'])
  log += '\nSpw, '
  for k, v in t['spw'].items():
      log += k +': %d.2%% - ' % (100.*v['flagged']/v['total'])
  log += '\nTotal: %d.2%%' % (100.*t['flagged']/t['total'])
  
  print (log.replace(' - \n','\n'))


def add_to_flag_log(start_line,flag_type):
	'''Creates new log file with only flag_stat and applycal output.
	Searches through entire log file currently created for key words, appends all lines from key word to end of task (i.e. until lots of #) into flagging only log file.'''

	output_file=open(flaglogfile,'a+')
	output_file.write('##########################################'+flag_type+'##########################################') #Write type of flagging 
	output_file.write('\n')
	output_file.close()

	logfile=open(loggername,'r')

	contval='False'	

	start_line_num=start_line #Set line that it should start appending from, i.e. don't append anything that has already been run through by a previous call of the function

	i=0
	
	for line in logfile:

		i+=1 #Line number		

		if 'flagdata' and 'mode="summary"' in line and i>start_line_num:
			contval='True'

		elif 'applycal' in line and i>start_line_num:
			contval='True'

		elif contval=='True' and '##########################################' not in line and i>start_line_num:
			contval='True'

		elif contval=='True' and '##########################################' in line:
			output_file=open(flaglogfile,'a+')
			output_file.write(line)
			output_file.close()
			contval='False'	#Append '#'s from last line of flagdata(mode="summary")/applycal but not beyond	
	
		elif contval=='False' and '##########################################' in line: #Don't append '#'s from start of task or from other casa tasks
			contval='False'
		
		else:
			contval='False'

		if contval=='True':
			output_file=open(flaglogfile,'a+')
			output_file.write(line) #Append whilst within flagdata(mode="summary") or applycal
			output_file.close()

	logfile.close()

	return i #Set new line number to start from on next function call

# delete and re-create MS-specific paths
if os.path.exists(path_plot): os.system('rm -r '+path_plot)
os.makedirs(path_plot)
if os.path.exists(path_data): os.system('rm -r '+path_data)
os.makedirs(path_data)
 
allspw_dataset   = path_data+'allspw.MS'
allspw_target    = path_data+'target.MS'

# listobs
if os.path.exists(basevis.replace('MS', 'listobs')): os.system('rm '+basevis.replace('MS', 'listobs'))
listobs( vis=basevis, listfile=basevis.replace('MS', 'listobs') )

print 'Original flagging'
flag_stat(basevis)
next_start_line=add_to_flag_log(next_start_line,'Original Flagging')

# remove setup scan
default("flagdata")
flagdata(vis=basevis, flagbackup=False, mode='manual', scan='1')
print 'Setup scan flagging'
flag_stat(basevis)
next_start_line=add_to_flag_log(next_start_line,'After Setup Scan Flagging')

# FLAG SHADOW
default("flagdata")
flagdata(vis=basevis,mode="shadow",flagbackup=False)
print 'Shadow flagging'
flag_stat(basevis)
next_start_line=add_to_flag_log(next_start_line,'After Shadow Flagging')

#Zero-amplitude data flagging
default(flagdata)
flagdata(vis=basevis,mode='clip',clipzeros=True,flagbackup=False)
print 'Zero amplitude flagging'
flag_stat(basevis)
next_start_line=add_to_flag_log(next_start_line,'After 0 Amplitude Flagging')

#Quack flagging (interval=sampling rate)
default(flagdata)
flagdata(vis=basevis,mode='quack',quackinterval=5.0,quackmode='beg',flagbackup=False)
print 'Quack flagging'
flag_stat(basevis)
next_start_line=add_to_flag_log(next_start_line,'After Quack Flagging')

full_hsmooth_file = path_data+'full_hsmooth.MS'

#apply the hanning smooth because of the RFI
default("hanningsmooth")
hanningsmooth(vis=basevis, datacolumn="data", outputvis=full_hsmooth_file)
flag_stat(full_hsmooth_file)
next_start_line=add_to_flag_log(next_start_line,'After Hanning Smooth')

#aoflagger on calibrators
subprocess.call('aoflagger -strategy '+flagging_strategy+' '+'-column DATA'+' '+'-fields '+calibrator_field_ids+' '+full_hsmooth_file+'>'+basevisnoMS+'_calibrator_flagging.log',shell=True)
print 'Calibrator flagging'
flag_stat(full_hsmooth_file)
next_start_line=add_to_flag_log(next_start_line,'After Calibrator Aoflagger Flagging')

for spw_id in spwlist:

  print ("Working on SPW", spw_id)

  path_plot_spw = path_plot+'/spw%02i/' % int(spw_id)
  os.makedirs(path_plot_spw)
  path_data_spw = path_data+'/spw%02i/' % int(spw_id)
  os.makedirs(path_data_spw)
  
  hsmooth_file = path_data_spw+'hsmooth.MS'
  outname_phasecal = path_data_spw+'phasecal.MS'
  outname_bandpass = path_data_spw+'bandpass.MS'
  outname_polang   = path_data_spw+'polang.MS'
  outname_instpol = path_data_spw+'instpol.MS'
  outname_target   = path_data_spw+'target.MS'
  
  calgceff       = path_data_spw+'cal.gceff'
  calantpos      = path_data_spw+'cal.antpos'
  calgcal1       = path_data_spw+'cal.gcal1' # P
  calbp1         = path_data_spw+'cal.bcal1'
  calgcal2       = path_data_spw+'cal.gcal2' # P
  calk1          = path_data_spw+'cal.kcal1'
  calbp2         = path_data_spw+'cal.bcal2'
  calgcal3       = path_data_spw+'cal.gcal3' # AP
  calkross1      = path_data_spw+'cal.kross1'
  caldf1         = path_data_spw+'cal.dcal1'
  calxf1         = path_data_spw+'cal.xcal1'
  calgcal_loc    = path_data_spw+'cal.gcalloc'

  default('split')
  split(vis=full_hsmooth_file,outputvis=hsmooth_file,spw=spw_id,datacolumn='data')
  
  # GENCAL GAIN CURVE 
  gencal(vis=hsmooth_file,caltable=calgceff,caltype="gceff",spw="",antenna="",pol="",parameter=[])
  
  # ANTPOS
  gencal(vis=hsmooth_file,caltable=calantpos,caltype='antpos',spw="",antenna="",pol="",parameter=[])
  
  if os.path.isdir(calantpos):  
     print ('Antenna position corrections found')
     applycal(vis=hsmooth_file, gaintable=[calgceff,calantpos], calwt=False, flagbackup=False)
     next_start_line=add_to_flag_log(next_start_line,'Antenna Position Correction + Gencal Gain Applycal Spw'+spw_id)	
  else:  
     print ('Antenna position corrections NOT found')
     applycal(vis=hsmooth_file, gaintable=calgceff, calwt=False, flagbackup=False)
     next_start_line=add_to_flag_log(next_start_line,'Gencal Gain Applycal Spw'+spw_id)	

  # BANDPASS
  default("split")
  split(vis=hsmooth_file,outputvis=outname_bandpass,datacolumn="corrected",field=bandpasscalname)
        
  # POLANG
  default("split")
  split(vis=hsmooth_file,outputvis=outname_polang,datacolumn="corrected",field=polangcalname)
 
  #Phasecal, + bandpass for fluxscale
  default("split")
  split(vis=hsmooth_file,outputvis=outname_phasecal,datacolumn="corrected",field=bandpasscalname+','+phasecalname)
  
  #Intrumental polarisation calibrator
  default("split")
  split(vis=hsmooth_file,outputvis=outname_instpol,datacolumn="corrected",field=instpolname)
    
  # TARGET
  default("split")
  split(vis=hsmooth_file,outputvis=outname_target,datacolumn="corrected",field=targetname) 
  
  os.system('rm -rf ' + hsmooth_file) 
  
  setjy(vis=outname_bandpass, field=bandpasscalname, spw="", usescratch=True)
  setjy(vis=outname_phasecal, field=bandpasscalname, spw="", usescratch=True)
  setjy(vis=outname_polang, field=polangcalname, spw="", usescratch=True)
  setjy(vis=outname_instpol, field=instpolname, spw="", usescratch=True)
 
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

  flag_stat(outname_bandpass)
  next_start_line=add_to_flag_log(next_start_line,'Starting Bpass Calibrator Flagging Spw'+spw_id)
  flag_stat(outname_polang)
  next_start_line=add_to_flag_log(next_start_line,'Starting Polang Calibrator Flagging Spw'+spw_id)

  applycal(vis=outname_bandpass, gaintable=[calgcal1,calbp1], interp=['linear','linear,linear'], calwt=False, flagbackup=False)
  next_start_line=add_to_flag_log(next_start_line,'Initial Bandpass On Bandpass Calibrator Applycal Spw'+spw_id)
  applycal(vis=outname_polang, gaintable=[calgcal1,calbp1], interp=['linear','linear,linear'], calwt=False, flagbackup=False)
  next_start_line=add_to_flag_log(next_start_line,'Initial Bandpass On Polarisation Angle Calibrator Applycal Spw'+spw_id)

  flag_stat(outname_bandpass)
  next_start_line=add_to_flag_log(next_start_line,'After Initial Bandpass On Bandpass Calibrator Applycal Spw'+spw_id)
  flag_stat(outname_polang)
  next_start_line=add_to_flag_log(next_start_line,'After Initial Bandpass On Polarisation Angle Calibrator Applycal Spw'+spw_id)

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
  gaincal(vis=outname_bandpass, caltable=calgcal3, spw="0:2~61", calmode='ap', solint="inf", refant=ref_ant, gaintable=gt, interp=['linear','nearest'], parang=True)

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
  polcal(vis=outname_instpol, caltable=caldf1, selectdata=True, uvrange=uvlimit, solint="inf", combine="scan", \
         refant=ref_ant, poltype="Df", gaintable=gt, interp=['linear','linear','linear','linear'], spwmap=[])
         
  figure_name='Df_amp'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=caldf1, xaxis='freq', yaxis='amp', subplot=221, \
              antenna=ant_plot, iteration='antenna', figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='Df_phase'
  for i in range(num_plots+1):
      ant_plot=str(i*4)+'~'+str(i*4+3)
      plotcal(caltable=caldf1, xaxis='freq', yaxis='phase', subplot=221, \
              antenna=ant_plot, iteration='antenna', figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  # Xf
  gt = [calgcal3, calbp2, calk1, calkross1, caldf1]
  polcal(vis=outname_polang, caltable=calxf1, selectdata=True, uvrange=uvlimit, solint="inf", combine="scan", \
         refant=ref_ant, poltype="Xf", gaintable=gt, interp=['linear','linear','linear','linear','nearest'])
      
  figure_name='Xf.png'
  plotcal(caltable=calxf1, xaxis='chan', yaxis='phase', plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name)

  ###################
  # Local cal (amp for elevation dependent problems) 
  gt = [calbp2,calk1,calkross1,caldf1,calxf1]
  gaincal(vis=outname_phasecal, caltable=calgcal_loc, spw="0:2~61", \
          selectdata=True, uvrange=uvlimit, solint="10min", refant=ref_ant, \
          gaintype="G", calmode="ap", gaintable=gt, interp=['linear,linear','nearest','nearest','nearest,linear'], parang=True)
          
  # FLUXSCALE (fitorder not "used" as we have only 1 spw)
  fluxscale(vis=outname_phasecal, caltable=calgcal_loc, fluxtable=calgcal_loc, reference=[bandpasscalname], transfer=phasecalname)
            
  plotcal(caltable=calgcal_loc, xaxis='time', yaxis='amp', figfile=path_plot_spw+'loc_amp.png')
  plotcal(caltable=calgcal_loc, xaxis='time', yaxis='phase', plotrange=[-1,-1,-180,180], figfile=path_plot_spw+'loc_ph.png')

  figure_name='loc_amp' #Plot local amplitude and phase solutions per antenna
  for i in range(num_plots+1):
  	ant_plot=str(i*4)+'~'+str(i*4+3)
  	plotcal(caltable=calgcal_loc, xaxis='time', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna', figfile=path_plot_spw+figure_name+str(i+1)+'.png')

  figure_name='loc_phase'
  for i in range(num_plots+1):
  	ant_plot=str(i*4)+'~'+str(i*4+3)
  	plotcal(caltable=calgcal_loc, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', plotrange=[-1,-1,-180,180], figfile=path_plot_spw+figure_name+str(i+1)+'.png')	
  
  # Apply to target
  gt = [calgcal_loc,calbp2,calk1,calkross1,caldf1,calxf1]
  flag_stat(outname_target)
  applycal(vis=outname_target, gaintable=gt, interp=['linear','linear,linearflag','nearest','nearest','nearest,linearflag','nearest,linearflag'], parang=True, calwt=False, flagbackup=False)
  next_start_line=add_to_flag_log(next_start_line,'Target Applycal Spw'+spw_id)
  flag_stat(outname_target)
  next_start_line=add_to_flag_log(next_start_line,'After Target Applycal Spw'+spw_id)

#put together all the spw in one dataset (that is still divided into no_spw spws)
all_spwlist = sorted(glob.glob(basevisnoMS+'_data/spw??/target.MS'))
concat(vis=all_spwlist, concatvis=allspw_dataset)
flag_stat(allspw_dataset)
next_start_line=add_to_flag_log(next_start_line,'After SPW concatenation')

#run AOFlagger on target field corrected data before self cal
subprocess.call('aoflagger -strategy '+flagging_strategy+' '+'-column CORRECTED_DATA'+' '+allspw_dataset+'>'+basevisnoMS+'_target_flagging.log',shell=True)
flag_stat(allspw_dataset)
next_start_line=add_to_flag_log(next_start_line,'After Target Aoflagger Flagging')

#split and averaging; we ignore the initial and final channel because they are bad (total of 16 channels)
split(vis=allspw_dataset, outputvis=allspw_target, datacolumn='corrected', field=targetname, width=4, timebin='10s')
flag_stat(allspw_target)
next_start_line=add_to_flag_log(next_start_line,'After Averaging')

# the allspw_target can now be used for the selfcal part
