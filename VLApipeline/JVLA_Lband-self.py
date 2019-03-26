# to be exectured in CASA
import os, sys
import subprocess

prefix = "C"
active_mss = ['C_target.MS']

# A conf
#cell = '0.3arcsec'
#size = 10000
# B conf
#cell = '1arcsec'
#size = 4000
# C conf
cell = '2.5arcsec'
size = 1500

noise_level='0.017mJy'

loggername=prefix+'self_cal.log'
casalog.setlogfile(loggername)
flaglogfile=prefix+'self_cal_flag_stat.log'


flagging_strategy='~/scripts/Optimised_JVLA_Lband-flag_with_original.rfis'

os.system('rm -r '+prefix+'_img '+prefix+'_cals '+prefix+'_plots')
os.makedirs(prefix+"_img")
os.makedirs(prefix+"_cals")
os.makedirs(prefix+"_plots")

num_ant = 26
num_plots = (num_ant/4)
ref_ant = 'ea20'
next_start_line = 0

def flag_stat(ms):
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
	output_file.write(flag_type) #Write type of flagging 
	output_file.close()

	contval='False'	

	logfile=open(loggername,'r')

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



def clean(active_mss, i):
    print ("Cleaning...")
    tclean(vis=active_mss, imagename=prefix+"_img/self"+i, imsize=size, cell=cell, gridder='wproject', wprojplanes=-1, \
        pblimit=0., deconvolver='mtmfs', nterms=3, weighting='briggs', robust=0.5, niter=10000, usemask='auto-multithresh', \
        pbcor=False, savemodel='modelcolumn')

flag_stat(active_mss)
next_start_line=add_to_flag_log(next_start_line,'Before Self-Cal Flagging')

clean(active_mss, '1')
# PHASE 5m
print ("5m P cal")
for active_ms in active_mss:
    gaintable = prefix+"_cals/"+active_ms+".Gp1"
    gaincal(vis=active_ms, caltable=gaintable, solint="5min", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="p")
    applycal(vis=active_ms, gaintable=[gaintable], calwt=False, flagbackup=False)

    figure_name=active_ms+'.G1p'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[-1,-1,-180,180], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png', showgui=False)

flag_stat(active_mss)
next_start_line=add_to_flag_log(next_start_line,'After 1st Self-Cal Flagging')

clean(active_mss, '2')
# PHASE 1m
print ("1m P cal")
for active_ms in active_mss:
    gaintable = prefix+"_cals/"+active_ms+".Gp2"
    gaincal(vis=active_ms, caltable=gaintable, solint="1min", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="p")
    applycal(vis=active_ms, gaintable=[gaintable], calwt=False, flagbackup=False)

    figure_name=active_ms+'.G2p'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[-1,-1,-180,180], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png', showgui=False)

flag_stat(active_mss)
next_start_line=add_to_flag_log(next_start_line,'After 2nd Self-Cal Flagging')

clean(active_mss, '3')
# PHASE 30s + AMP 10m
print ("30s P + 10m A cal")
for active_ms in active_mss:
    gaintable1 = prefix+"_cals/"+active_ms+".Gp3"
    gaincal(vis=active_ms, caltable=gaintable1, solint="30s", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="p")
    gaintable2 = prefix+"_cals/"+active_ms+".Gap3"
    gaincal(vis=active_ms, caltable=gaintable2, solint="10min", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="ap", gaintable=[gaintable1])
    applycal(vis=active_ms, gaintable=[gaintable1, gaintable2], calwt=False, flagbackup=False)

    figure_name=active_ms+'.G3p'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable1, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[-1,-1,-180,180], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png', showgui=False)

    figure_name=active_ms+'.G3a'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable2, xaxis='time', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png', showgui=False)

flag_stat(active_mss)
next_start_line=add_to_flag_log(next_start_line,'After 3rd Self-Cal Flagging')

clean(active_mss, '4')

#Plots per spw

spwlist=['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']
os.makedirs(prefix+'_plots/plots_per_spw')

for active_ms in active_mss:
	for spw_id in spwlist:

		spw_num='spw%02i' % int(spw_id)
		spw_dir=prefix+'_plots/plots_per_spw/'+spw_num+'/'
		os.makedirs(spw_dir)

		gaintable=prefix+'_cals/'+active_ms+'.Gp1'
		figure_name=spw_dir+active_ms+'.G1p'
        	for i in range(num_plots+1):
            		ant_plot=str(i*4)+'~'+str(i*4+3)
            		plotcal(caltable=gaintable, xaxis='time', yaxis='phase', spw=spw_id, subplot=221, antenna=ant_plot, iteration='antenna', plotrange=[-1,-1,-180,180], figfile=figure_name+str(i+1)+'.png', showgui=False)

		gaintable=prefix+'_cals/'+active_ms+'.Gp2'
		figure_name=spw_dir+active_ms+'.G2p'
        	for i in range(num_plots+1):
            		ant_plot=str(i*4)+'~'+str(i*4+3)
            		plotcal(caltable=gaintable, xaxis='time', yaxis='phase', spw=spw_id, subplot=221, antenna=ant_plot, iteration='antenna', plotrange=[-1,-1,-180,180], figfile=figure_name+str(i+1)+'.png', showgui=False)

		gaintable=prefix+'_cals/'+active_ms+'.Gp3'
		figure_name=spw_dir+active_ms+'.G3p'
        	for i in range(num_plots+1):
            		ant_plot=str(i*4)+'~'+str(i*4+3)
            		plotcal(caltable=gaintable, xaxis='time', yaxis='phase', spw=spw_id, subplot=221, antenna=ant_plot, iteration='antenna', plotrange=[-1,-1,-180,180], figfile=figure_name+str(i+1)+'.png', showgui=False) 

		gaintable=prefix+'_cals/'+active_ms+'.Gap3'
		figure_name=spw_dir+active_ms+'.G3a'
        	for i in range(num_plots+1):
            		ant_plot=str(i*4)+'~'+str(i*4+3)
            		plotcal(caltable=gaintable, xaxis='time', yaxis='amp', spw=spw_id, subplot=221, antenna=ant_plot, iteration='antenna', plotrange=[], figfile=figure_name+str(i+1)+'.png', showgui=False)


#Aoflagger on residuals
'''for active_ms in active_mss:
	uvsub(vis=active_ms,reverse=False) #Corrected - model to create residual dataset in corrected data column
	subprocess.call('aoflagger -strategy '+flagging_strategy+' -column CORRECTED_DATA '+active_ms+' > '+active_ms+'_residual_flagging.log',shell=True) #Flag residuals
	flag_stat(active_ms) #Print amount flagged after flagging residuals
	next_start_line=add_to_flag_log(next_start_line,'After Residuals  Aoflagger Flagging')
	gaintable1=prefix+"_cals/"+active_ms+".Gp3"
	gaintable2=prefix+"_cals/"+active_ms+".Gap3"
	applycal(vis=active_ms, gaintable=[gaintable1, gaintable2], calwt=False, flagbackup=False) #Apply most recent calibration tables to recreate corrected data column
clean(active_mss,'flagged_residuals') #Image again after flagging
'''
tclean(vis=active_mss, imagename=prefix+"_img/Cleaned_to_"+noise_level, imsize=size, cell=cell, gridder='wproject', wprojplanes=-1, \
        pblimit=-1., deconvolver='mtmfs', nterms=3, weighting='briggs', robust=0.5, niter=1000000, threshold=noise_level, usemask='auto-multithresh', \
        pbcor=False) #Image to set threshold (noise_level) 


'''
# then clean with all configurations combined and selfcal with 10min on Gap to aligh amp and phases
gaincal(vis='A1_final.MS', caltable='cals/A1.Gap',solint='10min',calmode='ap')
gaincal(vis='A2_final.MS', caltable='cals/A2.Gap',solint='10min',calmode='ap')
gaincal(vis='B1_final.MS', caltable='cals/B1.Gap',solint='10min',calmode='ap')
gaincal(vis='B2_final.MS', caltable='cals/B2.Gap',solint='10min',calmode='ap')
gaincal(vis='C1_final.MS', caltable='cals/C1.Gap',solint='10min',calmode='ap')
gaincal(vis='C2_final.MS', caltable='cals/C2.Gap',solint='10min',calmode='ap')
'''
