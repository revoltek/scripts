import subprocess
import os
spwlist=['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']
path_to_folder='/export/scratch/AG_deGasperin/alex/VLA/PSZ1_G096.89+24.17/fdg/lustre/aoc/ftp/e2earchive/stage/BI0004/obs_3_flag2_data/imaging_each_spw/'
cell = '2.5arcsec'
size = 1500
targetname="PSZ1G096"

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

for spw_id in spwlist:

	print('Imaging spw',spw_id)

	full_path=path_to_folder+'spw%02i/' % int(spw_id)+'/target.MS'
	
	split_output_vis=path_to_folder+'target_'+'spw%02i' % int(spw_id)+'_averaged.MS'
	#split(vis=full_path, outputvis=split_output_vis, datacolumn='corrected', field=targetname, width=4, timebin='10s')

	"""tclean(vis=split_output_vis, imagename=path_to_folder+'images/'+'spw%02i' % int(spw_id), imsize=size, cell=cell, gridder='wproject', wprojplanes=-1, \
        pblimit=0., deconvolver='mtmfs', nterms=2, weighting='briggs', robust=0.5, niter=1000, \
        pbcor=False)"""
	
	split_output_vis_flag=path_to_folder+'target_'+'spw%02i' % int(spw_id)+'_flagged.MS'
	os.system('cp -r '+full_path+' '+split_output_vis_flag)

	flag_stat(split_output_vis_flag)
	
	subprocess.call('aoflagger -strategy ~/scripts/JVLA_Lband-flag.rfis'+' '+'-column CORRECTED_DATA'+' '+split_output_vis_flag+'>'+'spw%02i' % int(spw_id)+'_target_flagging.log',shell=True)

	flag_stat(split_output_vis_flag)
	
	cal_flagged_averaged_vis=path_to_folder+'spw%02i' % int(spw_id)+'_pre-flagged_averaged_target.MS'
	split(vis=split_output_vis_flag, outputvis=cal_flagged_averaged_vis, datacolumn='corrected', field='PSZ1G096', width=4, timebin='10s')

	flag_stat(cal_flagged_averaged_vis)

	tclean(vis=cal_flagged_averaged_vis, imagename=path_to_folder+'images/flagged/'+'spw%02i' % int(spw_id), imsize=size, cell=cell, gridder='wproject', wprojplanes=-1, \
        pblimit=0., deconvolver='mtmfs', nterms=2, weighting='briggs', robust=0.5, niter=1000, \
        pbcor=False)
	
