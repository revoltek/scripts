# to be exectured in CASA
import os, sys

prefix = "A"
active_mss = ['A1_target.MS', 'A2_target.MS']

# A conf
cell = '0.3arcsec'
size = 10000
# B conf
#cell = '1arcsec'
#size = 4000
# C conf
#cell = '2.5arcsec'
#size = 1500

os.system('rm -r '+prefix+'_img '+prefix+'_cals '+prefix+'_plots')
os.makedirs(prefix+"_img")
os.makedirs(prefix+"_cals")
os.makedirs(prefix+"_plots")

num_ant = 27
num_plots = (num_ant/4)
ref_ant = 'ea20'

def clean(active_mss, i):
    print "Cleaning..."
    tclean(vis=active_mss, imagename=prefix+"_img/self"+i, imsize=size, cell=cell, gridder='wproject', wprojplanes=-1, \
        pblimit=0., deconvolver='mtmfs', nterms=3, weighting='briggs', robust=0.5, niter=10000, usemask='auto-multithresh', \
        pbcor=False, savemodel='modelcolumn')

    #clean(vis=active_mss, imagename=prefix+"_img/self1", mode="mfs", gridmode="widefield", wprojplanes=-1, \
    #niter=10000, gain=0.1, threshold="0.0mJy", psfmode="clark", imagermode="csclean", multiscale=[], mask="", \
    #imsize=size, cell=cell, weighting="briggs", robust=0.0, uvtaper=False, outertaper=[''], innertaper=['1.0'], \
    #pbcor=False, minpb=0.2, usescratch=True, nterms=3)

clean(active_mss, '1')
# PHASE 5m
print "5m P cal"
for active_ms in active_mss:
    gaintable = prefix+"_cals/"+active_ms+".Gp1"
    gaincal(vis=active_ms, caltable=gaintable, solint="5min", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="p")
    applycal(vis=active_ms, gaintable=[gaintable], calwt=False, flagbackup=False)

    figure_name=active_ms+'.G1p'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[-1,-1,-180,180], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png', showgui=False)

clean(active_mss, '2')
# PHASE 1m
print "1m P cal"
for active_ms in active_mss:
    gaintable = prefix+"_cals/"+active_ms+".Gp2"
    gaincal(vis=active_ms, caltable=gaintable, solint="1min", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="p")
    applycal(vis=active_ms, gaintable=[gaintable], calwt=False, flagbackup=False)

    figure_name=active_ms+'.G2p'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[-1,-1,-180,180], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png', showgui=False)

clean(active_mss, '3')
# PHASE 30s + AMP 10m
print "30s P + 10m A cal"
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

clean(active_mss, '4')

# then clean with all configurations combined and selfcal with 10min on Gap to aligh amp and phases
gaincal(vis='A1_final.MS', caltable='cals/A1.Gap',solint='10min',calmode='ap')
gaincal(vis='A2_final.MS', caltable='cals/A2.Gap',solint='10min',calmode='ap')
gaincal(vis='B1_final.MS', caltable='cals/B1.Gap',solint='10min',calmode='ap')
gaincal(vis='B2_final.MS', caltable='cals/B2.Gap',solint='10min',calmode='ap')
gaincal(vis='C1_final.MS', caltable='cals/C1.Gap',solint='10min',calmode='ap')
gaincal(vis='C2_final.MS', caltable='cals/C2.Gap',solint='10min',calmode='ap')
