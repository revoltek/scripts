# to be exectured in CASA
import os, sys

prefix = "C"
active_mss = ['C1_target-preflag.MS', 'C2_target-preflag.MS']

# A conf
#cell = '0.5'
#size = 5000
# B conf
#cell = '1.5arcsec'
#size = 2500
# C conf
cell = '2.5arcsec'
size = 1500

os.makedirs(prefix+"_img")
os.makedirs(prefix+"_cals")
os.makedirs(prefix+"_plots")

def clean(active_mss):
    tclean(vis=active_mss, magename=prefix+"_img/self1", imsize=size, cell=cell, gridder='wproject', wprojplanes=-1, \
        pblimit=-1, deconvolver='mtmfs', nterms=3, weighting='briggs', robust=0., niter=10000, usemask='auto-multithresh', \
        pbcor=False, savemodel='modelcolumn')

    #clean(vis=active_mss, imagename=prefix+"_img/self1", mode="mfs", gridmode="widefield", wprojplanes=-1, \
    #niter=10000, gain=0.1, threshold="0.0mJy", psfmode="clark", imagermode="csclean", multiscale=[], mask="", \
    #imsize=size, cell=cell, weighting="briggs", robust=0.0, uvtaper=False, outertaper=[''], innertaper=['1.0'], \
    #pbcor=False, minpb=0.2, usescratch=True, nterms=3)

clean(active_mss)
# PHASE 5m
for active_ms in active_mss
    gaintable = prefix+"_cals/"+active_ms+".Gp1"
    gaincal(vis=active_ms, caltable=gaintable, solint="5min", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="p")
    applycal(vis=active_ms, gaintable=[gaintable], calwt=False, flagbackup=False))

    figure_name='G1p'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[-1,-1,-180,180], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png')

clean(active_mss)
# PHASE 1m
for active_ms in active_mss
    gaintable = prefix+"_cals/"+active_ms+".Gp2"
    gaincal(vis=active_ms, caltable=gaintable, solint="1min", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="p")
    applycal(vis=active_ms, gaintable=[gaintable], calwt=False, flagbackup=False))

    figure_name='G2p'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[-1,-1,-180,180], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png')

clean(active_mss)
# PHASE 30s + AMP 10m
for active_ms in active_mss
    gaintable1 = prefix+"_cals/"+active_ms+".Gp3"
    gaincal(vis=active_ms, caltable=gaintable1, solint="30s", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="p")
    gaintable2 = prefix+"_cals/"+active_ms+".Gap3"
    gaincal(vis=active_ms, caltable=gaintable2, solint="10min", refant="", minblperant=4, minsnr=3.0, gaintype="G", calmode="ap", gaintable=[gaintable1])
    applycal(vis=active_ms, gaintable=[gaintable1, gaintable2], calwt=False, flagbackup=False))

    figure_name='G3p'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable1, xaxis='time', yaxis='phase', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[-1,-1,-180,180], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png')

    figure_name='G3a'
    for i in range(num_plots+1):
        ant_plot=str(i*4)+'~'+str(i*4+3)
        plotcal(caltable=gaintable2, xaxis='time', yaxis='amp', subplot=221, antenna=ant_plot, iteration='antenna', \
                plotrange=[], figfile=prefix+"_plots/"+figure_name+str(i+1)+'.png')

clean(active_mss)

# uvsub of outskirts and final clean
#uvsub()
#tclean(vis=active_mss, magename=prefix+"_img/final", imsize=size/4, cell=cell, gridder='wproject', wprojplanes=-1, \
#    pblimit=0.2, deconvolver='mtmfs', nterms=2, weighting='briggs', robust=0., niter=10000, usemask='auto-multithresh', \
#    pbcor=True)
