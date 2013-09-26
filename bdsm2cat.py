#!/usr/bin/python
"""
Script that takes a PyBDSM fits catalog and cross-
matches to VLSS, WENSS, and NVSS catalogs, derives
spectral indices, position offsets, and flux off-
sets (measured by PyBDSM vs. predicted from catalogs).
   
"""
import pyfits
import os
import sys
import subprocess
import numpy
import math
import pylab as pl
from scipy.optimize import leastsq   

def cone_search(fits_cat, vo_cat, r_search=0.0028, new_fits_cat=None, path_to_stilts=None, use_gsm=False, use_topcat=True):
    """Calls STILTS to do a multi-cone search on a VO catalog

       fits_cat = input PyBDSM fits catalog
       vo_cat = 'vlss', 'wenss', or 'nvss'
       r_search = search radius in degrees (default = 10 arcsec)
       new_fits_cat = name of new fits catalog with cross-matced results
                      (default = None => fits_cat_root + _vo_cat + '.fits'

       Available VO catalogs:
       vlss = http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=VIII/79A&
       wenss = http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=VIII/62&
       nvss = http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=VIII/65&

    """
    if path_to_stilts == None:
        path_to_stilts = 'stilts'
    if vo_cat in ['vlss', 'wenss', 'nvss']:
        if vo_cat == 'vlss':
            vo_url = 'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=VIII/79A&'
        if vo_cat == 'wenss':
            vo_url = 'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=VIII/62&'
        if vo_cat == 'nvss':
            vo_url = 'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=VIII/65&'

        if new_fits_cat == None:
            new_fits_cat = os.path.splitext(fits_cat)[0] + '_' + vo_cat + os.path.splitext(fits_cat)[1]

        # Call STILTS to perform search. All sources are kept. Ignore failed VO queries
        if use_topcat:
            cmd = ['topcat', '-stilts', 'coneskymatch', 'serviceurl='+vo_url, 'in='+fits_cat, 'ra=RA', 'dec=DEC', 'sr='+str(r_search), 'out='+new_fits_cat, 'find=each', 'erract=ignore']        
        else:
            cmd = [path_to_stilts, 'coneskymatch', 'serviceurl='+vo_url, 'in='+fits_cat, 'ra=RA', 'dec=DEC', 'sr='+str(r_search), 'out='+new_fits_cat, 'find=each', 'erract=ignore']
        p = subprocess.call(cmd)
        return new_fits_cat
    else:
        print 'Only "vlss", "wenss", and "nvss" queries are supported.'

        
def fit_spin(freq, flux, eflux, freq0, order=1):
    """Fits a simple power law to SED"""
    if isinstance(freq, list): freq = numpy.array(freq)
    if isinstance(flux, list): flux = numpy.array(flux)
    if isinstance(eflux, list): eflux = numpy.array(eflux)
    xfit = numpy.log10(freq/freq0); yfit = numpy.log10(flux); sig = numpy.abs(eflux/flux)/2.303
    p0 = numpy.array([0]*(order+1))
    p0[1] = (yfit[0]-yfit[-1])/(xfit[0]-xfit[-1])
    p0[0] = yfit[0]-p0[1]*xfit[0]
    if order==2:
        p0[2] = 0.0

    res = lambda p, xfit, yfit, sigfit: (yfit-poly(p, xfit))/sig
    (spin, cov, info, mesg, flag) = leastsq(res, p0, args=(xfit, yfit, sig), full_output=True)
    if cov != None:
        if numpy.sum(sig != 1.) > 0:
            espin = numpy.array([numpy.sqrt(cov[i,i]) for i in range(len(spin))])
        else:
            chisq = sum(info["fvec"]*info["fvec"])
            dof = len(info["fvec"])-len(spin)
            espin = numpy.array([numpy.sqrt(cov[i,i]*chisq/dof) for i in range(len(spin))])
    else:
        spin, espin = [numpy.nan, numpy.nan], [numpy.nan, numpy.nan]

    return spin, espin

def poly(c,x):
    """ y = Sum { c(i)*x^i }, i=0,len(c)"""
    y = numpy.zeros(len(x))
    for i in range(len(c)):
        y += c[i]*(x**i)
    return y

def angsep(ra1deg, dec1deg, ra2deg, dec2deg):
    """Returns angular separation between two coordinates (all in degrees)"""
    import math
    
    ra1rad=ra1deg*math.pi/180.0
    dec1rad=dec1deg*math.pi/180.0
    ra2rad=ra2deg*math.pi/180.0
    dec2rad=dec2deg*math.pi/180.0
    
    # calculate scalar product for determination
    # of angular separation
    x=math.cos(ra1rad)*math.cos(dec1rad)*math.cos(ra2rad)*math.cos(dec2rad)
    y=math.sin(ra1rad)*math.cos(dec1rad)*math.sin(ra2rad)*math.cos(dec2rad)
    z=math.sin(dec1rad)*math.sin(dec2rad)
    
    rad=math.acos(x+y+z)
        
    # use Pythargoras approximation if rad < 1 arcsec
    if rad<0.000004848:
        rad=math.sqrt((math.cos(dec1rad)*(ra1rad-ra2rad))**2+(dec1rad-dec2rad)**2)
        
    # Angular separation
    deg=rad*180/math.pi
    return deg

def plot_sed(freq, flux, eflux, spin, src_id, freq0, order):
    """Plots source SED and overlays best-fit powerlaw."""
    if isinstance(freq, list): freq = numpy.array(freq)
    if isinstance(flux, list): flux = numpy.array(flux)
    if isinstance(eflux, list): eflux = numpy.array(eflux)
    y = flux
    ey = eflux
    x = freq

    fig = pl.figure(figsize=(5.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('SED of Gaussian #'+str(src_id))
    ax1.plot(numpy.log10(x), numpy.log10(y), '*b')
    ax1.errorbar(numpy.log10(x), numpy.log10(y), ey/y)
    xsim = (numpy.array(range(100)) * max(x)/100.0) + min(x) 
    if order == 1:
        ax1.plot(numpy.log10(x), spin[0]+numpy.log10(x/freq0)*spin[1], '-g')
    else:
        ax1.plot(numpy.log10(xsim), spin[0]+numpy.log10(xsim/freq0)*spin[1]+
                 (numpy.log10(xsim/freq0))**2*spin[2], '-g')
    pl.xlabel('log Frequncy (Hz)')
    pl.ylabel('log Flux (Jy)')
    pl.savefig('Gaussian_'+str(src_id)+'_SED_fit.pdf', format='pdf')

def plot_extrap_to_meas_flux_ratio(fp, fp_err, fl, fl_err, ra, dec, ra_pc, 
                                   dec_pc, gaus_id, has_vlss, freq0, snr):
    import math
    ratio = fp/fl
    ratio_err = numpy.zeros(len(fp))
    sep = numpy.zeros(len(fp))
    for i in range(len(fp)):
        if fp_err[i] > fp[i]:
            fp_err[i] = fp[i]
        ratio_err[i] = ratio[i] * math.sqrt((fp_err[i]/fp[i])**2 + (fl_err[i]/fl[i])**2)
        sep[i] = angsep(ra[i], dec[i], ra_pc, dec_pc)
    high_snr = numpy.where(ratio_err/ratio < 1.0/snr)
    print high_snr
    fig = pl.figure(figsize=(7.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    ax1.errorbar(sep[high_snr], ratio[high_snr], yerr=ratio_err[high_snr], fmt='o')
    #ax1.set_yscale('log')
    ax1.set_ylim(0, 3)
    high_snr_and_vlss = numpy.where(has_vlss[high_snr] == True)
    ax1.plot(sep[high_snr][high_snr_and_vlss], ratio[high_snr][high_snr_and_vlss], 's', markersize=10, mfc='None', label='VLSS detection')
    pl.title('A2256 Field: %.1f MHz Flux Comparison\n(All sources with S/N or ratio > %.1f)' % (freq0/1e6, snr))
    pl.ylabel('Extrap. Flux (NVSS+WENSS+VLSS) / Meas. Flux')
    pl.xlabel('Distance from Phase Center (deg)')
    ax1.legend(numpoints=1)
    
    # Calculate mean ratio and std dev.
    mean_ratio = numpy.mean(ratio[high_snr])
    std = numpy.std(ratio[high_snr])
    print 'SNR > %.1f sources: Mean ratio = %.2f, std. dev. = %.2f' % (snr, mean_ratio, std)
#     ax1.plot([0.0, 8.0], [mean_ratio, mean_ratio], '--g')
#     ax1.plot([0.0, 8.0], [mean_ratio+std, mean_ratio+std], '-.g')
#     ax1.plot([0.0, 8.0], [mean_ratio-std, mean_ratio-std], '-.g')
    pl.savefig('extrap_meas_flux_ratio.pdf', format='pdf')
    mean_ratio = numpy.mean(ratio[high_snr][high_snr_and_vlss])
    std = numpy.std(ratio[high_snr][high_snr_and_vlss])
    print 'VLSS detections: Mean ratio = %.2f, std. dev. = %.2f' % (mean_ratio, std)
    
    # Next, plot ratio as function of position
    import matplotlib 
    import matplotlib.patches as mpatches
    from matplotlib.collections import PatchCollection 
    
    fig2=pl.figure() 
    ax=fig2.add_subplot(111) 
    
    resolution = 50 # the number of vertices 
    N = 20 
    x       = ra[high_snr]
    y       = dec[high_snr]
    radii   = ratio_err[high_snr]/ratio[high_snr]*10.0
    colors  = ratio[high_snr]
    verts   = [] 
    for x1,y1,r in zip(x, y, radii): 
        circle = mpatches.Circle((x1,y1), r, ec="none")
        verts.append(circle) 
    
    p = PatchCollection(verts, cmap=matplotlib.cm.jet)
    p.set_array(pl.array(colors)) 
    ax.add_collection(p) 
    pl.colorbar(p) 
    
    ax.axis('equal') 
    pl.xlabel('RA (deg)')
    pl.ylabel('Dec (deg)')
    pl.savefig('extrap_meas_flux_ratio_ra_dec.pdf', format='pdf')

def plot_flux_ratios(cat):
    col_names = cat.dtype.names
    flux_ratio_vlss = []
    flux_ratio_wenss = []
    flux_ratio_nvss = []
    dist_vlss = []
    dist_wenss = []
    dist_nvss = []

    for i, entry in enumerate(cat):
        ra_pybdsm = entry[col_names.index('RA')]
        dec_pybdsm = entry[col_names.index('DEC')]
        flux_pybdsm = entry[col_names.index('RA')]

    flux_pybdsm_vlss = [] # flux in Jy
    eflux_pybdsm_vlss = [] # error on flux in Jy
    flux_pybdsm_nvss = [] # flux in Jy
    eflux_pybdsm_nvss = [] # error on flux in Jy
    flux_pybdsm_wenss = [] # flux in Jy
    eflux_pybdsm_wenss = [] # error on flux in Jy
    flux_vlss = [] # flux in Jy
    eflux_vlss = [] # error on flux in Jy
    flux_wenss = [] # flux in Jy
    eflux_wenss = [] # error on flux in Jy
    flux_nvss = [] # flux in Jy
    eflux_nvss = [] # error on flux in Jy

    for i, entry in enumerate(cat):
        # First, put LOFAR fluxes in array

        # Check for VLSS fluxes
        if 'Si' in col_names:
            fl = entry[col_names.index('Si')]
            efl = entry[col_names.index('e_Si')]
            if not numpy.isnan(fl):
                if fl > 0.0 and efl > 0.0:
                    flux_vlss.append(fl) # Jy
                    eflux_vlss.append(efl) # Jy
                    flux_pybdsm_vlss.append(entry['Total_flux']) # Jy
                    eflux_pybdsm_vlss.append(0.0) # Jy

        # Check for WENSS fluxes
        if 'Sint' in col_names:
            fl = entry[col_names.index('Sint')]
            if not numpy.isnan(fl):
                if fl > 0.0:
                    flux_wenss.append(fl/1000.0) # Jy
                    eflux_wenss.append(fl/1000.0*0.1) # No errors! Assume 10%
                    flux_pybdsm_wenss.append(entry['Total_flux']) # Jy
                    eflux_pybdsm_wenss.append(0.0) # Jy

        # Check for NVSS fluxes
        if 'S1.4' in col_names:
            fl = entry[col_names.index('S1.4')]
            efl = entry[col_names.index('e_S1.4')]
            if not numpy.isnan(fl):
                if fl > 0.0 and efl > 0.0:
                    flux_nvss.append(fl/1000.0) # Jy
                    eflux_nvss.append(efl/1000.0) # Jy
                    flux_pybdsm_nvss.append(entry['Total_flux']) # Jy
                    eflux_pybdsm_nvss.append(0.0) # Jy

    fig = pl.figure(figsize=(5.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('Fluxes: PyBDSM vs. VLSS')
    ax1.plot(flux_vlss, flux_pybdsm_vlss, '*b')
    xmin, xmax, ymin, ymax = pl.axis()
    ax1.plot([xmin, xmax], [xmin, xmax], '-g')
    pl.ylabel('Flux PyBDSM (Jy)')
    pl.xlabel('Flux VLSS (Jy)')
    pl.savefig('flux_vlss.png')

    fig = pl.figure(figsize=(5.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('Fluxes: PyBDSM vs. WENSS')
    ax1.plot(flux_wenss, flux_pybdsm_wenss, '*b')
    xmin, xmax, ymin, ymax = pl.axis()
    ax1.plot([xmin, xmax], [xmin, xmax], '-g')
    pl.ylabel('Flux PyBDSM (Jy)')
    pl.xlabel('Flux WENSS (Jy)')
    pl.savefig('flux_wenss.png')

    fig = pl.figure(figsize=(5.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('Fluxes: PyBDSM vs. NVSS')
    ax1.plot(flux_nvss, flux_pybdsm_nvss, '*b')
    xmin, xmax, ymin, ymax = pl.axis()
    ax1.plot([xmin, xmax], [xmin, xmax], '-g')
    pl.ylabel('Flux PyBDSM (Jy)')
    pl.xlabel('Flux NVSS (Jy)')
    pl.savefig('flux_nvss.png')

def plot_offsets(cat):
    """Plots RA and DEC offsets and vectors in arcsec for 
    PyBDSM sources relative to VLSS, WENSS, and NVSS."""
    col_names = cat.dtype.names
    ra_pybdsml = []
    dec_pybdsml = []
    flux_pybdsml = []
    ra_vlssl = []
    dec_vlssl = []
    ra_wenssl = []
    dec_wenssl = []
    ra_nvssl = []
    dec_nvssl = []
    ra_offsets_vlss = []
    dec_offsets_vlss = []
    ra_offsets_wenss = []
    dec_offsets_wenss = []
    ra_offsets_nvss = []
    dec_offsets_nvss = []
    for i, entry in enumerate(cat):
        ra_pybdsm = entry[col_names.index('RA')]
        ra_pybdsml.append(ra_pybdsm)
        dec_pybdsm = entry[col_names.index('DEC')]
        dec_pybdsml.append(dec_pybdsm)
        flux_pybdsml.append(entry[col_names.index('Total_flux')])
        ra_vlss = entry[col_names.index('_RAJ2000_0')]
        ra_vlssl.append(ra_vlss)
        dec_vlss = entry[col_names.index('_DEJ2000_0')]
        dec_vlssl.append(dec_vlss)
        ra_wenss = entry[col_names.index('_RAJ2000_1')]
        ra_wenssl.append(ra_wenss)
        dec_wenss = entry[col_names.index('_DEJ2000_1')]
        dec_wenssl.append(dec_wenss)
        ra_nvss = entry[col_names.index('_RAJ2000')]
        ra_nvssl.append(ra_nvss)
        dec_nvss = entry[col_names.index('_DEJ2000')]
        dec_nvssl.append(dec_nvss)

        if not numpy.isnan(ra_vlss):
            if ra_pybdsm >= ra_vlss:
                sign = 1.0
            else:
                sign = -1.0
            ra_offsets_vlss.append(sign * angsep(ra_pybdsm, dec_pybdsm, ra_vlss, dec_pybdsm) * 3600.0) # arcsec
            if dec_pybdsm >= dec_vlss:
                sign = 1.0
            else:
                sign = -1.0
            dec_offsets_vlss.append(sign * angsep(ra_pybdsm, dec_pybdsm, ra_pybdsm, dec_vlss) * 3600.0) # arcsec  
        if not numpy.isnan(ra_wenss):
            if ra_pybdsm >= ra_wenss:
                sign = 1.0
            else:
                sign = -1.0
            ra_offsets_wenss.append(sign * angsep(ra_pybdsm, dec_pybdsm, ra_wenss, dec_pybdsm) * 3600.0) # arcsec
            if dec_pybdsm >= dec_wenss:
                sign = 1.0
            else:
                sign = -1.0
            dec_offsets_wenss.append(sign * angsep(ra_pybdsm, dec_pybdsm, ra_pybdsm, dec_wenss) * 3600.0) # arcsec  
        if not numpy.isnan(ra_nvss):
            if ra_pybdsm >= ra_nvss:
                sign = 1.0
            else:
                sign = -1.0
            ra_offsets_nvss.append(sign * angsep(ra_pybdsm, dec_pybdsm, ra_nvss, dec_pybdsm) * 3600.0) # arcsec
            if dec_pybdsm >= dec_nvss:
                sign = 1.0
            else:
                sign = -1.0
            dec_offsets_nvss.append(sign * angsep(ra_pybdsm, dec_pybdsm, ra_pybdsm, dec_nvss) * 3600.0) # arcsec  
            
    fig = pl.figure(figsize=(5.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('Positional offsets: PyBDSM - VLSS')
    ax1.plot(ra_offsets_vlss, dec_offsets_vlss, '*b')
    xmin, xmax, ymin, ymax = pl.axis()
    ax1.plot([xmin, xmax], [0.0, 0.0], '-g')
    ax1.plot([0.0, 0.0], [ymin, ymax], '-g')
    pl.xlabel('RA Offset (arcsec)')
    pl.ylabel('Dec Offset (arcsec)')
    pl.savefig('offsets_vlss.png')
    fig = pl.figure(figsize=(6.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('Positional offset vectors: PyBDSM -> VLSS')
    minra = min(ra_pybdsml)
    mindec = min(dec_pybdsml)
    for i in range(len(ra_offsets_vlss)): 
        if not numpy.isnan(ra_vlssl[i]):
            ax1.arrow((ra_pybdsml[i]-minra)*3600.0, (dec_pybdsml[i]-mindec)*3600.0, -1.0*ra_offsets_vlss[i]*50.0,
                      -1.0*dec_offsets_vlss[i]*50.0, width=500.0*numpy.log10(flux_pybdsml[i]))
            fluxstr = '%.2f Jy' % (flux_pybdsml[i])
            ax1.text((ra_pybdsml[i]-minra)*3600.0, (dec_pybdsml[i]-mindec)*3600.0, fluxstr, size='x-small')
    xmin, xmax, ymin, ymax = pl.axis()
    ax1.plot([xmin, xmax], [0.0, 0.0], '-g')
    ax1.plot([0.0, 0.0], [ymin, ymax], '-g')
    ax1.set_ylim(0, (max(dec_pybdsml)-mindec)*3600.0)
    ax1.set_xlim((max(ra_pybdsml)-minra)*3600.0, 0)
    pl.xlabel('RA (arcsec)')
    pl.ylabel('Dec (arcsec)')
    pl.savefig('offset_vectors_vlss.pdf', format='pdf')
    
    fig = pl.figure(figsize=(5.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('Positional offsets: PyBDSM - WENSS')
    ax1.plot(ra_offsets_wenss, dec_offsets_wenss, '*b')
    xmin, xmax, ymin, ymax = pl.axis()
    ax1.plot([xmin, xmax], [0.0, 0.0], '-g')
    ax1.plot([0.0, 0.0], [ymin, ymax], '-g')
    pl.xlabel('RA Offset (arcsec)')
    pl.ylabel('Dec Offset (arcsec)')
    pl.savefig('offsets_wenss.png')
    fig = pl.figure(figsize=(6.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('Positional offset vectors: PyBDSM -> WENSS')
    minra = min(ra_pybdsml)
    mindec = min(dec_pybdsml)
    for i in range(len(ra_offsets_wenss)): 
        if not numpy.isnan(ra_wenssl[i]):
            ax1.arrow((ra_pybdsml[i]-minra)*3600.0, (dec_pybdsml[i]-mindec)*3600.0, -1.0*ra_offsets_wenss[i]*50.0,
                      -1.0*dec_offsets_wenss[i]*50.0, width=500.0)
            fluxstr = '%.2f Jy' % (flux_pybdsml[i])
            ax1.text((ra_pybdsml[i]-minra)*3600.0, (dec_pybdsml[i]-mindec)*3600.0, fluxstr, size='x-small')
    xmin, xmax, ymin, ymax = pl.axis()
    ax1.plot([xmin, xmax], [0.0, 0.0], '-g')
    ax1.plot([0.0, 0.0], [ymin, ymax], '-g')
    ax1.set_ylim(0, (max(dec_pybdsml)-mindec)*3600.0)
    ax1.set_xlim((max(ra_pybdsml)-minra)*3600.0, 0)
    pl.xlabel('RA (arcsec)')
    pl.ylabel('Dec (arcsec)')
    pl.savefig('offset_vectors_wenss.pdf', format='pdf')

    fig = pl.figure(figsize=(5.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('Positional offsets: PyBDSM - NVSS')
    ax1.plot(ra_offsets_nvss, dec_offsets_nvss, '*b')
    print "Mean RA offsets NVSS: ", numpy.mean(ra_offsets_nvss)
    print "Mean DEC offsets NVSS: ", numpy.mean(dec_offsets_nvss)
    xmin, xmax, ymin, ymax = pl.axis()
    ax1.plot([xmin, xmax], [0.0, 0.0], '-g')
    ax1.plot([0.0, 0.0], [ymin, ymax], '-g')
    pl.xlabel('RA Offset (arcsec)')
    pl.ylabel('Dec Offset (arcsec)')
    pl.savefig('offsets_nvss.png')

    fig = pl.figure(figsize=(6.0,5.0))
    ax1 = pl.subplot(1, 1, 1)
    pl.title('Positional offset vectors: PyBDSM -> NVSS')
    minra = min(ra_pybdsml)
    mindec = min(dec_pybdsml)
    for i in range(len(ra_offsets_nvss)): 
        if not numpy.isnan(ra_nvssl[i]):
            ax1.arrow((ra_pybdsml[i]-minra)*3600.0, (dec_pybdsml[i]-mindec)*3600.0, -1.0*ra_offsets_nvss[i]*50.0,
                      -1.0*dec_offsets_nvss[i]*50.0, width=500.0)
    xmin, xmax, ymin, ymax = pl.axis()
    ax1.plot([xmin, xmax], [0.0, 0.0], '-g')
    ax1.plot([0.0, 0.0], [ymin, ymax], '-g')
    ax1.set_ylim(0, (max(dec_pybdsml)-mindec)*3600.0)
    ax1.set_xlim((max(ra_pybdsml)-minra)*3600.0, 0)
    pl.xlabel('RA (arcsec)')
    pl.ylabel('Dec (arcsec)')
    pl.savefig('offset_vectors_nvss.pdf', format='pdf')


if __name__=='__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <PyBDSM fits catalog>')
    parser.add_option('--output_cat', dest='output_cat', help='Name of output fits catalog contianing spectral indices and predicted fluxes; default = None => infile_root + "_specindex.fits"', metavar='VAL', default=None)
    parser.add_option('--sr', dest='sr', help='Search radius in arcsec; default = 10 arcsec', metavar='VAL', default='10.0')
    parser.add_option('--snr', dest='snr', help='Signal to noise of ratio of extrapolated to measured flux; default = 10', metavar='VAL', default='10.0')
    parser.add_option('--order', dest='order', help='Order of polynomial to fit; default = 1', metavar='VAL', default='1')
    parser.add_option('--freq0', dest='freq0', help='Reference frequency of LOFAR catalog fluxes (Hz); default = None => get from catalog', metavar='VAL', default=None)
    parser.add_option('--stilts_path', dest='spath', help='Full path to STILTS script; default = None => assume in current path', metavar='VAL', default=None)
    parser.add_option('--phase_center', dest='pc', help='Phase center of parent LOFAR image as "RA, Dec" in degrees', metavar='VAL', default=None)
    parser.add_option('-t', action='store_true', dest='topcat', help='Use TOPCAT instead of STILTS', default=False)
    parser.add_option('-c', action='store_true', dest='clobber', help='Clobber any existing files', default=False)
    parser.add_option('-p', action='store_true', dest='plot', help='Skip plotting of SEDs and offsets', default=False)
    parser.add_option('-m', action='store_true', dest='match', help='Skip cross matching of catalogs', default=False)
    parser.add_option('-k', action='store_true', dest='keep_cat', help='Keep fits catalog with cross matches', default=False)
    (options, args) = parser.parse_args()
    if len(args) == 1:
        infile = args[0]
        outfile = options.output_cat
        if outfile == None:
            outfile = os.path.splitext(infile)[0] + '_specindex' + os.path.splitext(infile)[1]
        sr = options.sr
        snr = float(options.snr)
        order = int(options.order)
        r_search = float(sr)/3600.0 # degrees
        freq0 = options.freq0
        path_to_stilts = options.spath
        clobber = options.clobber
        do_plots = not options.plot
        do_match = not options.match
        keep_cat = options.keep_cat
        use_topcat = options.topcat
        phase_center = options.pc
        
        # Cross match sources with VLSS, WENSS, and NVSS
        if do_match:
            print 'Cross matching PyBDSM Gaussians with VLSS, WENSS, and NVSS catalogs...'
            infile1 = infile
            infile2 = cone_search(infile1, 'vlss', r_search=r_search, new_fits_cat=None, path_to_stilts=path_to_stilts, use_topcat=use_topcat)
            infile3 = cone_search(infile2, 'wenss', r_search=r_search, new_fits_cat=None, path_to_stilts=path_to_stilts, use_topcat=use_topcat)
            cm_cat_file = cone_search(infile3, 'nvss', r_search=r_search, new_fits_cat=None, path_to_stilts=path_to_stilts, use_topcat=use_topcat)
            print '...done.'
        else:
            cm_cat_file = os.path.splitext(infile)[0] + '_vlss_wenss_nvss' + os.path.splitext(infile)[1]
            infile2 = os.path.splitext(infile)[0] + '_' + 'vlss' + os.path.splitext(infile)[1]
            infile3 = os.path.splitext(infile2)[0] + '_' + 'nvss' + os.path.splitext(infile2)[1]


        # Construct SEDs
        print 'Fitting for spectral indices...'
        cat = pyfits.getdata(cm_cat_file, 1)
        hdr = pyfits.getheader(infile, 1)
        cdtype = cat.dtype
        col_names = cdtype.names
        if freq0 == None:
            try:
                freq0 = hdr['FREQ0'] # Hz
		print 'Frequency set to [Hz]: '+str(freq0)
            except:
                sys.exit("No frequency found in header. Please specify with --freq0=VAL.")
        else:
            freq0 = float(freq0)
        if do_plots:
            plot_offsets(cat)
            plot_flux_ratios(cat)
      

        gaul_id = numpy.zeros(len(cat), dtype=int)
        spin = numpy.zeros(len(cat))
        espin = numpy.zeros(len(cat))
        spin2 = numpy.zeros(len(cat))
        espin2 = numpy.zeros(len(cat))
        spnorm = numpy.zeros(len(cat))
        espnorm = numpy.zeros(len(cat))
        spinL = numpy.zeros(len(cat))
        espinL = numpy.zeros(len(cat))
        spnormL = numpy.zeros(len(cat))
        espnormL = numpy.zeros(len(cat))
        has_vlss = numpy.zeros(len(cat), dtype=bool)
        for i, entry in enumerate(cat):
            freq = [] # frequency in Hz
            flux = [] # flux in Jy
            eflux = [] # error on flux in Jy

            # First, put LOFAR fluxes in array
            try:
                gaul_id[i] = entry['Gaus_id']
                is_srl = False
            except:
                gaul_id[i] = entry['Source_id']
                is_srl = True
            flux.append(entry['Total_flux']) # Jy *0.85
            eflux.append(entry['E_Total_flux']) # Jy
            freq.append(freq0)

            # Check for VLSS fluxes
            if 'Si' in col_names:
                fl = entry[col_names.index('Si')]
                efl = entry[col_names.index('e_Si')]
                if not numpy.isnan(fl):
                    if fl > 0.0 and efl > 0.0:
                        flux.append(fl) # Jy
                        eflux.append(efl) # Jy
                        freq.append(7.4e7)
                        has_vlss[i] = True

            # Check for WENSS fluxes
            # Flags     : First flag is source type (S,M,E,C).
            #       S-Single component source.
            #       M-Multicomponent source.
            #       C-Component of a multicomponent source.
            #       E-Extende source (more than four components). 
            # The second flag is a warning flag ('*'), iindicating problems in fitting 
            # the source. Source parameters were obtained using 'aperture' integration.
            if 'Sint' in col_names:
                fl = entry[col_names.index('Sint')]
                flag = entry[col_names.index('flg1')]
                if not numpy.isnan(fl):
                    if fl > 0.0 and flag == 'S':
                        flux.append(fl/1000.0) # Jy
                        eflux.append(fl/1000.0*0.1) # No errors! Assume 10%
                        freq.append(3.30e8)

            # Check for NVSS fluxes
            if 'S1.4' in col_names:
                fl = entry[col_names.index('S1.4')]
                efl = entry[col_names.index('e_S1.4')]
                if not numpy.isnan(fl):
                    if fl > 0.0 and efl > 0.0:
                        flux.append(fl/1000.0) # Jy
                        eflux.append(efl/1000.0) # Jy
                        freq.append(1.4e9)
        
            # Fit SEDs to get global spectral indices (incl. LOFAR)
#             if len(freq) > 1:
#                 sp, esp = fit_spin(freq, flux, eflux, freq0)
#                 spnormL[i] = math.pow(10.0, sp[0])
#                 espnormL[i] = abs(sp[0]*math.log(10.0)*esp[0])
#                 spinL[i] = sp[1]
#                 espinL[i] = esp[1]
#                 if do_plots and not numpy.isnan(spin[i]):
#                     plot_sed(freq, flux, eflux, sp, gaul_id[i], freq0)
#             else:
#                 spnormL[i] = numpy.nan
#                 espnormL[i] = numpy.nan
#                 spinL[i] = numpy.nan
#                 espinL[i] = numpy.nan

            # Fit SEDs to get global spectral indices (excl. LOFAR)
            order_fit=1
            if len(freq) > 2:
                if len(freq) <= 3:
                    order_fit=1
                if len(freq) > 3 and order == 2:
                    order_fit=2
                sp, esp = fit_spin(freq[1:], flux[1:], eflux[1:], freq0, order_fit)
                spnorm[i] = math.pow(10.0, sp[0])
                espnorm[i] = abs(sp[0]*math.log(10.0)*esp[0])
                spin[i] = sp[1]
                espin[i] = esp[1]
                if order_fit==2:
                    spin2[i] = sp[2]
                    espin2[i] = esp[2]      
                else:
                    spin2[i] = 0.0
                    espin2[i] = 0.0
                if do_plots and not numpy.isnan(spin[i]):
                    plot_sed(freq, flux, eflux, sp, gaul_id[i], freq0, order_fit)
            else:
                spnorm[i] = numpy.nan
                espnorm[i] = numpy.nan
                spin[i] = numpy.nan
                espin[i] = numpy.nan
        print '...done.'
               
        # Make a fits table of spectral index and predicted flux
        c0 = pyfits.Column(name='Gaus_id', format='J', array=gaul_id)
        c1 = pyfits.Column(name='Global_sp_in_excl', format='D', array=spin)
        c2 = pyfits.Column(name='Err_global_sp_in_excl', format='D', array=espin)
        c3 = pyfits.Column(name='Global_sp_in_incl', format='D', array=spinL)
        c4 = pyfits.Column(name='Err_global_sp_in_incl', format='D', array=espinL)
        c5 = pyfits.Column(name='Pred_total_flux', format='D', unit='Jy', array=spnorm)
        c6 = pyfits.Column(name='Err_pred_total_flux', format='D', unit='Jy', array=espnorm)
        tbhdu = pyfits.new_table([c0, c1, c2, c3, c4, c5, c6])
        tbhdu.writeto(outfile, clobber=clobber)
        print 'Wrote '+outfile
        
        # Make plot of extrapolated/measured flux vs. distance from phase center
        if phase_center != None:
            ra_cen = float(phase_center.split(',')[0])
            dec_cen = float(phase_center.split(',')[1])
            flux = [] # flux in Jy
            eflux = [] # error on flux in Jy
            ra = []
            dec = []
            gaus_id = []
            for i, entry in enumerate(cat):
                flux.append(entry['Total_flux']) # Jy *0.85
                eflux.append(entry['E_Total_flux']) # Jy
                ra.append(entry['RA'])
                dec.append(entry['DEC'])
                if is_srl:
                    gaus_id.append(entry['Source_id'])            
                else:
                    gaus_id.append(entry['Gaus_id'])

            plot_extrap_to_meas_flux_ratio(spnorm, espnorm, numpy.array(flux),
                                           numpy.array(eflux), numpy.array(ra),
                                           numpy.array(dec), ra_cen, dec_cen,
                                           gaus_id, has_vlss, freq0, snr)
            
            
        # Clean up
        cmd = ['rm', '-f', infile2]
        p = subprocess.call(cmd)
        cmd = ['rm', '-f', infile3]
        p = subprocess.call(cmd)
        if not keep_cat:
            cmd = ['rm', '-f', cm_cat_file]
            p = subprocess.call(cmd)
     
    else:
        parser.print_help()
        
