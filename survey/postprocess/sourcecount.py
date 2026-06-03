#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
from matplotlib.ticker import ScalarFormatter, MaxNLocator, LogLocator
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyregion

# LoTSS Deep Fields source count coefficients from Mandal et al. (2021)
COEFFS = [1.655, -0.1150, 0.2272, 0.51788, -0.449661, 0.160265, -0.028541, 0.002041]
MAX_SEPARATION = 4. # in degrees, to match the primary beam radius
PHASECENTER = SkyCoord(ra=242.7519347, dec=54.94888887, unit=u.deg) # ELAIS-N1 phase center
FIGSIZE = 7
FONTSIZE = 12
MOCKS_EXIST = False

def get_frac_spline(bin_centers):
    from scipy.interpolate import Akima1DInterpolator
    det_flxs, meas_flxs, cat_flxs = get_matched_flux_histogram(max_sep=MAX_SEPARATION)
    fracs = np.nanmedian(meas_flxs/cat_flxs, axis=0)
    fracs[np.isnan(fracs)] = 1.
    frac_spline = Akima1DInterpolator(bin_centers[np.isfinite(fracs)], fracs[np.isfinite(fracs)], method="makima")
    return frac_spline

def draw_random_histvals_data(fluxes, flux_scale_error=0.1, uncertainty=0, bins=np.linspace(-1, 5, 30), draws=10):
    uncertainty = np.sqrt((flux_scale_error * fluxes)**2 + uncertainty**2)
    sampled_fluxes = np.zeros((len(fluxes), draws))
    result = np.zeros((draws, len(bins_data)-1))
    for i, flux in enumerate(fluxes):
        sampled_fluxes[i] = np.random.normal(loc=flux, scale=uncertainty[i], size=draws)
    for d in range(draws):
        result[d] = np.histogram(sampled_fluxes[:,d], bins=10**bins_data)[0]
    return result

def source_counts(log_flux: float, coeff:list) -> float:
    """Return the differential source counts dN/dS at a given flux (in mJy)"""
    #log_flux = np.log10(flux)
    log_dnds = 0
    for i in range(len(coeff)):
        log_dnds += coeff[i] * log_flux**i
    return log_dnds 


srl = fits.open('pybdsf_output_wideDDP-merged-MFS-image-pb/wideDDP-merged-MFS-image-pb_catalog.srl')

table = Table(srl[1].data)
separation = PHASECENTER.separation(SkyCoord(ra=table['RA'], dec=table['DEC'], unit=u.deg))
table = table[separation < MAX_SEPARATION*u.deg]
fluxes = np.asarray(table['Total_flux'])
fluxes_err = np.asarray(table['E_Total_flux'])

bins_data = np.linspace(-4, 2, 50)
histvals, bin_edges = np.histogram(fluxes, bins=10**bins_data)
bin_centers_data = np.asarray([bin_edges[i] + (bin_edges[i+1]-bin_edges[i])/2 for i in range(len(bin_edges)-1)])

area_sr = 2 * np.pi * (1 - np.cos(np.radians(MAX_SEPARATION*u.deg))) * u.sr
histvals_enc = histvals / (area_sr.value * (bins_data[1]-bins_data[0]) * np.log(10) * bin_centers_data**-1.5)


# MC sample generation to estimate error or binning.
n_draws = 4000
result = draw_random_histvals_data(fluxes, uncertainty=fluxes_err, flux_scale_error=0.1, bins=bins_data, draws=n_draws)

result_sr = result / (area_sr.value * (bins_data[1]-bins_data[0]) * np.log(10))
result_eucl_norm = np.zeros_like(result)
for d in range(n_draws):
    result_eucl_norm[d] = result_sr[d] / bin_centers_data**-1.5

percentiles = np.percentile(result_eucl_norm, [50-34.1, 50, 50+34.1], axis=0)
vals = percentiles[1]
err = (percentiles[2]-percentiles[0]) / 2


# correct for incompleteness using the mocks
if MOCKS_EXIST:
    frac_spline = get_frac_spline(bin_centers_data)
    corrected_vals = histvals_enc / frac_spline(np.log10(bin_centers_data *(150./50)**-0.7))
    

# LoLSS dN/dS from de Gasperin et al. (2023)
import pandas as pd
lolss_dNdS = pd.read_csv("nesc_lolss.csv", delimiter=', ')
lolss_dNdS
lolss_dNdS_binscenters = 10**(np.log10(lolss_dNdS['max_flux'].values) - np.diff(np.log10(lolss_dNdS['max_flux'].values))[0]/2)
lolss_dNdS_errors = lolss_dNdS['dNdS_err'].values
lolss_dNdS = lolss_dNdS['dNdS'].values


# figure
fig, ax = plt.subplots(figsize=(FIGSIZE, FIGSIZE*3./4))

x_arr = np.linspace(-2,4.7,100)
source_counts_vals = 10**source_counts(x_arr, COEFFS)
x_arr = x_arr-3
source_counts_vals *= (10**x_arr)**-1.5 * ((x_arr[1]-x_arr[0]) * np.log(10)) # convert to actual counts
x_arr += np.log10((50./150)**-0.7) # to 50 MHz
source_counts_vals /= ((x_arr[1]-x_arr[0]) * np.log(10) * (10**x_arr)**-1.5)
ax.plot(10**x_arr*1e3, source_counts_vals, color='black', label='Mandal et al. (2021)')

ax.errorbar(lolss_dNdS_binscenters*1e3, lolss_dNdS, yerr=lolss_dNdS_errors, fmt='s', markersize=5, color='darkorange', label='LoLSS (de Gasperin et al. 2023)')
ax.errorbar(bin_centers_data*1e3, histvals_enc, yerr=err, fmt='.', markersize=12, markeredgecolor='mediumseagreen', color='none', ecolor='mediumseagreen', elinewidth=1.2, markeredgewidth=1.2, label='Raw')
if MOCKS_EXIST:
    frac_spline = get_frac_spline(bin_centers_data)
    ax.errorbar(bin_centers_data*1e3, histvals_enc / frac_spline(np.log10(bin_centers_data)), yerr=err/ frac_spline(np.log10(bin_centers_data)), fmt='.', markersize=12, color='mediumseagreen', elinewidth=1.2, markeredgewidth=1.2, label='Corrected')

ax.set_xlabel(r'$S_{50\;\mathrm{MHz}}$ (mJy beam$^{-1}$)', fontsize=FONTSIZE+2)
ax.set_ylabel('$S^{2.5}dN/dS$ (Jy$^{1.5}$sr$^{-1}$)', fontsize=FONTSIZE+2)
ax.tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=FONTSIZE+2) 
ax.set_xlim(10**bins_data[0]*1e3, 10**bins_data[-2]*1e3)
ax.legend(loc='lower right', frameon=False, fontsize=FONTSIZE)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(7e0, 4e4)
ax.set_xlim(6e-1,1.5e4)
#ax.yaxis.set_major_formatter(ScalarFormatter())
#ax.yaxis.set_major_locator(LogLocator(base=10, numticks=5))
plt.savefig('plots/source_counts.pdf', bbox_inches='tight')