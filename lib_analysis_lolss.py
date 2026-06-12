import os, sys, glob
import scipy
import numpy as np
from scipy.stats import norm
from astropy.table import Table, Column, vstack
import astropy.units as u
from itertools import cycle
from scipy.stats import gaussian_kde
from astropy.stats import median_absolute_deviation
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Ellipse

def isolated(tab, dist):
    # reduce to isolated sources - nothing closer than $dist arcsec
    _, sep, _ = match_coordinates_sky(SkyCoord(tab['RA'], tab['DEC']),
                                      SkyCoord(tab['RA'], tab['DEC']), nthneighbor=2)
    idx_match = np.arange(0, len(tab))[sep > dist*u.arcsec]
    print('Removing %i out of %i sources' % (len(tab) - len(idx_match), len(tab)))
    return tab[idx_match]

def cross_match(tab1, tab2, dist, columns_to_keep=None, plotname=None):
    """
    add to tab1 the matching sources of tab 2
    columns_to_keep shound be the column name of e.g. the flux density in tab2 that will be added in tab1
    tab1 with new columns is eventually returned
    """
    tab1 = tab1.copy()  # avoid modifying original table

    idx_match, sep, _ = match_coordinates_sky(SkyCoord(tab1['RA'], tab1['DEC']),\
                                          SkyCoord(tab2['RA'], tab2['DEC']))
    if columns_to_keep is None:
        columns_to_keep = []
    idx_match_2 = idx_match[sep<dist*u.arcsec]
    idx_match_1 = np.arange(0,len(tab1))[sep<dist*u.arcsec]

    # make histogram of separations
    if plotname is not None:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'separation (arcsec)')
        ax.hist(sep[sep<dist*u.arcsec].arcsec, bins=30, orientation='vertical', color='black', alpha=0.5, histtype='bar', ec='black')
        fig.savefig(plotname, bbox_inches='tight')

    # add column
    for column_to_keep in columns_to_keep:
        c = Column(name=column_to_keep, dtype=tab2[column_to_keep].dtype, length=len(tab1))
        tab1.add_column(c)
        tab1[column_to_keep][idx_match_1] = tab2[column_to_keep][idx_match_2]

    return tab1

def cross_match_multicomp(tab1, tab2, dist, flux_col, columns_to_keep=None):
    """
    Cross-match two radio catalogues with different resolutions.
    
    tab1: lower-resolution catalogue (one source may match multiple in tab2)
    tab2: higher-resolution catalogue (multiple components summed into one match)
    dist: matching radius in arcsec
    flux_col: name of the flux density column in tab2 to be summed
    columns_to_keep: other columns from tab2 to keep (nearest-match value used, not summed)
    
    Returns tab1 with new columns added.
    """
    tab1 = tab1.copy()  # avoid modifying original table

    if columns_to_keep is None:
        columns_to_keep = []

    coord1 = SkyCoord(tab1['RA'], tab1['DEC'], unit='deg')
    coord2 = SkyCoord(tab2['RA'], tab2['DEC'], unit='deg')

    # search_around_sky returns ALL pairs within dist, not just nearest
    # idx1: indices into tab1, idx2: indices into tab2
    idx1, idx2, sep, _ = search_around_sky(coord1, coord2, dist * u.arcsec)

    # --- summed flux column ---
    flux_summed = Column(name=flux_col+'_sum', data=np.full(len(tab1), np.nan))
    tab1.add_column(flux_summed)

    # --- number of components matched ---
    try:
        n_comp = Column(name='n_comp1', data=np.zeros(len(tab1), dtype=int))
        tab1.add_column(n_comp)
    except:
        n_comp = Column(name='n_comp2', data=np.zeros(len(tab1), dtype=int))
        tab1.add_column(n_comp)

    # --- other columns (nearest match) ---
    for col in columns_to_keep:
        c = Column(name=col, dtype=tab2[col].dtype, length=len(tab1))
        tab1.add_column(c)

    # group tab2 matches by their tab1 counterpart
    for i1 in np.unique(idx1):
        mask = idx1 == i1
        matched_idx2 = idx2[mask]
        matched_sep  = sep[mask]

        # sum flux of all tab2 components matched to this tab1 source
        tab1[flux_col+'_sum'][i1] = np.sum(tab2[flux_col][matched_idx2])
        tab1[n_comp.name][i1] = len(matched_idx2)

        # for other columns, take the value from the nearest match
        nearest = matched_idx2[np.argmin(matched_sep)]
        for col in columns_to_keep:
            tab1[col][i1] = tab2[col][nearest]

    return tab1

# extimate errors and accept errors on ydata
def linearfit(x, y, yerr=None, tolog=False):
    # Using OLS (X|Y)" # for more algo read: Isobe et al 1990
    # tolog : convert in log space x, y, and yerr before doing linear regression
    from scipy.optimize import curve_fit

    def f(x, B0, B1):
        return B0*x + B1

    if tolog:
        if yerr is not None: yerr = 0.434*np.array(yerr)/y
        x=np.log10(x)
        y=np.log10(y)
    out = curve_fit(f, x, y, [-1.5 ,1], yerr)
    # return B0, B1, errB0, errB1 (err are in std dev)
    if isinstance(out[1], np.ndarray):
        return (out[0][0], out[0][1], np.sqrt(out[1][0][0]), np.sqrt(out[1][1][1]))
    else:
        return (out[0][0], out[0][1])

def get_err(flux, rms, scale_err=.1):
    return np.sqrt( (flux*scale_err)**2 + rms**2 )

def median_err(data):
    sigma = np.std(data)
    n = len(data)
    sigma_median = 1.253*sigma/np.sqrt(n)
    return sigma_median
