import os, sys, glob
import scipy
import numpy as np
from scipy.stats import norm
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, vstack
import astropy.units as u
from itertools import cycle
from scipy.stats import gaussian_kde
from astropy.stats import median_absolute_deviation

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Ellipse

def isolated(tab, dist):
    # reduce to isolated sources - nothing closer than $dist arcsec
    idx_match, sep, _ = match_coordinates_sky(SkyCoord(tab['RA'], tab['DEC']),\
                                              SkyCoord(tab['RA'], tab['DEC']), nthneighbor=2)
    idx_match = np.arange(0,len(tab))[sep>dist*u.arcsec]
    print('Removing %i out of %i source' % (sum(sep<dist*u.arcsec), len(sep)))
    return tab[idx_match]

def cross_match(tab1, tab2, dist, columns_to_keep=[], plotname=None):
    """
    add to tab1 the matching sources of tab 2
    columns_to_keep shound be the column name of e.g. the flux density in tab2 that will be added in tab1
    tab1 with new columns is eventually returned
    """
    idx_match, sep, _ = match_coordinates_sky(SkyCoord(tab1['RA'], tab1['DEC']),\
                                          SkyCoord(tab2['RA'], tab2['DEC']))
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
    #if yerr is None: yerr = np.ones(len(y))
    #for i,e in enumerate(yerr):
    #    if e == 0: yerr[i] = 1
    out = curve_fit(f, x, y, [-1.5 ,1], yerr)
    # return B0, B1, errB0, errB1 (err are in std dev)
    if type(out[1]) is np.ndarray:
        return (out[0][0], out[0][1], np.sqrt(out[1][0][0]), np.sqrt(out[1][1][1]))
    else:
        return (out[0][0], out[0][1])

def get_err(flux, rms, scale_err=.1):
    return np.sqrt( (flux*scale_err)**2 + rms**2 )

def median_err(data):
     sigma=np.std(data)
     n=len(data)
     sigma_median=1.253*sigma/np.sqrt(n)
     return sigma_median

def keep_sources_pointing(cat, pointing_name):
    """
    remove from the table all sources that have a pointing closest than pointing_name
    t needs to have "ra" and "dec" columns
    """
    grid = Table.read('../allsky-grid.fits')
    cat['ra'].unit = u.deg
    cat['dec'].unit = u.deg
    idx_match, sep, _ = match_coordinates_sky(SkyCoord(cat['ra'], cat['dec']), SkyCoord(grid['ra'], grid['dec']))
    return cat[idx_match == list(grid['name']).index(pointing_name.upper())]
