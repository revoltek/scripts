"""
Helpers shared by the analysis scripts that were ported from the notebooks in
imagemakers/ (survey_compactness.py, survey_astrometry.py,
survey_spidx_crossmatch.py, survey_dynamicrange_plot.py,
survey_falsepositive_plot.py).
"""

import os
import sys
import json

import numpy as np

# lib_analysis_lolss.py lives in the root of the scripts repository
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

DIR_SURVEY     = os.path.expanduser('~/storage/surveytgts/')
DIR_ANALYSIS   = DIR_SURVEY + 'analysis/'
DIR_CATALOGUES = DIR_SURVEY + 'catalogues/'
DIR_EXTCAT     = os.path.expanduser('~/storage/catalogues/')

FILE_LOLSS_CAT       = DIR_CATALOGUES + 'LoLSS_DR2-cat.fits'
FILE_LOTSS_CAT       = DIR_EXTCAT + 'LoTSS_DR3_v1.0.srl.fits'
FILE_NVSS_CAT        = DIR_EXTCAT + 'NVSS.fits'
FILE_FIRST_CAT       = os.path.expanduser('~/storage/FIRST_14dec17.fits.gz')
FILE_COMPACTNESS_FIT = DIR_ANALYSIS + 'compactness_fit.json'

# Envelope of unresolved sources: ln(S_int/S_peak) < A + B / (1 + (SNR/C)**D), SNR = S_int / sigma_S_int
# LoTSS: Shimwell et al. (2022)
COMPACT_LOTSS = dict(A=0.42, B=1.08, C=96.57, D=2.49)
# LoLSS: values used so far in the notebooks, superseded by the output of survey_compactness.py
COMPACT_LOLSS_DEFAULT = dict(A=0.14, B=1.02, C=161., D=1.25)


def compact_envelope(snr, coeff):
    """Maximum ln(S_int/S_peak) of an unresolved source with the given S/N."""
    return coeff['A'] + coeff['B'] / (1 + (np.asarray(snr)/coeff['C'])**coeff['D'])


def is_compact(total_flux, peak_flux, e_total_flux, coeff):
    """True for the sources below the envelope of unresolved sources."""
    snr = np.asarray(total_flux) / np.asarray(e_total_flux)
    return compact_envelope(snr, coeff) > np.log(np.asarray(total_flux) / np.asarray(peak_flux))


def load_compact_lolss(filename=FILE_COMPACTNESS_FIT):
    """Coefficients of the LoLSS envelope from survey_compactness.py, or the default ones."""
    if filename is not None and os.path.exists(filename):
        with open(filename) as f:
            coeff = json.load(f)
        print('LoLSS compactness envelope from %s: A=%.3f B=%.3f C=%.1f D=%.3f (percentile %g)'
              % (filename, coeff['A'], coeff['B'], coeff['C'], coeff['D'], coeff.get('percentile', np.nan)))
        return coeff
    print('WARNING: %s not found, using the default LoLSS compactness envelope %s' % (filename, COMPACT_LOLSS_DEFAULT))
    return dict(COMPACT_LOLSS_DEFAULT)


def fast_kde(xy, max_fit=5000, seed=0):
    """Point density for colouring a scatter plot: gaussian_kde fitted on a subsample for speed."""
    from scipy.stats import gaussian_kde
    n = xy.shape[1]
    if n > max_fit:
        idx = np.random.default_rng(seed).choice(n, max_fit, replace=False)
        return gaussian_kde(xy[:, idx])(xy)
    return gaussian_kde(xy)(xy)


def colour_histogram(n, patches, cmap):
    """Colour the bars of a histogram according to their height."""
    col = (n - n.min()) / (n.max() - n.min())
    for c, p in zip(col, patches):
        p.set_facecolor(cmap(c))
