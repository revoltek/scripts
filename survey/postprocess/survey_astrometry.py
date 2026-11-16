#!/usr/bin/env python3
"""
survey_astrometry.py  (from imagemakers/catalogue.ipynb)

Astrometric accuracy of the LoLSS catalogue against FIRST, using isolated and
compact sources in both catalogues.

Outputs (in --outdir): astrometry.pdf, astrometry.json
"""

import os, json, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import norm
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.stats import median_absolute_deviation
import astropy.units as u

import lolss_analysis_common as common
from lib_analysis_lolss import isolated


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cat', default=common.FILE_LOLSS_CAT, help='LoLSS source catalogue (default: %(default)s)')
    parser.add_argument('--ref', default=common.FILE_FIRST_CAT, help='FIRST catalogue (default: %(default)s)')
    parser.add_argument('--compactness-fit', default=common.FILE_COMPACTNESS_FIT,
                        help='Output of survey_compactness.py (default: %(default)s)')
    parser.add_argument('--beam', type=float, default=15., help='LoLSS beam FWHM [arcsec] (default: 15)')
    parser.add_argument('--outdir', default=common.DIR_ANALYSIS, help='Output directory (default: %(default)s)')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ref_t = Table.read(args.ref)
    img_t = Table.read(args.cat)
    print('Initial len: %i (ref:%i)' % (len(img_t), len(ref_t)))

    # reduce to isolated sources
    img_t = isolated(img_t, 3*args.beam)
    ref_t = isolated(ref_t, args.beam)
    print('After isolated sources len: %i (ref:%i)' % (len(img_t), len(ref_t)))

    # restrict to compact
    pnt = common.is_compact(img_t['Total_flux'], img_t['Peak_flux'], img_t['E_Total_flux'],
                            common.load_compact_lolss(args.compactness_fit))
    print('Is pnt %f%%' % (100*np.sum(pnt)/len(img_t)))
    img_t = img_t[pnt]
    ref_t = ref_t[(ref_t['FINT']/ref_t['FPEAK']) < 1.2]
    print('After compact source len: %i (ref:%i)' % (len(img_t), len(ref_t)))

    # cross match, keeping the pairs within 10 times the (iterated) MAD of the separation
    idx_match, sep, _ = match_coordinates_sky(SkyCoord(img_t['RA'], img_t['DEC']),
                                              SkyCoord(ref_t['RA'], ref_t['DEC']))
    sep_mad = median_absolute_deviation(sep[np.where(sep < 100*u.arcsec)])
    sep_mad_old = 0
    while sep_mad != sep_mad_old:
        print('Sep MAD: %f"' % sep_mad.arcsec)
        sep_mad_old = sep_mad
        sep_mad = median_absolute_deviation(sep[np.where(sep < 10 * sep_mad)])
    sep_max = 10*sep_mad
    idx_match_ref = idx_match[sep < sep_max]
    img_t = img_t[sep < sep_max]
    print('After match source len: %i (within %.2f")' % (len(img_t), sep_max.arcsec))

    # offsets in arcsec; the RA difference is wrapped before projecting it on the sky
    ddec = np.array(ref_t['DEC'][idx_match_ref] - img_t['DEC'])*3600
    dra = np.array(ref_t['RA'][idx_match_ref] - img_t['RA'])
    dra[dra > 180] -= 360
    dra[dra < -180] += 360
    dra = dra*np.cos(np.deg2rad(np.array(img_t['DEC'])))*3600

    meandra, meanddec = np.mean(dra), np.mean(ddec)
    stddra, stdddec = np.std(dra), np.std(ddec)
    print('mean = %.2f" %.2f" - std = %.2f" %.2f"' % (meandra, meanddec, stddra, stdddec))

    # gaussian fit using the 90% of points closest to the mean
    thresh_ra = np.percentile(np.abs(dra - meandra), 90)
    thresh_dec = np.percentile(np.abs(ddec - meanddec), 90)
    ra_mu, ra_sigma = norm.fit(dra[np.abs(dra - meandra) <= thresh_ra])
    dec_mu, dec_sigma = norm.fit(ddec[np.abs(ddec - meanddec) <= thresh_dec])
    print('gaussian fit (central 90%%): RA %.3f" +- %.3f" - Dec %.3f" +- %.3f"' % (ra_mu, ra_sigma, dec_mu, dec_sigma))

    with open(os.path.join(args.outdir, 'astrometry.json'), 'w') as f:
        json.dump(dict(n_sources=len(img_t), match_radius_arcsec=float(sep_max.arcsec),
                       mean_dra=float(meandra), mean_ddec=float(meanddec), std_dra=float(stddra), std_ddec=float(stdddec),
                       fit_dra=[float(ra_mu), float(ra_sigma)], fit_ddec=[float(dec_mu), float(dec_sigma)]), f, indent=1)

    # plot
    cm_inferno = plt.colormaps['inferno']
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes([left, bottom, width, height])
    ax.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes([left, bottom + height, width, 0.2])
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes([left + width, bottom, 0.2, height])
    ax_histy.tick_params(direction='in', labelleft=False)
    ax.set_xlabel(r'Differential RA (arcsec)')
    ax.set_ylabel(r'Differential Dec (arcsec)')

    z = common.fast_kde(np.vstack([dra, ddec]))
    ax.scatter(dra, ddec, marker='.', c=z, cmap=cm_inferno, alpha=.5, rasterized=True)
    ax.set_xlim(xmin=-7, xmax=7)
    ax.set_ylim(ymin=-7, ymax=7)
    ax.add_artist(Ellipse(xy=[meandra, meanddec], width=2*stddra, height=2*stdddec, angle=0, alpha=.7,
                          color='white', fill=False, ls='--', lw=2))

    bins = np.arange(-8, 8 + 0.5, 0.15)
    n, _, patches = ax_histx.hist(dra, bins=bins, color='black', alpha=0.9, histtype='bar', density=True)
    common.colour_histogram(n, patches, cm_inferno)
    n, _, patches = ax_histy.hist(ddec, bins=bins, orientation='horizontal', color='black', alpha=0.9, histtype='bar', density=True)
    common.colour_histogram(n, patches, cm_inferno)
    ax_histx.axvline(meandra, color='r', linestyle='dashed', linewidth=1)
    ax_histy.axhline(meanddec, color='r', linestyle='dashed', linewidth=1)
    ax_histx.plot(bins, norm.pdf(bins, ra_mu, ra_sigma), 'k-', linewidth=2)
    ax_histy.plot(norm.pdf(bins, dec_mu, dec_sigma), bins, 'k-', linewidth=2)
    ax_histx.set_yticks([], [])
    ax_histy.set_xticks([], [])
    ax_histx.set_xlim(ax.get_xlim())
    ax_histy.set_ylim(ax.get_ylim())
    ax.axvline(meandra, ls=":", color='r')
    ax.axhline(meanddec, ls=":", color='r')
    fig.savefig(os.path.join(args.outdir, 'astrometry.pdf'), bbox_inches='tight', facecolor='w', dpi=200)


if __name__ == '__main__':
    main()
