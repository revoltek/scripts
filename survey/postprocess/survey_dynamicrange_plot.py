#!/usr/bin/env python3
"""
survey_dynamicrange_plot.py  (from imagemakers/dynamic range.ipynb)

Plot the rms noise measured by survey_dynamicrange.py in annuli around isolated
sources, as a function of distance and in bins of source flux density.

Output (in --outdir): dynrange.png
"""

import os, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from astropy.table import Table
from astropy.stats import median_absolute_deviation

import lolss_analysis_common as common
from survey_dynamicrange import ANNULI_INNER, ANNULI_WIDTH


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dynrange', default=common.DIR_ANALYSIS + 'dynrange.txt',
                        help='Output of survey_dynamicrange.py (default: %(default)s)')
    parser.add_argument('--groups', type=int, default=7, help='Number of flux density bins (default: 7)')
    parser.add_argument('--flux-min', type=float, default=20., help='Minimum flux density [mJy] (default: 20)')
    parser.add_argument('--flux-max', type=float, default=5000., help='Maximum flux density [mJy] (default: 5000)')
    parser.add_argument('--rms', type=float, default=None,
                        help='Reference rms of the survey to mark with a dashed line [mJy/beam], e.g. the median of the rms map')
    parser.add_argument('--outdir', default=common.DIR_ANALYSIS, help='Output directory (default: %(default)s)')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cat_dr = Table.read(args.dynrange, format='ascii')
    print("tot sources: %i" % len(cat_dr))
    cat_dr.rename_column('col1', 'Source_name')
    cat_dr.rename_column('col2', 'Total_flux')
    n_annuli = len(cat_dr.colnames) - 2
    if n_annuli != len(ANNULI_INNER):
        raise SystemExit('%s has %i annuli, survey_dynamicrange.py defines %i' % (args.dynrange, n_annuli, len(ANNULI_INNER)))
    for i in range(n_annuli):
        cat_dr.rename_column('col%i' % (i+3), 'rms_%i' % (i+1))
    cat_dr.sort('Total_flux')
    cat_dr = cat_dr[(cat_dr['Total_flux'] > args.flux_min) & (cat_dr['Total_flux'] < args.flux_max)]

    # mid radius of each annulus, arcmin
    dists = (ANNULI_INNER + ANNULI_WIDTH/2.) / 60.

    fig, ax = plt.subplots(nrows=1, figsize=(5, 5))
    ax.set_xlabel(r'Distance (arcmin)')
    ax.set_ylabel(r'Rms noise (mJy/b)')

    colors = plt.cm.inferno_r(np.linspace(0, 1, args.groups))
    for s, subset in enumerate(np.array_split(cat_dr, args.groups)):
        min_flux, max_flux = np.min(subset['Total_flux']), np.max(subset['Total_flux'])
        n = len(subset)
        print("%i: Min flux: %f mJy - Max flux: %f" % (n, min_flux, max_flux))
        rmss = np.array([np.median(subset['rms_%i' % i]) for i in range(1, n_annuli+1)])*1e3
        # standard error of the median: sqrt(pi/2) * sigma / sqrt(n), with sigma estimated as 1.4826 * MAD
        rmss_err = np.array([median_absolute_deviation(subset['rms_%i' % i]) for i in range(1, n_annuli+1)])*1e3 \
            * 1.4826 * np.sqrt(np.pi/2) / np.sqrt(n)
        ax.plot(dists, rmss, marker='o', label='%.0f - %.0f mJy' % (min_flux, max_flux), c=colors[s], alpha=.8)
        ax.fill_between(dists, rmss-rmss_err, rmss+rmss_err, color=colors[s], alpha=.1)
    ax.legend(loc=1, prop={'size': 8}, ncol=2)

    if args.rms is not None:
        ax.axhline(args.rms, ls="--", color='k')
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.set_xlim(xmin=0, xmax=dists[-1] + 1)
    fig.savefig(os.path.join(args.outdir, 'dynrange.png'), bbox_inches='tight')


if __name__ == '__main__':
    main()
