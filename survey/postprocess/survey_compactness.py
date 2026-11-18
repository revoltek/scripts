#!/usr/bin/env python3
"""
survey_compactness.py  (from imagemakers/catalogue.ipynb)

Derive the envelope that separates unresolved from resolved sources in the
LoLSS catalogue, in the plane ln(S_int/S_peak) vs S/N.

The envelope is the --percentile of ln(S_int/S_peak), in bins of S/N, of the
LoLSS sources that are isolated, single-Gaussian, smaller than 30" and matched
to an unresolved LoTSS source. It is fitted with
    R = A + B / (1 + (SNR/C)**D)
and the coefficients are written to compactness_fit.json, which is read by
survey_spidx_crossmatch.py.

Outputs (in --outdir): compactness_fit.json, inttopeak.png
"""

import os, json, argparse
import numpy as np
import scipy.optimize
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import lolss_analysis_common as common
from lib_analysis_lolss import isolated, cross_match
from astropy.table import Table
from astropy.stats import median_absolute_deviation


def binning(x, y, nbins, percentile):
    """Median, MAD and percentile of y in nbins-1 bins of x equally spaced in log."""
    bins = np.logspace(np.log10(min(x)), np.log10(max(x)), num=nbins)
    bin_med, bin_edges, _ = stats.binned_statistic(x, y, statistic='median', bins=bins)
    bin_std, _, _ = stats.binned_statistic(x, y, statistic=median_absolute_deviation, bins=bins)
    bin_pct, _, _ = stats.binned_statistic(x, y, statistic=lambda v: np.percentile(v, percentile), bins=bins)
    bin_width = bin_edges[1:] - bin_edges[:-1]
    bin_centers = bin_edges[1:] - bin_width/2
    return bin_centers, bin_width, bin_med, bin_std, bin_pct


def model(x, coeffs):
    return abs(coeffs[3] + (coeffs[0] + coeffs[3])/(1+(x/coeffs[2])**(coeffs[1])))


def residuals(coeffs, y, x):
    return y - model(x, coeffs)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cat', default=common.FILE_LOLSS_CAT, help='LoLSS source catalogue (default: %(default)s)')
    parser.add_argument('--lotss', default=common.FILE_LOTSS_CAT, help='LoTSS source catalogue (default: %(default)s)')
    parser.add_argument('--outdir', default=common.DIR_ANALYSIS, help='Output directory (default: %(default)s)')
    parser.add_argument('--percentile', type=float, default=99.,
                        help='Percentile of ln(S_int/S_peak) of the unresolved sources that defines the envelope (default: 99)')
    parser.add_argument('--snr-min', type=float, default=4., help='Minimum S/N of the sources used for the envelope (default: 4)')
    parser.add_argument('--no-plot', action='store_true', help='Only fit the envelope')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cat_lolss = isolated(Table.read(args.cat), dist=15*2)
    cat_lotss = isolated(Table.read(args.lotss), dist=15*2)
    cat_lotss.rename_column('Total_flux', 'flux_lotss')
    cat_lotss['is_pnt'] = common.is_compact(cat_lotss['flux_lotss'], cat_lotss['Peak_flux'],
                                            cat_lotss['E_Total_flux'], common.COMPACT_LOTSS)
    print('Is pnt (LoTSS) %f%%' % (100*np.sum(cat_lotss['is_pnt'])/len(cat_lotss)))

    cat_full = cross_match(cat_lolss, cat_lotss, dist=15, columns_to_keep=['flux_lotss', 'is_pnt'])
    del cat_lotss
    print('Match LoTSS: %i' % np.sum(cat_full['flux_lotss'] != 0))
    print('un-Match LoTSS: %i' % np.sum(cat_full['flux_lotss'] == 0))

    is_pnt = np.array((cat_full['is_pnt']) & (cat_full['S_Code'] == 'S') & (cat_full['Maj'] < 30))
    print('Point source selection: %i / %i (%.1f%%)' % (is_pnt.sum(), len(cat_full), 100*is_pnt.sum()/len(cat_full)))

    int_to_peak = np.array(np.log(cat_full['Total_flux']/cat_full['Peak_flux']))
    snratio = np.array(cat_full['Total_flux']/cat_full['E_Total_flux'])

    # envelope of the unresolved sources and fit
    selection = is_pnt & (snratio > args.snr_min)
    bin_centers, bin_width, bin_med, bin_std, bin_pct = binning(snratio[selection], int_to_peak[selection],
                                                                nbins=15, percentile=args.percentile)
    x0 = [0.7, 2.0, 100.0, 0.4]
    xcoeff, flag = scipy.optimize.leastsq(residuals, x0, args=(bin_pct, bin_centers))
    coeff = dict(A=float(xcoeff[3]), B=float(xcoeff[0] + xcoeff[3]), C=float(xcoeff[2]), D=float(xcoeff[1]),
                 percentile=args.percentile, snr_min=args.snr_min, n_sources=int(selection.sum()), catalogue=args.cat)
    print('Fit coeff:', xcoeff, flag)
    print('R_{%g} = %f + \\left( \\frac{%f}{1+\\left( \\frac{\\rm SNR}{%f} \\right)^{%f}} \\right)'
          % (args.percentile, coeff['A'], coeff['B'], coeff['C'], coeff['D']))
    with open(os.path.join(args.outdir, 'compactness_fit.json'), 'w') as f:
        json.dump(coeff, f, indent=1)

    resolved = int_to_peak > model(snratio, xcoeff)
    print('Sources above the pntsrc line: %i (%.2f%%)' % (np.sum(resolved), np.sum(resolved)/len(int_to_peak)*100))

    if args.no_plot:
        return

    cm_inferno = plt.colormaps['inferno']
    left, width = 0.1, 1.0
    bottom, height = 0.1, 0.5
    fig = plt.figure(figsize=(5, 7))
    ax1 = plt.axes([left, bottom + height, width, height])
    ax1_histy = plt.axes([left + width, bottom + height, 0.2, height])
    ax2 = plt.axes([left, bottom, width, height])
    ax2_histy = plt.axes([left + width, bottom, 0.2, height])
    for ax in (ax1, ax2):
        ax.tick_params(direction='in', top=True, right=True)
    for ax in (ax1_histy, ax2_histy):
        ax.tick_params(direction='in', labelleft=False)
    ax2.set_xlabel(r'Signal-to-noise ratio')
    ax1.set_ylabel(r'$\ln(S_I / S_P)$')
    ax2.set_ylabel(r'$\ln(S_I / S_P)$')

    print('Computing KDE (all sources, n=%i)...' % len(snratio))
    z = common.fast_kde(np.vstack([snratio, int_to_peak]))
    ax1.scatter(snratio, int_to_peak, marker='.', c=z, cmap=cm_inferno, alpha=.3, rasterized=True)
    print('Computing KDE (point sources, n=%i)...' % is_pnt.sum())
    z = common.fast_kde(np.vstack([snratio[is_pnt], int_to_peak[is_pnt]]))
    ax2.scatter(snratio[is_pnt], int_to_peak[is_pnt], marker='.', c=z, cmap=cm_inferno, alpha=.3, rasterized=True)

    # median and scatter in bins
    for ax, sel in ((ax1, np.ones(len(snratio), dtype=bool)), (ax2, is_pnt)):
        c, w, med, std, _ = binning(snratio[sel], int_to_peak[sel], nbins=10, percentile=args.percentile)
        ax.errorbar(c, med, xerr=w/2., yerr=std, color='r', zorder=100, ls='', lw=2, alpha=0.8)
    ax2.scatter(bin_centers, bin_pct, color='orange', zorder=100, marker='x', alpha=0.8)
    xpoints = np.logspace(np.log10(2.5), np.log10(1000))
    ax2.plot(xpoints, model(xpoints, xcoeff), color='orange', zorder=100, ls=':', marker='', alpha=0.8)

    bins = np.arange(-1, 4 + 0.1, 0.05)
    for ax, values in ((ax1_histy, int_to_peak), (ax2_histy, int_to_peak[is_pnt])):
        n, _, patches = ax.hist(values, bins=bins, orientation='horizontal', color='black', alpha=0.9, histtype='bar')
        ax.set_xticks([], minor=[])
        common.colour_histogram(n, patches, cm_inferno)

    for ax in (ax1, ax2):
        ax.set_xlim(xmin=2.5, xmax=1000)
        ax.set_ylim(ymin=-0.5, ymax=3)
        ax.set_xscale('log')
        ax.axhline(0, ls="--", color='c')
        ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])
    ax1_histy.set_ylim(ax1.get_ylim())
    ax2_histy.set_ylim(ax2.get_ylim())
    fig.savefig(os.path.join(args.outdir, 'inttopeak.png'), bbox_inches='tight', dpi=200)


if __name__ == '__main__':
    main()
