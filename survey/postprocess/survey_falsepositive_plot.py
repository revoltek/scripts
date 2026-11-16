#!/usr/bin/env python3
"""
survey_falsepositive_plot.py  (from imagemakers/completeness.ipynb)

False-positive rate of the catalogue: sources detected in the inverted mosaics
(survey_falsepositive.py + survey_catalogues_concat.py) compared with the real
detections, in total and as a function of flux density.

Outputs (in --outdir): falsepositive_hist.pdf, falsepositive.txt
"""

import os, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from astropy.table import Table

import lolss_analysis_common as common


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cat', default=common.FILE_LOLSS_CAT, help='LoLSS source catalogue (default: %(default)s)')
    parser.add_argument('--cat-inv', default=common.DIR_CATALOGUES + 'catalogues-inv/LoLSS_DR2-inv-cat.fits',
                        help='Catalogue extracted from the inverted mosaics (default: %(default)s)')
    parser.add_argument('--outdir', default=common.DIR_ANALYSIS, help='Output directory (default: %(default)s)')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    t = Table.read(args.cat)
    t_fp = Table.read(args.cat_inv)
    print("Total fake sources:", len(t_fp))
    print("False positive: %.2f%%" % (100*len(t_fp)/len(t)))

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    ax.tick_params(direction='in', top=True, right=True)
    ax.set_xlabel(r'Flux density (mJy)')
    ax.set_ylabel(r'Number of sources')
    ax2.set_ylabel(r'Fraction of false positives')
    ax.set_yticklabels([])

    bins = np.arange(3, 600, 20)
    n_real, _, _ = ax.hist(t['Total_flux'], bins=bins, orientation='vertical', ec='black', histtype='bar', alpha=1, color='orange')
    n_false, _, _ = ax.hist(t_fp['Total_flux'], bins=bins, orientation='vertical', ec='black', histtype='bar', alpha=1, color='gray')
    centres = 0.5*(bins[1:] + bins[:-1])
    with np.errstate(divide='ignore', invalid='ignore'):
        frac = n_false/n_real
    ax2.plot(centres, frac, color='red')

    Table({'flux_lo': bins[:-1], 'flux_hi': bins[1:], 'n_real': n_real.astype(int), 'n_false': n_false.astype(int),
           'fraction': np.round(frac, 5)}).write(os.path.join(args.outdir, 'falsepositive.txt'),
                                                 format='ascii.fixed_width', delimiter=None, overwrite=True)

    ax.tick_params(axis=u'both', which=u'both', length=0)
    ax.grid(alpha=.2)
    ax.set_yscale('log')
    ax2.set_ylim(ymin=0, ymax=.05)
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red')
    handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in ['orange', 'gray']]
    ax.legend(handles, ["Detections", "False positives"])
    fig.savefig(os.path.join(args.outdir, 'falsepositive_hist.pdf'), bbox_inches='tight', facecolor='w')


if __name__ == '__main__':
    main()
