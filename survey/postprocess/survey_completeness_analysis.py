#!/usr/bin/env python3
"""
survey_completeness_analysis.py

Analyse the injection files written by survey_completeness.py (one
HPxxxx-mosaicI.injection.txt per residual mosaic) and derive the survey
completeness:

  - fraction of injected sources detected as a function of flux density
    (differential) and above a flux density limit (integral);
  - flux densities at which the completeness reaches the requested levels,
    interpolated from the curves;
  - the same as a function of the local signal-to-noise ratio, using the rms
    maps of the catalogue, and of flux density in bins of local rms.

Only the sources injected in the HEALPix pixel of their own mosaic are used,
as it is done for the catalogue.

The integral completeness depends on the flux distribution of the injected
sources (the -I slope of survey_completeness.py); the differential one does not.

Outputs (in --outdir):
  completeness_injections.fits   all the injected sources, with local rms
  completeness_curve.txt         differential and integral completeness
  completeness_curve.png         completeness vs flux density
  completeness_snr.png           completeness vs local S/N and per rms bin

Usage:
    survey_completeness_analysis.py [--injdir DIR] [--rmsdir DIR] [--outdir DIR]
"""

import os
import glob
import argparse
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS, FITSFixedWarning
from astropy_healpix import HEALPix
import astropy.units as u

warnings.filterwarnings('ignore', category=FITSFixedWarning)


def read_injections(injdir, nside):
    """Read all the injection files, keeping the sources in the pixel of their own mosaic."""
    hp = HEALPix(nside=nside)
    tables = []
    for datafile in sorted(glob.glob(os.path.join(injdir, 'HP*injection.txt'))):
        pix = int(os.path.basename(datafile)[2:6])
        t = Table.read(datafile, format='ascii')
        if len(t) == 0:
            continue
        t['HEALPIX'] = hp.lonlat_to_healpix(t['ra'].data*u.deg, t['dec'].data*u.deg)
        t = t[t['HEALPIX'] == pix]
        t.meta = {}
        tables.append(t)
    if not tables:
        raise SystemExit(f'No injection file with sources found in {injdir}')
    t = vstack(tables)
    t.rename_columns(['fv', 'rfv', 'efv'], ['f_inj', 'f_meas', 'f_meas_err'])
    for col in ('f_inj', 'f_meas', 'f_meas_err'):
        t[col] *= 1e3
        t[col].unit = u.mJy
    t['detected'] = t['f_meas'] != 0
    return t


def add_local_rms(t, rmsdir):
    """Add the value of the catalogue rms map (mJy/beam) at the position of each injected source."""
    t['rms'] = np.nan
    t['rms'].unit = u.mJy/u.beam
    missing = 0
    for pix in np.unique(t['HEALPIX']):
        rmsfile = os.path.join(rmsdir, f'HP{pix:04d}-mosaicI.rms.fits')
        if not os.path.exists(rmsfile):
            missing += 1
            continue
        sel = np.where(t['HEALPIX'] == pix)[0]
        with fits.open(rmsfile, memmap=True) as hdul:
            wcs = WCS(hdul[0].header).celestial
            ny, nx = hdul[0].data.shape[-2:]
            x, y = wcs.all_world2pix(t['ra'][sel], t['dec'][sel], 0)
            x, y = np.rint(x).astype(int), np.rint(y).astype(int)
            inside = (x >= 0) & (x < nx) & (y >= 0) & (y < ny)
            # read one pixel at a time: the map stays on disk
            data = hdul[0].data.reshape(ny, nx)
            t['rms'][sel[inside]] = [float(data[j, i])*1e3 for i, j in zip(x[inside], y[inside])]
    if missing:
        print(f'WARNING: rms map missing for {missing} pixels')
    return t


def binomial_error(k, n):
    """Half-width of the 68% Wilson interval for k successes out of n."""
    k, n = np.asarray(k, dtype=float), np.asarray(n, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        p = k/n
        err = np.sqrt(p*(1-p)/n + 0.25/n**2) / (1 + 1/n)
    return err


def fraction_in_bins(x, detected, edges):
    """Detected fraction (with error) of the sources in each bin of x."""
    n_inj = np.histogram(x, bins=edges)[0]
    n_det = np.histogram(x[detected], bins=edges)[0]
    with np.errstate(divide='ignore', invalid='ignore'):
        frac = n_det / n_inj
    return frac, binomial_error(n_det, n_inj), n_inj


def fraction_above(x, detected, limits):
    """Detected fraction (with error) of the sources with x above each limit."""
    order = np.argsort(x)
    xs = x[order]
    # number of injected / detected sources with x > limit
    n_inj = len(xs) - np.searchsorted(xs, limits, side='right')
    cum_det = np.concatenate([[0], np.cumsum(detected[order])])
    n_det = cum_det[-1] - cum_det[len(xs) - n_inj]
    with np.errstate(divide='ignore', invalid='ignore'):
        frac = n_det / n_inj
    return frac, binomial_error(n_det, n_inj), n_inj


def level_crossing(x, frac, level):
    """
    Smallest x above which the curve stays at or above level, interpolating
    linearly in log(x). NaN if the level is never reached or always exceeded.
    """
    ok = np.isfinite(frac)
    x, frac = np.asarray(x)[ok], np.asarray(frac)[ok]
    below = np.where(frac < level)[0]
    if len(below) == 0 or below[-1] == len(x) - 1:
        return np.nan
    i = below[-1]
    w = (level - frac[i]) / (frac[i+1] - frac[i])
    return 10**(np.log10(x[i]) + w*(np.log10(x[i+1]) - np.log10(x[i])))


def main():
    survey = os.path.expanduser('~/storage/surveytgts/')
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--injdir', default=survey+'mosaics/healpix-res/',
                        help='Directory with the HP*injection.txt files (default: %(default)s)')
    parser.add_argument('--rmsdir', default=survey+'catalogues/',
                        help='Directory with the HP*-mosaicI.rms.fits maps (default: %(default)s)')
    parser.add_argument('--outdir', default=survey+'analysis/',
                        help='Output directory (default: %(default)s)')
    parser.add_argument('--nside', type=int, default=16,
                        help='HEALPix nside of the mosaics (default: 16)')
    parser.add_argument('--levels', type=float, nargs='+', default=[0.5, 0.9, 0.95],
                        help='Completeness levels to report (default: 0.5 0.9 0.95)')
    parser.add_argument('--fmax', type=float, default=100.,
                        help='Maximum flux density of the curves [mJy] (default: 100)')
    parser.add_argument('--nbins', type=int, default=25,
                        help='Number of logarithmic bins of the differential curves (default: 25)')
    parser.add_argument('--no-rms', action='store_true',
                        help='Do not read the rms maps (skip the S/N analysis)')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    t = read_injections(args.injdir, args.nside)
    n_cycles = len(np.unique(t['cycle']))
    print(f'{len(t)} injected sources in their own pixel, {np.sum(t["detected"])} detected, '
          f'{len(np.unique(t["HEALPIX"]))} pixels, {n_cycles} cycle(s)')
    if not args.no_rms:
        t = add_local_rms(t, args.rmsdir)
    t.write(os.path.join(args.outdir, 'completeness_injections.fits'), overwrite=True)

    f_inj = np.array(t['f_inj'])
    detected = np.array(t['detected'])
    fmin = f_inj.min()

    # ── completeness vs flux density ───────────────────────────────────────────
    edges = np.logspace(np.log10(fmin), np.log10(args.fmax), args.nbins + 1)
    centres = np.sqrt(edges[1:]*edges[:-1])
    diff, diff_err, diff_n = fraction_in_bins(f_inj, detected, edges)
    limits = np.logspace(np.log10(fmin), np.log10(args.fmax), 200)
    integ, integ_err, integ_n = fraction_above(f_inj, detected, limits)

    curve = Table({'flux_lo': edges[:-1], 'flux_hi': edges[1:], 'flux': centres,
                   'n_inj': diff_n, 'frac_diff': diff, 'frac_diff_err': diff_err,
                   'frac_above_flux_lo': np.interp(edges[:-1], limits, integ)})
    for col in curve.colnames:
        if curve[col].dtype.kind == 'f':
            curve[col].format = '.4g'
    curve.write(os.path.join(args.outdir, 'completeness_curve.txt'), format='ascii.fixed_width',
                delimiter=None, overwrite=True)

    print('\nCompleteness levels [mJy]')
    print('  level   integral (>S)   differential (at S)')
    crossings = {}
    for level in args.levels:
        crossings[level] = (level_crossing(limits, integ, level), level_crossing(centres, diff, level))
        print(f'  {level:5.0%}   {crossings[level][0]:13.2f}   {crossings[level][1]:19.2f}')
    print(f'(integral values assume the injected flux distribution; injected range {fmin:.2f}-{f_inj.max():.0f} mJy)')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.fill_between(limits, integ - integ_err, integ + integ_err, color='C3', alpha=0.2, lw=0)
    ax.plot(limits, integ, color='C3', label='above flux density limit')
    ax.errorbar(centres, diff, yerr=diff_err, color='C0', marker='o', ms=3, lw=1,
                label='at a given flux density')
    for level in args.levels:
        s_int = crossings[level][0]
        if np.isfinite(s_int):
            ax.plot([s_int, s_int], [0, level], color='k', ls=':', lw=0.8)
            ax.plot([fmin, s_int], [level, level], color='k', ls=':', lw=0.8)
            ax.annotate(f'{level:.0%}: {s_int:.1f} mJy', xy=(fmin, level), xytext=(4, 2),
                        textcoords='offset points', fontsize=7)
    ax.set_xscale('log')
    ax.set_xlim(fmin, args.fmax)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel('Integrated flux density (mJy)')
    ax.set_ylabel('Completeness')
    ax.tick_params(direction='in', top=True, right=True, which='both')
    ax.grid(alpha=0.2)
    ax.legend(loc='lower right', fontsize=8, frameon=False)
    fig.savefig(os.path.join(args.outdir, 'completeness_curve.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)

    # ── completeness vs local S/N and in bins of local rms ─────────────────────
    if args.no_rms:
        return
    rms = np.array(t['rms'])
    good = np.isfinite(rms) & (rms > 0)
    print(f'\nLocal rms available for {good.sum()}/{len(t)} sources: median {np.median(rms[good]):.2f} mJy/beam')
    snr = f_inj[good]/rms[good]
    det = detected[good]

    snr_edges = np.logspace(np.log10(max(snr.min(), 1.)), np.log10(100.), args.nbins + 1)
    snr_centres = np.sqrt(snr_edges[1:]*snr_edges[:-1])
    snr_frac, snr_err, _ = fraction_in_bins(snr, det, snr_edges)
    print('Completeness levels in local S/N (injected flux density / local rms)')
    for level in args.levels:
        print(f'  {level:5.0%}   {level_crossing(snr_centres, snr_frac, level):6.1f}')

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(11, 4))
    ax1.errorbar(snr_centres, snr_frac, yerr=snr_err, color='C0', marker='o', ms=3, lw=1)
    ax1.set_xscale('log')
    ax1.set_xlabel('Injected flux density / local rms')
    ax1.set_ylabel('Fraction detected')

    # same sources split in three bins of local rms
    q = np.quantile(rms[good], [0, 1/3, 2/3, 1])
    for i, color in enumerate(['C2', 'C1', 'C3']):
        sel = (rms[good] >= q[i]) & (rms[good] <= q[i+1])
        frac, err, _ = fraction_in_bins(f_inj[good][sel], det[sel], edges)
        ax2.errorbar(centres, frac, yerr=err, color=color, marker='o', ms=3, lw=1,
                     label=f'rms {q[i]:.2f}-{q[i+1]:.2f} mJy/beam')
        print(f'  rms {q[i]:.2f}-{q[i+1]:.2f} mJy/beam: 50% at {level_crossing(centres, frac, 0.5):.2f} mJy, '
              f'95% at {level_crossing(centres, frac, 0.95):.2f} mJy (differential)')
    ax2.set_xscale('log')
    ax2.set_xlim(fmin, args.fmax)
    ax2.set_xlabel('Integrated flux density (mJy)')
    ax2.set_ylabel('Fraction detected')
    ax2.legend(loc='lower right', fontsize=8, frameon=False)
    for ax in (ax1, ax2):
        ax.set_ylim(-0.02, 1.05)
        ax.tick_params(direction='in', top=True, right=True, which='both')
        ax.grid(alpha=0.2)
    fig.savefig(os.path.join(args.outdir, 'completeness_snr.png'), dpi=200, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
