#!/usr/bin/env python3
"""
prepare_pointingsets4mosaic.py

For a list of primary beam FITS files, determine which fine HEALPix pixels
(nside=128 by default) are covered above a given beam threshold, then aggregate
by coarse HEALPix pixels (nside=16) and write the list of pointings that need
to be mosaiced together for each coarse pixel.

Usage:
    ~/scripts/survey/postprocess/prepare_pointingsets4mosaic.py ~/storage/surveytgts/done/*/primarybeam.fits
"""

import os
import sys
import glob
import argparse
import pickle
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy_healpix import HEALPix
from astropy.wcs import FITSFixedWarning
from scipy.ndimage import map_coordinates

from LiLF.surveys_db import SurveysDB

warnings.filterwarnings('ignore', category=FITSFixedWarning)


def get_beam_values(beam_data, wcs2d, sky_coords):
    """
    Bilinearly sample beam_data at sky positions given by sky_coords.

    Parameters
    ----------
    beam_data  : 2D ndarray
    wcs2d      : astropy.wcs.WCS (2D celestial)
    sky_coords : SkyCoord (N,)

    Returns
    -------
    values : ndarray (N,), 0.0 where outside the image or NaN in beam
    """
    ra  = sky_coords.ra.deg
    dec = sky_coords.dec.deg
    # wcs_world2pix returns 0-based pixel coords; origin=0
    xy  = wcs2d.wcs_world2pix(np.column_stack([ra, dec]), 0)
    x, y = xy[:, 0], xy[:, 1]

    # replace NaN in the beam (masked pixels) with 0 so they don't propagate
    beam_clean = np.where(np.isfinite(beam_data), beam_data, 0.0)

    # map_coordinates with mode='constant', cval=0 handles out-of-bounds coords
    values = map_coordinates(
        beam_clean,
        [y, x],
        order=1,
        mode='constant',
        cval=0.0,
    )
    return values


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Build per-coarse-HEALPix pointing lists from primary beam FITS files.'
        )
    )
    parser.add_argument(
        'beam_files', nargs='+',
        help='Primary beam FITS files (glob patterns are expanded automatically)',
    )
    parser.add_argument(
        '--threshold', type=float, default=0.4,
        help='Minimum beam value for a pixel to be considered covered (default: 0.4)',
    )
    parser.add_argument(
        '--search_radius', type=float, default=10.0,
        help='Search radius around pointing centre [deg] (default: 10.0)',
    )
    parser.add_argument(
        '--nside_fine', type=int, default=128,
        help='Fine HEALPix nside for per-pixel coverage map (default: 128)',
    )
    parser.add_argument(
        '--nside_coarse', type=int, default=16,
        help='Coarse HEALPix nside for output pointing sets (default: 16)',
    )
    parser.add_argument(
        '--output', default='/homes/fdg/storage/surveytgts/analysis/pointing_sets.txt',
        help='Output text file (default: /homes/fdg/storage/surveytgts/analysis/pointing_sets.txt)',
    )
    parser.add_argument(
        '--save_dict', default=None,
        help='If given, save the coarse-pixel->pointings dict as a pickle to this path',
    )
    args = parser.parse_args()

    # expand any glob patterns passed as arguments
    beam_files = []
    for pattern in args.beam_files:
        expanded = glob.glob(pattern)
        beam_files.extend(expanded if expanded else [pattern])
    beam_files = sorted(set(beam_files))

    if not beam_files:
        print('No beam files found.', file=sys.stderr)
        sys.exit(1)

    hp_fine   = HEALPix(nside=args.nside_fine,   order='ring', frame='icrs')
    hp_coarse = HEALPix(nside=args.nside_coarse, order='ring', frame='icrs')

    print(f'Fine   HEALPix: nside={args.nside_fine},   '
          f'{hp_fine.npix} pixels, '
          f'pixel size {hp_fine.pixel_resolution.to(u.arcmin):.1f}')
    print(f'Coarse HEALPix: nside={args.nside_coarse}, '
          f'{hp_coarse.npix} pixels, '
          f'pixel size {hp_coarse.pixel_resolution.to(u.deg):.2f}')
    print(f'Beam threshold: {args.threshold},  '
          f'search radius: {args.search_radius} deg\n')

    # overlay pointing status from the survey DB
    with SurveysDB(survey='lba', readonly=True) as sdb:
        sdb.execute('SELECT id, status FROM fields')
        rows = sdb.cur.fetchall()
    field_status = {row['id']: row['status'] for row in rows}

    # fine_pixel_index -> set of pointing names above threshold
    fine_to_pointings: dict[int, set] = {}
    # fine_pixel_index -> cumulative sum of beam values across all pointings
    fine_beam_sum: dict[int, float] = {}

    for beam_file in beam_files:
        # derive pointing name from the parent directory (matches LiLF convention)
        pointing_name = os.path.basename(os.path.dirname(os.path.abspath(beam_file)))

        status = field_status.get(pointing_name)
        if status != 'Done':
            print(f'Skipping {pointing_name} (status: {status})')
            continue

        print(f'Processing {pointing_name}  ({beam_file})')

        with fits.open(beam_file) as hdul:
            data   = hdul[0].data.copy()
            header = hdul[0].header.copy()

        # collapse any degenerate leading axes (Stokes, frequency, …)
        while data.ndim > 2:
            data = data[0]

        # 2D celestial WCS; wcs.crval gives the pointing centre (RA, Dec)
        wcs2d           = WCS(header).celestial
        ra_cen, dec_cen = wcs2d.wcs.crval
        centre          = SkyCoord(ra=ra_cen * u.deg, dec=dec_cen * u.deg, frame='icrs')

        # fine HEALPix pixels within the search cone around the pointing centre
        pixel_ids = hp_fine.cone_search_skycoord(centre, radius=args.search_radius * u.deg)
        if len(pixel_ids) == 0:
            print('  -> no pixels in cone, skipping')
            continue

        pixel_coords = hp_fine.healpix_to_skycoord(pixel_ids)
        beam_values  = get_beam_values(data, wcs2d, pixel_coords)

        covered_mask = beam_values >= args.threshold
        covered      = pixel_ids[covered_mask]
        covered_vals = beam_values[covered_mask]
        print(f'  -> {len(covered)}/{len(pixel_ids)} fine pixels above threshold')

        for pid, val in zip(covered, covered_vals):
            fine_to_pointings.setdefault(int(pid), set()).add(pointing_name)
            fine_beam_sum[int(pid)] = fine_beam_sum.get(int(pid), 0.0) + float(val)

    if not fine_to_pointings:
        print('\nNo pixels found above threshold. Check beam files and threshold value.',
              file=sys.stderr)
        sys.exit(1)

    print(f'\nFine pixels with coverage: {len(fine_to_pointings)}')

    # ------------------------------------------------------------------
    # aggregate fine-pixel coverage to coarse pixels:
    # each coarse pixel collects all pointings from any fine pixel inside it
    # ------------------------------------------------------------------
    coarse_to_pointings: dict[int, set] = {}

    fine_ids    = np.array(list(fine_to_pointings.keys()))
    fine_coords = hp_fine.healpix_to_skycoord(fine_ids)
    coarse_ids  = hp_coarse.skycoord_to_healpix(fine_coords)

    for coarse_id, fine_id in zip(coarse_ids, fine_ids):
        coarse_to_pointings.setdefault(int(coarse_id), set()).update(
            fine_to_pointings[int(fine_id)]
        )

    print(f'Coarse pixels covered: {len(coarse_to_pointings)}')

    # ------------------------------------------------------------------
    # write output: one line per coarse pixel
    # ------------------------------------------------------------------
    with open(args.output, 'w') as f:
        f.write(
            f'# coarse_healpix_id (nside={args.nside_coarse})'
            f'  pointing_names (beam >= {args.threshold})\n'
        )
        for coarse_id in sorted(coarse_to_pointings):
            pointings_str = ' '.join(sorted(coarse_to_pointings[coarse_id]))
            f.write(f'HP{coarse_id:04d}  {pointings_str}\n')

    print(f'Output written to {args.output}')

    if args.save_dict:
        with open(args.save_dict, 'wb') as f:
            pickle.dump(coarse_to_pointings, f)
        print(f'Coarse-pixel dict saved to {args.save_dict}')

    # ------------------------------------------------------------------
    # plots
    # ------------------------------------------------------------------
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    stem = os.path.splitext(args.output)[0]

    # --- 1. HEALPix Mollweide map of cumulative beam sum per fine pixel --
    beam_map = np.full(hp_fine.npix, np.nan)
    for pid, val in fine_beam_sum.items():
        beam_map[pid] = val

    all_ids    = np.where(~np.isnan(beam_map))[0]
    all_coords = hp_fine.healpix_to_skycoord(all_ids)
    all_ra_deg  = all_coords.ra.deg
    all_dec_deg = all_coords.dec.deg
    all_vals    = beam_map[all_ids]

    # keep only northern hemisphere
    north       = all_dec_deg >= 0
    all_ra_deg  = all_ra_deg[north]
    all_dec_deg = all_dec_deg[north]
    all_vals    = all_vals[north]

    # polar axes: theta = RA (radians, 0 at top, increasing clockwise to match
    # the astronomical convention of RA increasing to the left when north is up)
    # r = 90 - Dec  → Dec=90° at centre, Dec=0° at edge
    ra_rad  = np.deg2rad(all_ra_deg)
    r_polar = 90.0 - all_dec_deg           # co-declination [0, 90]

    fig = plt.figure(figsize=(8, 8))
    ax  = fig.add_subplot(111, projection='polar')
    ax.set_theta_zero_location('N')        # RA=0h at the top
    ax.set_theta_direction(-1)             # RA increases clockwise
    sc  = ax.scatter(
        ra_rad, r_polar,
        c=all_vals,
        s=0.3,
        cmap='viridis',
        rasterized=True,
    )
    ax.set_rlim(0, 90)                     # Dec from 90° (centre) to 0° (edge)
    # relabel radial ticks as Dec values
    rticks = [0, 15, 30, 45, 60, 75, 90]
    ax.set_yticks(rticks)
    ax.set_yticklabels([f'{90 - r}°' for r in rticks], fontsize=7)
    # label RA hours on the azimuthal axis
    ax.set_xticks(np.deg2rad([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]))
    ax.set_xticklabels(['0h','2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h'])
    ax.grid(True, alpha=0.4, linestyle='--', linewidth=0.5)

    # overlay pointing status from the survey DB
    with SurveysDB(survey='lba''', readonly=True) as sdb:
        sdb.execute('SELECT id, status, ra, decl FROM fields')
        rows = sdb.cur.fetchall()
    # rows are dicts: {'id': ..., 'status': ..., 'ra': ..., 'decl': ...}
    err_ra,  err_r  = [], []   # status == 'Error'  → red X
    warn_ra, warn_r = [], []   # anything else (not 'Done') → yellow X
    n_done = 0
    for row in rows:
        status  = row['status']
        ra_deg  = float(row['ra'])
        dec_deg = float(row['decl'])
        if dec_deg < 0:
            continue
        r_p = 90.0 - dec_deg
        ra_r = np.deg2rad(ra_deg)
        if status == 'Not observed':
            continue
        elif status == 'Error':
            err_ra.append(ra_r);  err_r.append(r_p)
        elif status == 'Done':
            n_done += 1
        else:
            warn_ra.append(ra_r); warn_r.append(r_p)
    handles = []
    if err_ra:
        h = ax.scatter(err_ra, err_r, marker='x', s=30, linewidths=1.5,
                        color='red', zorder=5, label=f'Error ({len(err_ra)})')
        handles.append(h)
    if warn_ra:
        h = ax.scatter(warn_ra, warn_r, marker='x', s=30, linewidths=1.5,
                        color='yellow', zorder=5,
                        label=f'Other / not done ({len(warn_ra)})')
        handles.append(h)
    if handles:
        pass  # legend rebuilt below after A-team markers
    print(f'DB overlay: {len(err_ra)} Error, {len(warn_ra)} other, {n_done} Done (not marked)')

    # bright A-team sources: Cassiopeia A and Cygnus A
    ateam = [
        ('Cas A', 350.866, 58.812),   # 23h 23m 27.9s  +58° 48' 42"
        ('Cyg A', 299.868, 40.734),   # 19h 59m 28.3s  +40° 44' 02"
    ]
    for name, ra_deg, dec_deg in ateam:
        h = ax.scatter(
            np.deg2rad(ra_deg), 90.0 - dec_deg,
            marker='D', s=70, color='black', zorder=6,
        )
        ax.annotate(
            name,
            xy=(np.deg2rad(ra_deg), 90.0 - dec_deg),
            xytext=(5, 5), textcoords='offset points',
            fontsize=7, color='black', zorder=6,
        )
    handles.append(
        ax.scatter([], [], marker='D', s=70, color='black', label='A-team')
    )
    ax.legend(handles=handles, loc='lower right', fontsize=8,
              bbox_to_anchor=(1.18, -0.05))

    fig.colorbar(sc, ax=ax, label='Sum of beam values', shrink=0.7, pad=0.08)
    ax.set_title(
        f'Fine HEALPix beam sum map — north polar view\n'
        f'(nside={args.nside_fine}, threshold={args.threshold})',
        pad=15,
    )
    fig.tight_layout()
    map_png = stem + f'_beam_map_thr{args.threshold}.png'
    fig.savefig(map_png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Beam map saved to {map_png}')

    # --- 2. Histogram of per-pixel beam sums (all sky) -------------------
    all_vals_full = beam_map[~np.isnan(beam_map)]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(all_vals_full, bins=50, color='steelblue', edgecolor='white', linewidth=0.4)
    ax.set_xlabel('Sum of beam values per fine pixel')
    ax.set_ylabel('Number of pixels')
    ax.set_title(
        f'Distribution of fine-pixel beam sums  (nside={args.nside_fine}, threshold={args.threshold})')
    ax.axvline(np.nanmedian(all_vals_full), color='orange', linestyle='--',
               label=f'Median = {np.nanmedian(all_vals_full):.2f}')
    ax.legend()
    fig.tight_layout()
    hist_png = stem + f'_beam_hist_thr{args.threshold}.png'
    fig.savefig(hist_png, dpi=150)
    plt.close(fig)
    print(f'Beam histogram saved to {hist_png}')


if __name__ == '__main__':
    main()

