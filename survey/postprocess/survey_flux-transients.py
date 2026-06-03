#!/usr/bin/env python3
"""
survey_flux-transients.py

For each pointing directory under ~/storage/surveytgts/done/, run PyBDSF
source extraction on:
  - wideDDS-TCXX-MFS-image.fits  (XX = 00-99, time-chunk images)
  - wideDDS-c0-MFS-image.fits    (DDS MFS continuum)

For every image the primary beam file is expected at the same path with
the suffix changed to '-pb.fits'  (e.g. wideDDS-TC00-MFS-image-pb.fits).

Usage:
    survey_flux-transients.py [--done-dir PATH] [--redo]
"""

import os
import re
import glob
import logging
import warnings
import argparse
import contextlib
from pathlib import Path
from typing import Optional

from astropy.wcs import FITSFixedWarning
warnings.filterwarnings('ignore', category=FITSFixedWarning)

import numpy as np
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy.io import fits as _fits
from astropy.wcs import WCS
import astropy.units as u
from scipy.ndimage import map_coordinates

import bdsf

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger(__name__)

# ── configuration ──────────────────────────────────────────────────────────────
DONE_DIR    = os.path.expanduser('~/storage/surveytgts/done-test/')
PB_SUFFIX   = '-pb.fits'   # appended to image stem to locate the PB image
BDSF_LOG    = os.path.expanduser('~/storage/surveytgts/analysis/bdsf.log')
CAT_DIR     = os.path.expanduser('~/storage/surveytgts/analysis/catalogues/')

# glob patterns used to discover images inside each pointing directory
IMAGE_PATTERNS = [
    'wideDDS-TC[0-9][0-9]-c0-MFS-image.fits',   # time-chunk images
    'wideDDS-c0-MFS-image.fits',              # DDS MFS continuum
]

# PyBDSF keyword arguments (tune as needed)
BDSF_KWARGS = dict(
    thresh_isl        = 4.0,
    thresh_pix        = 5.0,
    rms_box           = (150, 15),
    rms_map           = True,
    mean_map          = 'zero',
    ini_method        = 'intensity',
    adaptive_rms_box  = True,
    adaptive_thresh   = 50,
    rms_box_bright    = (60, 15),
    group_by_isl      = False,
    group_tol         = 10.0,
    output_opts       = True,
    output_all        = False,
    atrous_do         = False,
    atrous_jmax       = 4,
    flagging_opts     = True,
    flag_maxsize_fwhm = 0.5,
    advanced_opts     = True,
    blank_limit       = None,
    quiet             = True,
    detection_image    = None,  # to be set per run
)


# ── helpers ────────────────────────────────────────────────────────────────────

def pb_path_for(image: Path) -> Path:
    """Primary-beam FITS expected alongside the image."""
    return image.with_name(image.stem + PB_SUFFIX)

def catalog_path_for(image: Path) -> Path:
    """Output source-list catalogue: CAT_DIR/pointingname_imagename.fits"""
    return Path(CAT_DIR) / f'{image.parent.name}_{image.stem}.fits'


def filter_catalog_by_pb(cat_path: Path, pb_path: Path, threshold: float):
    """Read cat_path, drop sources where primarybeam < threshold, overwrite."""
    if not pb_path.exists():
        log.warning('    primarybeam not found: %s — skipping PB filter', pb_path)
        return
    try:
        t = Table.read(cat_path)
        with _fits.open(pb_path) as hdul:
            data  = np.squeeze(hdul[0].data).astype(float)
            wcs2d = WCS(hdul[0].header).celestial
        sky  = SkyCoord(ra=np.array(t['RA']) * u.deg,
                        dec=np.array(t['DEC']) * u.deg)
        px, py = wcs2d.world_to_pixel(sky)
        beam_vals = map_coordinates(data, [py, px], order=1,
                                    mode='constant', cval=0.0)
        keep = beam_vals >= threshold
        log.info('    PB filter: kept %d/%d sources (thr=%.2f) in %s',
                 keep.sum(), len(t), threshold, cat_path.name)
        t[keep].write(cat_path, overwrite=True)
    except Exception as exc:
        log.warning('    PB filter failed for %s: %s', cat_path.name, exc)

# ── classes ───────────────────────────────────────────────────────────────────

class BDSFRun:
    """Wraps a single PyBDSF source-extraction run for one image."""

    def __init__(self, image: Path, image_pb: Path):
        self.image    = image
        self.image_pb = image_pb
        self.catalog  = catalog_path_for(image)

    def run(self, redo: bool = False, pb_threshold: float = 0.4) -> bool:
        """Run PyBDSF. Returns True on success. Skips if catalogue already exists."""
        if not redo and self.catalog.exists():
            log.info('    catalogue exists, skipping: %s', self.catalog.name)
            return True

        if not self.image.exists():
            log.warning('    image not found: %s', self.image)
            return False

        kwargs = dict(BDSF_KWARGS)
        kwargs['detection_image'] = str(self.image)

        log.info('    bdsf: %s', self.image_pb.name)
        try:
            os.makedirs(os.path.dirname(BDSF_LOG), exist_ok=True)
            os.makedirs(CAT_DIR, exist_ok=True)
            _root_logger   = logging.getLogger()
            _handlers_before = list(_root_logger.handlers)
            with open(BDSF_LOG, 'a') as _bdsf_logf:
                with contextlib.redirect_stdout(_bdsf_logf), contextlib.redirect_stderr(_bdsf_logf):
                    img = bdsf.process_image(str(self.image_pb), **kwargs)
            img.write_catalog(
                outfile      = str(self.catalog),
                catalog_type = 'srl',
                format       = 'fits',
                clobber      = True,
            )
            img.write_catalog(
                outfile      = str(self.catalog.with_suffix('').with_suffix('.gaul.fits')),
                catalog_type = 'gaul',
                format       = 'fits',
                clobber      = True,
            )
            # remove any handlers bdsf added that now point to a closed file
            for _h in list(_root_logger.handlers):
                if _h not in _handlers_before:
                    _root_logger.removeHandler(_h)
                    _h.close()
            # filter both catalogues by the pointing primary beam
            _pb = self.image.parent / 'primarybeam.fits'
            filter_catalog_by_pb(self.catalog, _pb, pb_threshold)
            _gaul = self.catalog.with_suffix('').with_suffix('.gaul.fits')
            filter_catalog_by_pb(_gaul, _pb, pb_threshold)
            return True
        except Exception as exc:
            log.error('    bdsf FAILED on %s: %s', self.image.name, exc)
            return False


class Pointing:
    """
    Represents one survey pointing directory.
    Discovers all relevant images and drives PyBDSF runs.
    """

    def __init__(self, directory: Path):
        self.dir  = directory
        self.name = directory.name
        self.runs: list[BDSFRun] = []
        self._discover()

    def _discover(self):
        for pattern in IMAGE_PATTERNS:
            for img in sorted(self.dir.glob(pattern)):
                self.runs.append(BDSFRun(img, pb_path_for(img)))

    def process(self, redo: bool = False, pb_threshold: float = 0.4):
        log.info('Pointing %-20s  %d image(s)', self.name, len(self.runs))
        for run in self.runs:
            run.run(redo=redo, pb_threshold=pb_threshold)

    def catalogs(self) -> list[Path]:
        return [r.catalog for r in self.runs if r.catalog.exists()]


class SurveyProcessor:
    """Top-level controller: iterates all pointing directories."""

    def __init__(self, done_dir: str = DONE_DIR):
        self.done_dir  = Path(done_dir)
        self.pointings: list[Pointing] = []
        self._discover()

    def _discover(self):
        subdirs = sorted(p for p in self.done_dir.iterdir() if p.is_dir())
        self.pointings = [Pointing(d) for d in subdirs]
        log.info('Found %d pointing directory/ies in %s',
                 len(self.pointings), self.done_dir)

    def process_all(self, redo: bool = False, pb_threshold: float = 0.4):
        for pointing in self.pointings:
            pointing.process(redo=redo, pb_threshold=pb_threshold)

    def all_catalogs(self) -> list[Path]:
        cats = []
        for p in self.pointings:
            cats.extend(p.catalogs())
        return cats

# ── catalogue collection ───────────────────────────────────────────────────────

def collect_catalogs(cat_paths: list[Path]) -> dict[str, dict]:
    """
    Load all source-list catalogues and split them by image type.

    Returns a dict with keys:
      'tc'  : dict[pointing_name -> list[Table]]  (one Table per TC epoch)
      'dds' : dict[pointing_name -> Table]  (one DDS MFS catalogue per pointing)
    """
    tc_by_pointing:  dict[str, list[Table]] = {}
    dds_by_pointing: dict[str, Table]       = {}

    for cat_path in cat_paths:
        try:
            t = Table.read(cat_path)
        except Exception as exc:
            log.warning('Cannot read catalogue %s: %s', cat_path, exc)
            continue

        # tag each source with its origin
        # filename format: {pointing}_{image_stem}.fits — all cats live in CAT_DIR
        pointing = cat_path.stem.split('_wideDDS')[0]
        name     = cat_path.name
        _tc_match = re.search(r'wideDDS-TC(\d{2})', name)
        t['cat_type']    = f'TC{_tc_match.group(1)}' if _tc_match else 'DDS'
        t['pointing']    = pointing

        if _tc_match:
            tc_by_pointing.setdefault(pointing, []).append(t)
        elif 'wideDDS-c0' in name:
            if pointing in dds_by_pointing:
                log.warning('Duplicate DDS catalogue (%s) for pointing %s, overwriting', cat_path, pointing)
            dds_by_pointing[pointing] = t

    # keep TC tables per pointing as a list (one Table per epoch, not stacked)
    n_tc_src    = sum(len(t) for tlist in tc_by_pointing.values() for t in tlist)
    n_tc_epochs = sum(len(tlist) for tlist in tc_by_pointing.values())
    log.info('Collected %6d sources across %d epoch(s) in %d pointing(s)  [tc]',
             n_tc_src, n_tc_epochs, len(tc_by_pointing))

    n_dds = sum(len(v) for v in dds_by_pointing.values())
    log.info('Collected %6d sources across %d pointing(s)  [dds]',
             n_dds, len(dds_by_pointing))

    return {'tc': tc_by_pointing, 'dds': dds_by_pointing}

# ── main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--done-dir', default=DONE_DIR,
                        help='Path to the done/ directory (default: %(default)s)')
    parser.add_argument('--redo', action='store_true',
                        help='Re-run PyBDSF even if a catalogue already exists')
    parser.add_argument('--pb-threshold', type=float, default=0.4,
                        help='Primary beam threshold for source filtering (default: 0.4)')
    args = parser.parse_args()

    # ── source extraction ──────────────────────────────────────────────────────
    processor = SurveyProcessor(done_dir=args.done_dir)
    processor.process_all(redo=args.redo, pb_threshold=args.pb_threshold)

    # ── collect catalogues ─────────────────────────────────────────────────────
    all_cat_paths = processor.all_catalogs()
    log.info('\nTotal catalogues found: %d', len(all_cat_paths))
    cats = collect_catalogs(all_cat_paths)

    cat_tc  = cats['tc']    # dict[pointing_name -> list[Table]]  (one Table per TC epoch)
    cat_dds = cats['dds']   # dict[pointing_name -> Table]  (one DDS MFS per pointing)

    # ── TC quality filter ─────────────────────────────────────────────────────
    # Reject TC epochs whose median Isl_rms deviates more than 5σ from the
    # global median across all TC epochs of all fields.  σ = 1.4826 * MAD.
    _tc_medians: list[tuple[str, int, float]] = []
    for _pname, _tc_list in cat_tc.items():
        for _ei, _tc in enumerate(_tc_list):
            _tc_medians.append((_pname, _ei, float(np.median(_tc['Isl_rms']))))

    if _tc_medians:
        _rms_vals   = np.array([m[2] for m in _tc_medians])
        _global_med = np.median(_rms_vals)
        _sigma      = 1.4826 * np.median(np.abs(_rms_vals - _global_med))
        _threshold  = 4.0 * _sigma
        log.info('TC Isl_rms filter: global median=%.3e  MAD-sigma=%.3e  4sigma threshold=%.3e',
                 _global_med, _sigma, _threshold)

        _keep = {(_pname, _ei) for _pname, _ei, _med in _tc_medians
                 if abs(_med - _global_med) <= _threshold}
        _n_before = sum(len(v) for v in cat_tc.values())
        cat_tc = {
            _pname: [_tc for _ei, _tc in enumerate(_tc_list) if (_pname, _ei) in _keep]
            for _pname, _tc_list in cat_tc.items()
        }
        cat_tc = {k: v for k, v in cat_tc.items() if v}
        _n_after = sum(len(v) for v in cat_tc.values())
        log.info('TC quality filter: kept %d/%d epochs', _n_after, _n_before)

    # ── cross-match section ────────────────────────────────────────────────────
    # PyBDSF source-list column names assumed: RA, DEC, Total_flux, E_Total_flux,
    # Peak_flux, E_Peak_flux.  Adjust if your PyBDSF version differs.

    # load pointing grid (name, ra[deg], dec[deg])
    _grid_path = os.path.expanduser('~/storage/allsky-grid.fits')
    _grid_tbl  = Table.read(_grid_path)
    grid: dict[str, tuple[float, float]] = {
        str(row['name']).strip(): (float(row['ra']), float(row['dec']))
        for row in _grid_tbl
    }
    log.info('Loaded %d pointings from %s', len(grid), _grid_path)

    MATCH_RADIUS   = 5  * u.arcsec
    MIN_DETECTIONS = 2            # source must appear in at least this many TC epochs
    SN_MIN         = 5.0          # S/N cut for transient candidates
    RATIO_HIGH     = 2.0          # upper flux-ratio threshold for transients
    RATIO_LOW      = 0.5          # lower flux-ratio threshold for transients
    RADIAL_STEP    = 0.5          # deg — radial bin width for dispersion plot
    FLUX_BINS      = [0, 1e-3, 3e-3, 1e-2, 1e-1, np.inf]  # Jy edges

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    all_transients: list[Table] = []
    # dispersion data collected across all pointings for the final plot
    # each entry: (r_mid_deg, flux_bin_idx, dispersion_value)
    disp_records: list[tuple[float, int, float]] = []

    for pointing, tc_list in cat_tc.items():
        if len(tc_list) < 2:
            log.info('Pointing %s: only %d TC epoch(s), skipping cross-match',
                     pointing, len(tc_list))
            continue
        if pointing not in cat_dds:
            log.warning('Pointing %s: no DDS catalogue, skipping', pointing)
            continue

        dds  = cat_dds[pointing]
        dds_coords = SkyCoord(ra=dds['RA'], dec=dds['DEC'])

        # phase centre from the survey grid
        if pointing[:-1] in grid:
            _ra0, _dec0 = grid[pointing[:-1]]  # remove trailing 'o' or 's' if present
        else:
            log.warning('Pointing %s not found in allsky-grid.fits, using DDS median', pointing)
            _ra0, _dec0 = float(np.median(dds['RA'])), float(np.median(dds['DEC']))
        phase_centre = SkyCoord(ra=_ra0* u.deg, dec=_dec0* u.deg)

        # --- match every TC epoch to the DDS catalogue ---
        # ratio_table columns: dds_idx, tc_epoch, ratio, sn, r_deg
        per_source_ratios: dict[int, list[float]] = {}   # dds_row_idx -> [ratio, ...]
        per_source_sn:     dict[int, float]       = {}   # dds_row_idx -> DDS S/N
        per_source_r:      dict[int, float]       = {}   # dds_row_idx -> dist from centre [deg]

        for epoch_idx, tc in enumerate(tc_list):
            tc_coords = SkyCoord(ra=tc['RA'], dec=tc['DEC'])
            idx, sep, _ = tc_coords.match_to_catalog_sky(dds_coords)
            matched = sep < MATCH_RADIUS

            for tc_row, dds_row in zip(np.where(matched)[0], idx[matched]):
                ratio = float(tc['Total_flux'][tc_row]) / float(dds['Total_flux'][dds_row])
                per_source_ratios.setdefault(int(dds_row), []).append(ratio)
                if int(dds_row) not in per_source_sn:
                    sn  = float(dds['Peak_flux'][dds_row]) / float(dds['Isl_rms'][dds_row])
                    r   = dds_coords[dds_row].separation(phase_centre).deg
                    per_source_sn[int(dds_row)]  = sn
                    per_source_r[int(dds_row)]   = r

        # --- filter: detected in at least MIN_DETECTIONS epochs ---
        multi = {k: v for k, v in per_source_ratios.items() if len(v) >= MIN_DETECTIONS}
        log.info('Pointing %s: %d/%d DDS sources detected in %d TC epochs',
                 pointing, len(multi), len(dds), len(tc_list))

        if not multi:
            continue

        # --- average ratio per source ---
        dds_rows   = np.array(sorted(multi.keys()))
        avg_ratio  = np.array([np.mean(multi[k]) for k in dds_rows])
        sn_vals    = np.array([per_source_sn[k]  for k in dds_rows])
        r_vals     = np.array([per_source_r[k]   for k in dds_rows])

        # --- transient candidates ---
        transient_mask = (sn_vals > SN_MIN) & ((avg_ratio > RATIO_HIGH) | (avg_ratio < RATIO_LOW))
        if transient_mask.any():
            t_rows = dds[dds_rows[transient_mask]].copy()
            t_rows['avg_ratio']  = avg_ratio[transient_mask]
            t_rows['sn']         = sn_vals[transient_mask]
            t_rows['r_deg']      = r_vals[transient_mask]
            t_rows['pointing']   = pointing
            all_transients.append(t_rows)
            log.info('  -> %d transient candidate(s)', transient_mask.sum())

        # --- dispersion per (radial bin, flux bin) ---
        r_max = r_vals.max()
        r_edges = np.arange(0, r_max + RADIAL_STEP, RADIAL_STEP)
        flux_vals = dds['Total_flux'][dds_rows]

        for ri in range(len(r_edges) - 1):
            r_lo, r_hi = r_edges[ri], r_edges[ri + 1]
            r_mid = 0.5 * (r_lo + r_hi)
            in_r = (r_vals >= r_lo) & (r_vals < r_hi)
            if not in_r.any():
                continue
            for fi in range(len(FLUX_BINS) - 1):
                f_lo, f_hi = FLUX_BINS[fi], FLUX_BINS[fi + 1]
                in_f = (flux_vals[in_r] >= f_lo) & (flux_vals[in_r] < f_hi)
                if in_f.sum() < 2:
                    continue
                disp = float(np.std(avg_ratio[in_r][in_f]))
                disp_records.append((r_mid, fi, disp))

    # --- save transient catalogue ---
    if all_transients:
        cat_transients = vstack(all_transients)
        cat_transients.write('/homes/fdg/storage/surveytgts/analysis/transients/transients.fits', overwrite=True)
        log.info('Transient catalogue saved: transients.fits  (%d sources)', len(cat_transients))
    else:
        log.info('No transient candidates found.')

    # ── plot section ───────────────────────────────────────────────────────────
    if disp_records:
        flux_labels = []
        for i in range(len(FLUX_BINS) - 1):
            lo = FLUX_BINS[i]   * 1e3   # convert to mJy for label
            hi = FLUX_BINS[i+1] * 1e3
            flux_labels.append(f'{lo:.0f}–{hi:.0f} mJy' if hi < 1e6 else f'>{lo:.0f} mJy')

        n_fbins = len(FLUX_BINS) - 1
        cmap    = cm.get_cmap('tab10', n_fbins)

        fig, ax = plt.subplots(figsize=(9, 6))
        for fi in range(n_fbins):
            pts = [(r, d) for r, fbi, d in disp_records if fbi == fi]
            if not pts:
                continue
            r_arr = np.array([p[0] for p in pts])
            d_arr = np.array([p[1] for p in pts])
            ax.scatter(r_arr, d_arr, s=18, alpha=0.7,
                       color=cmap(fi), label=flux_labels[fi])

        ax.set_xlabel('Radial distance from phase centre [deg]')
        ax.set_ylabel('Std dev of TC/DDS flux ratio')
        ax.set_title('Flux-ratio dispersion vs. radial distance\n'
                     '(one point per radial-bin per field per flux bin)')
        ax.legend(title='DDS flux bin', fontsize=8, title_fontsize=8,
                  loc='upper left', bbox_to_anchor=(1.01, 1))
        fig.tight_layout()
        fig.savefig('ratio_dispersion.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        log.info('Dispersion plot saved: ratio_dispersion.png')
    else:
        log.info('No dispersion data to plot.')


if __name__ == '__main__':
    main()

