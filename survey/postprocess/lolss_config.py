"""
Shared configuration for the LoLSS post-processing scripts.

Paths, HEALPix tiling, beam threshold and PyBDSF parameters used to be repeated
in each script; they are defined here once so that mosaics, catalogues,
completeness and false-positive runs stay consistent with each other.

The scripts are run directly from this directory, so they can simply do
    import lolss_config as cfg
The survey root can be changed without editing this file by setting the
LOLSS_DIR environment variable.
"""

import os

# ── paths ──────────────────────────────────────────────────────────────────────
DIR_SURVEY         = os.path.join(os.environ.get('LOLSS_DIR', '/homes/fdg/storage/surveytgts/'), '')
DIR_DONE           = DIR_SURVEY + 'done/'             # one directory per processed pointing
DIR_MOSAICS        = DIR_SURVEY + 'mosaics/healpix/'  # HPxxxx-mosaicI.fits
DIR_CATALOGUES     = DIR_SURVEY + 'catalogues/'       # per-mosaic PyBDSF products and concatenated catalogues
DIR_ANALYSIS       = DIR_SURVEY + 'analysis/'
FILE_POINTING_SETS = DIR_ANALYSIS + 'pointing_sets.txt'
FILE_GRID          = os.path.expanduser('~/storage/allsky-grid.fits')

# file names inside each pointing directory
FILE_IMAGE = 'wideDDS-c0-MFS-image.fits'
FILE_BEAM  = 'primarybeam.fits'

# ── survey / mosaic geometry ───────────────────────────────────────────────────
HEALPIX_NSIDE      = 16    # coarse tiling: one mosaic and one catalogue per pixel
HEALPIX_NSIDE_FINE = 512   # fine tiling used to decide which pointings cover a coarse pixel
BEAM_THRESHOLD     = 0.4   # primary-beam level below which a pointing is not used
RESTFRQ            = 54e6  # Hz
RESOLUTION         = 15    # arcsec, FWHM of the restoring beam of pointings and mosaics
CELLSIZE           = 4     # arcsec, mosaic pixel size

# ── source extraction ──────────────────────────────────────────────────────────
# PyBDSF parameters of the survey catalogue. Use bdsf_kwargs() to get a copy,
# optionally with run-specific additions (frequency, rmsmean_map_filename, ...).
BDSF_KWARGS = dict(
    thresh_isl        = 4.0,
    thresh_pix        = 5.0,
    rms_box           = (150, 15),
    rms_map           = True,
    mean_map          = 'zero',
    ini_method        = 'intensity',
    adaptive_rms_box  = True,
    adaptive_thresh   = 50,
    rms_box_bright    = (50, 10),
    group_by_isl      = False,
    group_tol         = 10.0,
    output_opts       = True,
    output_all        = False,
    atrous_do         = False,
    flagging_opts     = True,
    flag_maxsize_fwhm = 0.5,
    advanced_opts     = True,
    blank_limit       = None,
)


def bdsf_kwargs(**overrides):
    """Return a copy of BDSF_KWARGS updated with the given keyword arguments."""
    kwargs = dict(BDSF_KWARGS)
    kwargs.update(overrides)
    return kwargs
