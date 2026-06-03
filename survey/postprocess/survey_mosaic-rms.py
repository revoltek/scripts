#!/usr/bin/env python3

import os, sys, glob, logging, subprocess, math, re, warnings
import numpy as np
from astropy.io import fits as afits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy_healpix import HEALPix

warnings.filterwarnings('ignore', category=FITSFixedWarning)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

os.chdir(os.path.expanduser('~/storage/surveytgts/'))

dir_raw     = 'catalogues'  # input: *.rms.fits from survey_make_catalogues.py
dir_proj    = 'mosaics-rms/projected'
dir_blank   = 'mosaics-rms/blanked'
hdr_final   = 'mosaics-rms/mosaic-rms.hdr'
tbl_proj    = 'mosaics-rms/mosaic-rms.tbl'
fits_output = 'mosaics-rms/mosaic-rms.fits'
healpix_nside = 16
hp_grid = HEALPix(nside=healpix_nside, order='ring', frame='icrs')

if not os.path.exists(dir_raw):
    logger.error('Missing directory with raw files: %s', dir_raw)
    sys.exit(1)

os.makedirs(dir_proj, exist_ok=True)
os.makedirs(dir_blank, exist_ok=True)
# northern sky above dec +24, 10 arcmin pixels — ZEA centered on north pole
PIXSIZE_DEG = 10.0 / 60.0
DEC_MIN     = 24.0
CRVAL1      = 0.0    # RA of projection centre (north pole, any RA works)
CRVAL2      = 90.0   # DEC of projection centre

# ZEA projection radius (degrees) for a point at DEC_MIN:
#   R = (180/pi) * 2 * sin(theta/2),  theta = angular distance from pole
theta_max   = 90.0 - DEC_MIN          # 66 deg
R_max_deg   = math.degrees(2.0 * math.sin(math.radians(theta_max / 2.0)))
R_max_pix   = R_max_deg / PIXSIZE_DEG
NAXIS1      = 2 * math.ceil(R_max_pix) + 2   # square image with small margin
NAXIS2      = NAXIS1
CRPIX1      = NAXIS1 / 2.0 + 0.5
CRPIX2      = NAXIS2 / 2.0 + 0.5

with open(hdr_final, 'w') as _f:
    _f.write("SIMPLE  =                    T\n")
    _f.write("BITPIX  =                  -64\n")
    _f.write("NAXIS   =                    2\n")
    _f.write(f"NAXIS1  =                 {NAXIS1:d}\n")
    _f.write(f"NAXIS2  =                 {NAXIS2:d}\n")
    _f.write("CTYPE1  = 'RA---ZEA'\n")
    _f.write("CTYPE2  = 'DEC--ZEA'\n")
    _f.write(f"CRVAL1  =           {CRVAL1:12.6f}\n")
    _f.write(f"CRVAL2  =           {CRVAL2:12.6f}\n")
    _f.write(f"CRPIX1  =           {CRPIX1:12.6f}\n")
    _f.write(f"CRPIX2  =           {CRPIX2:12.6f}\n")
    _f.write(f"CDELT1  =           {-PIXSIZE_DEG:12.8f}\n")
    _f.write(f"CDELT2  =           { PIXSIZE_DEG:12.8f}\n")
    _f.write("CROTA2  =             0.000000\n")
    _f.write("LONPOLE =           180.000000\n")
    _f.write("END\n")
logger.info('Written %s: %dx%d px, %.0f arcmin/px, ZEA north pole, DEC >= %.0f deg',
            hdr_final, NAXIS1, NAXIS2, PIXSIZE_DEG * 60, DEC_MIN)

# project all files onto the same pixel grid
for fits_raw in sorted(glob.glob(dir_raw + '/HP*-mosaicI.rms.fits')):
    fits_proj  = dir_proj + '/' + os.path.basename(fits_raw).replace('.fits', '-proj.fits')
    fits_blank = dir_blank + '/' + os.path.basename(fits_raw).replace('.fits', '-blank.fits')
    if not os.path.exists(fits_proj):

        # mask on the unprojected file: keep only pixels closest to this HEALPix pixel
        _m = re.search(r'HP(\d+)', os.path.basename(fits_raw))
        if _m is not None:
            _hp_id = int(_m.group(1))
            with afits.open(fits_raw) as _hdul:
                _data = _hdul[0].data.squeeze().astype(np.float32)
                _hdr  = _hdul[0].header.copy()
            _wcs = WCS(_hdr).celestial
            _ny, _nx = _data.shape
            _xg, _yg = np.meshgrid(np.arange(_nx), np.arange(_ny))
            _ra, _dec = _wcs.all_pix2world(_xg.ravel(), _yg.ravel(), 0)
            _owner = hp_grid.skycoord_to_healpix(
                SkyCoord(ra=_ra * u.deg, dec=_dec * u.deg, frame='icrs')
            ).reshape(_ny, _nx)
            _mask = _owner != _hp_id
            _frac = _mask.sum() / _mask.size
            _data[_mask] = np.nan
            afits.writeto(fits_blank, data=_data, header=_hdr, overwrite=True)
            logger.info('  Masked to HEALPix pixel %04d (%.1f%% pixels masked) -> %s',
                        _hp_id, 100.0 * _frac, os.path.basename(fits_blank))
            fits_to_project = fits_blank
        else:
            fits_to_project = fits_raw

        logger.info('Projecting %s', fits_to_project)
        subprocess.run(['mProject', '-X', '-z', '0.5', fits_to_project, fits_proj, hdr_final], check=True)

# make image table for mAdd
subprocess.run(['mImgtbl', dir_proj, tbl_proj], check=True)
# add projected images
subprocess.run(['mAdd', '-p', dir_proj, tbl_proj, hdr_final, fits_output], check=True)

# histogram of all finite pixels with logarithmically spaced bins
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

with afits.open(fits_output) as _hdul:
    _pixels = _hdul[0].data.ravel()

_pixels = _pixels[np.isfinite(_pixels) & (_pixels > 0)]

_bins = np.logspace(np.log10(_pixels.min()), np.log10(10e-3), 100)
fig, ax = plt.subplots()
_n, _bedges, _ = ax.hist(_pixels, bins=_bins, edgecolor='black', linewidth=0.5)
_peak_idx = np.argmax(_n)
_peak_x = np.sqrt(_bedges[_peak_idx] * _bedges[_peak_idx + 1])
ax.axvline(_peak_x, color='red', linestyle='--', label=f'peak: {_peak_x:.3e} Jy/beam')
ax.legend()
ax.set_xscale('log')
ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, numticks=20))
ax.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10, subs='all', numticks=100))
ax.xaxis.set_major_formatter(matplotlib.ticker.LogFormatter(base=10, labelOnlyBase=False))
ax.xaxis.set_minor_formatter(matplotlib.ticker.LogFormatter(base=10, labelOnlyBase=False, minor_thresholds=(10, 0.4)))
ax.tick_params(axis='x', which='both', rotation=45)
ax.set_xlabel('rms noise (Jy/beam)')
ax.set_ylabel('Count')
ax.set_title('RMS mosaic pixel distribution')
fig.tight_layout()
_hist_out = fits_output.replace('.fits', '-hist.png')
fig.savefig(_hist_out, dpi=150)
plt.close(fig)
logger.info('Histogram saved to %s', _hist_out)

logger.info('Done: %s', fits_output)