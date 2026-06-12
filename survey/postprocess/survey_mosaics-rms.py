#!/usr/bin/env python3

import os, sys, glob, logging, subprocess, math, re, warnings, argparse
from multiprocessing import Pool
import numpy as np
from astropy.io import fits as afits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy_healpix import HEALPix

warnings.filterwarnings('ignore', category=FITSFixedWarning)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Make RMS mosaic from catalogues')
parser.add_argument('-j', '--njobs', dest='njobs', type=int, default=12,
                    help='Number of parallel mProject processes (default: 12)')
args = parser.parse_args()

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
def project_one(fits_raw):
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
    else:
        logger.info('Already projected: %s', fits_proj)

with Pool(processes=args.njobs) as pool:
    pool.map(project_one, sorted(glob.glob(dir_raw + '/HP*-mosaicI.rms.fits')))

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
_pixels *= 1e3   # convert to mJy/beam

# total sky area covered: one HEALPix pixel per input file
_hp_ids_used = set()
for _f in glob.glob(dir_raw + '/HP*-mosaicI.rms.fits'):
    _m = re.search(r'HP(\d+)', os.path.basename(_f))
    if _m is not None:
        _hp_ids_used.add(int(_m.group(1)))
totarea = len(_hp_ids_used) * hp_grid.pixel_area.to(u.deg**2).value
logger.info('Total area: %.1f deg^2 from %d HEALPix pixels (nside=%d)',
            totarea, len(_hp_ids_used), healpix_nside)

# plot rms
fig = plt.figure(figsize=(6, 4))
fig.subplots_adjust(wspace=0, hspace=0)
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.tick_params(direction='in', top=True, right=True)

ax.set_xlabel(r'rms noise (mJy beam$^{-1}$)')
ax.set_yticklabels([])

_bins = np.logspace(np.log10(_pixels.min()), np.log10(10), 500)
n, bins, patches = ax.hist(_pixels, bins=_bins, orientation='vertical', histtype='bar', alpha=.9)

cm = plt.colormaps['inferno']
col = (n-n.min())/(n.max()-n.min())
for c, p in zip(col, patches):
    p.set_facecolor(cm(c))

ax2.hist(_pixels, bins=_bins, density=True, histtype='step', cumulative=True, label='Empirical', color='gray')
ax2.set_ylabel(r'Cumulative area (deg$^2$)')
ax2.set_yticks([0.5,1.0])
ax2.set_yticklabels([str(int(x))+' '+y for x,y in zip(np.array([0.5,1.0])*totarea, ['\n[50%]', '\n[100%]'])])

#ax.set_ylim(ymin=0,ymax=4)
median = np.median(_pixels)
print('Median rms: %f' % median)
thirty = np.quantile(_pixels,0.3)
print('30%% quantile: %f' % thirty)

fifty = np.quantile(_pixels,0.5)
print('50%% quantile: %f' % fifty)
commonrms = bins[np.argmax(n)]
print('Most common value: %f' % commonrms)
ax.set_xlim(xmin=0.8,xmax=3.8)
ax.axvline(commonrms, ls="--", color='r')
ax.axvline(fifty, ls="--", color='c')
ax2.axhline(0.5, ls="--", color='c')
ax.tick_params(axis=u'both', which=u'both',length=0)
_hist_out = fits_output.replace('.fits', '-hist.png')
fig.savefig(_hist_out, bbox_inches='tight', facecolor='w')

logger.info('Done: %s', fits_output)