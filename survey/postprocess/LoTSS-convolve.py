#! /usr/bin/env python3
# Convolve LoTSS mosaics to 15" resolution

import os, glob
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel, convolve_fft
from reproject import reproject_interp
import numpy as np

dir_mosaics_lotss = "/homes/fdg/storage/catalogues/LoTSS-DR3/healpix/"
dir_mosaics_lotss_out = "/homes/fdg/storage/catalogues/LoTSS-DR3/healpix_convolved/"
os.makedirs(dir_mosaics_lotss_out, exist_ok=True)
dir_mosaics_lolss = "/homes/fdg/storage/surveytgts/mosaics/healpix/"

infiles = sorted(glob.glob(dir_mosaics_lotss + '*.fits'))

for infile in infiles:
    if os.path.exists(os.path.join(dir_mosaics_lotss_out, os.path.basename(infile))):
        print('Output file already exists, skipping %s' % infile)
        continue

    print('Working on %s' % infile)
    with fits.open(infile) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    # calculate kernel size for convolution to 15" resolution
    pixscale = np.abs(header['CDELT2'] * 3600)  # arcsec/pixel
    fwhm_current = header['BMAJ'] * 3600  # arcsec
    fwhm_target = 15.0  # arcsec
    if fwhm_current >= fwhm_target:
        print('Already at or above target resolution, skipping convolution')
        continue
    sigma_current = fwhm_current / (2 * np.sqrt(2 * np.log(2)))
    sigma_target = fwhm_target / (2 * np.sqrt(2 * np.log(2)))
    sigma_kernel = np.sqrt(sigma_target**2 - sigma_current**2) / pixscale  # in pixels

    kernel = Gaussian2DKernel(sigma_kernel)
    data_convolved = convolve_fft(data, kernel, boundary='fill', fill_value=0, allow_huge=True)

    # regrid to the WCS of the matching LoLSS file
    lolss_file = os.path.join(dir_mosaics_lolss, os.path.basename(infile))
    with fits.open(lolss_file) as hdul_lolss:
        target_header = hdul_lolss[0].header
        target_wcs = WCS(target_header)
        target_shape = hdul_lolss[0].data.shape

    data_regridded, _ = reproject_interp((data_convolved, WCS(header)), target_wcs, shape_out=target_shape)

    out_header = target_header.copy()
    for key in ('BMAJ', 'BMIN', 'BPA'):
        if key in header:
            out_header[key] = fwhm_target / 3600.0 if key in ('BMAJ', 'BMIN') else 0.0

    outname = os.path.join(dir_mosaics_lotss_out, os.path.basename(infile))
    fits.writeto(outname, data_regridded.astype(np.float32), out_header, overwrite=True)