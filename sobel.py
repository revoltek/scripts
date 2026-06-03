#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#./sobel.py name.img

import sys, os, argparse
from astropy.io import fits
from scipy import ndimage
import numpy as np

parser = argparse.ArgumentParser(description='Edge/ridge filter for FITS images.')
parser.add_argument('fitsfile', help='Input FITS file')
parser.add_argument('--ridge', action='store_true',
                    help='Use Hessian ridge filter instead of Sobel (better for filament spines)')
parser.add_argument('--sigmas', nargs='+', type=float, default=[1, 2, 4],
                    help='Smoothing scales for the ridge filter (default: 1 2 4 pixels)')
args = parser.parse_args()
img_f = args.fitsfile

# img: a 2-D array
def sobel(img):
    dx = ndimage.sobel(img, 0)  # horizontal derivative
    dy = ndimage.sobel(img, 1)  # vertical derivative
    mag = np.hypot(dx, dy)  # magnitude
    max_mag = np.max(mag)
    if max_mag > 0:
        mag *= 100. / max_mag  # normalize
    return mag


def ridge(img, sigmas=(1, 2, 4)):
    """
    Multi-scale Hessian ridge filter.
    Highlights the central spine of bright elongated structures (filaments)
    instead of their two edges.  For each scale sigma the Hessian H is
    computed via Gaussian-derivative filtering; the smaller eigenvalue lam2
    is strongly negative at a bright ridge and near-zero elsewhere.
    The response is max-pooled across scales so filaments of different
    widths are all captured.
    """
    result = np.zeros_like(img, dtype=float)
    for sigma in sigmas:
        Hxx = ndimage.gaussian_filter(img, sigma=sigma, order=[2, 0])
        Hyy = ndimage.gaussian_filter(img, sigma=sigma, order=[0, 2])
        Hxy = ndimage.gaussian_filter(img, sigma=sigma, order=[1, 1])
        # Eigenvalues of the symmetric 2x2 Hessian at each pixel
        tmp = np.sqrt((Hxx - Hyy) ** 2 + 4 * Hxy ** 2)
        lam2 = 0.5 * (Hxx + Hyy - tmp)  # smaller (most negative) eigenvalue
        # Bright ridge: lam2 < 0  →  ridgeness = -lam2
        ridgeness = np.where(lam2 < 0, -lam2, 0.0)
        result = np.maximum(result, ridgeness)  # keep max response across scales
    max_val = np.max(result)
    if max_val > 0:
        result *= 100.0 / max_val
    return result


img = fits.open(img_f)
original_shape = img[0].data.shape
pixels = np.squeeze(img[0].data)
print("Image shape:", pixels.shape)

#pixels = ndimage.gaussian_filter(pixels, sigma=3)
if args.ridge:
    print("Applying Hessian ridge filter (sigmas=%s)" % args.sigmas)
    pixels = ridge(pixels, sigmas=args.sigmas)
    suffix = '.ridge'
else:
    pixels = sobel(pixels)
    suffix = '.sobel'

# Write filtered pixel data (restore original shape)
img[0].data = pixels.reshape(original_shape)
img.writeto(img_f + suffix, overwrite=True)
print("Written:", img_f + suffix)
