#!/usr/bin/python3

# list of hips:
# https://aladin.u-strasbg.fr/hips/list

import urllib
import urllib.error

import argparse, os
import sys
import logging
import PIL
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astroquery.hips2fits import hips2fits
from astroquery.skyview import SkyView
from astropy.visualization import (SqrtStretch, PercentileInterval,
                                   LinearStretch, LogStretch,
                                   ImageNormalize, AsymmetricPercentileInterval)
from lib_plot import addRegion, addCbar, addBeam, addScalebar
from lib_fits import flatten


def get_overlay_image(coord: SkyCoord, size, wcs: WCS) -> [np.array]:
    """
    Returns an image from the SDSS image query in bytes format. Can be used
    to set as RGB background for FITSimages in aplpy.

    :param coord: Coordinate object for the center
    :param size: float, size of the cutout in deg
    :param wcs: WCS of the fits file
    :return: img, io.BytesIO object, visual image of the chosen region
    """
    image: np.array = None
    try:
        print("Trying SDSS")
        # Check if image is in SDSS
        SkyView.get_images(position=f"{coord.ra.degree}, {coord.dec.degree}", radius=size.to('deg'), survey=["SDSSi"])
        # Load the Image
        image = PIL.Image.fromarray(hips2fits.query_with_wcs("CDS/P/SDSS9/color-alt", wcs=wcs, format="png"))
    except urllib.error.HTTPError as e:
        print("Not found in SDSS")

    if image is None:
        try:
            print("Trying DSS")
            # Check if image is in DSS
            SkyView.get_images(position=f"{coord.ra.degree}, {coord.dec.degree}", radius=size.to('deg'), survey=["DSS2 Blue"])
            # Load the Image
            image = PIL.Image.fromarray(hips2fits.query_with_wcs("CDS/P/DSS2/color", wcs=wcs, format="png"))
        except urllib.error.HTTPError as e:
            print("Not found in DSS, abort", e)
    return image


parser = argparse.ArgumentParser(description='Plotting script to overlay fits contours on SDSS image')
parser.add_argument('image', help='fits image to plot.')
parser.add_argument('target', help='Name of target (Messier, NGC, VCC, IC...) or "ra,dec" in deg')
parser.add_argument('-s', '--size', type=float, default=8., help='size in arcmin')
parser.add_argument('-z', '--redshift', type=float, default=None, help='redshift.')
parser.add_argument('-n', '--noise', type=float, default=0.01, help='Hardcode noise level in Jy/beam.')
parser.add_argument('--noisemap', type=str, help='Noise map in Jy/beam.')
parser.add_argument('--overwrite', action='store_true', help='Overwrite?')
parser.add_argument('-o', '--outfile', default=None, help='prefix of output image')

args = parser.parse_args()
if args.image == None:
    logging.error('No input image found.')
    sys.exit()

# Usage:
fontsize = 12
name = args.target  # NGC ... and M.. names definitely work
outfile = args.outfile if args.outfile else args.target
if os.path.exists(outfile + '.png') and not args.overwrite:
    print(f'{outfile}.png exists - exiting.')
    sys.exit(0)
size = args.size * u.arcmin  # Size of the image in arcmin (so 10'x10')
# Load FITS file
header, data = flatten(args.image)
wcs = WCS(header)
if ',' in name:
    coord = [float(x) for x in name.split(',')]
    coord = SkyCoord(coord[0]*u.deg, coord[1]*u.deg)
else:
    coord = SkyCoord.from_name(name)
# Cutout central region
cutout = Cutout2D(data, coord, size=size, wcs=wcs)
noise = args.noise
if args.noisemap:
    print(f"Using {args.noisemap} to calculate noise...")
    # Load FITS file
    _, data_n = flatten(args.noisemap)
    cutout_n = Cutout2D(data_n, coord, size=size, wcs=wcs)
    noise = np.nanstd(cutout_n.data)
    print(f"Found background rms: {noise:.3f}mJy/beam.")

wcs=cutout.wcs
image = get_overlay_image(coord, size=size, wcs=wcs)

# Plot image with e.g. matplotlib
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(1, 1, 1, projection=wcs, slices=('x', 'y'))
lon = ax.coords['ra']
lat = ax.coords['dec']
ax.imshow(image, origin="lower", cmap="gray", interpolation='kaiser')

contour_limits = 3 * 2 ** np.arange(20) * noise
print('Make contours...', contour_limits)
vmin, vmax = contour_limits[0], np.max(cutout.data)
norm = ImageNormalize(data, vmin=vmin, vmax=vmax, stretch=SqrtStretch())
ax.contour(cutout.data, levels=contour_limits, cmap='Reds', alpha=0.8, norm=norm, linewidths=0.5)
#ax.contour(cutout.data, levels=-contour_limits[::-1], cmap='Reds', alpha=0.8, linewidths=1, linestyles='dashed', norm=norm)

if args.redshift:
    addScalebar(ax, wcs, args.redshift, 10, fontsize, color='white')
addBeam(ax, header, edgecolor='white')

lon.set_axislabel('Right Ascension (J2000)', fontsize=fontsize)
lat.set_axislabel('Declination (J2000)', fontsize=fontsize)
lon.set_ticklabel(size=fontsize)
lat.set_ticklabel(size=fontsize)

# small img
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm')
lat.set_ticklabel(rotation=90)  # to turn dec vertical

plt.title(f"{args.target}, "+r"$\sigma_{rms}=$"+f"{noise:.3f}mJy/beam", size=fontsize+2)
plt.savefig(outfile + '.png', bbox_inches='tight')
