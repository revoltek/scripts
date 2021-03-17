#!/usr/bin/env python3
#
# Script to trace the flux-density evolution of a path in one or more fits-files.
# Required input: a ds9 region file which contains an ordered sequence of points, fits image(s).

import os, sys, argparse
import logging as log
import numpy as np
from astropy import units
from astropy.coordinates import SkyCoord, Angle
from scipy import interpolate#, integrate, optimize
import pandas as pd
import pyregion

import lib_fits
log.root.setLevel(log.INFO)


def fit_path_to_regions(region, image, z, n):
    """
    Fit a path to a number of ds9 point regions.
    The region file must be ordered!
    Parameters
    ----------
    region: string, filename of ds9 region file
    image: object, lib_fits.Image()
    z: float, redshift
    n: int, number of points to sample
    Returns
    -------
    xy : (n,2) array of type int
        Interpolated path in pixel coordinates
    l: (n,) array of floats
        Length along the path in kpc
    """
    approxres = 1000
    # Load region, region must be sorted!
    trace = pyregion.open(region)
    trace = np.array([p.coord_list for p in trace.as_imagecoord(image.img_hdr)])

    # Linear interpolation
    distance_lin = np.cumsum(np.linalg.norm(np.diff(trace, axis=0), axis=1))
    distance_lin = np.insert(distance_lin,0,0)
    k = min([3, len(trace)-1]) # order spline
    tck, u = interpolate.splprep([trace[:,0], trace[:,1]], u=distance_lin, s=0, k=k)

    # Cubic spline interpolation of linear interpolated data to get correct distances
    # Calculate a lot of points on the spline and then use the accumulated distance from points n to point n+1 as the integrated path
    xy = interpolate.splev(np.linspace(0,u[-1],approxres), tck, ext=2)
    distance_cube = np.cumsum(np.linalg.norm(np.diff(xy, axis=1), axis=0))
    distance_cube = np.insert(distance_cube,0,0)
    tck, u = interpolate.splprep([xy[0], xy[1]], s=0)
    length = np.linspace(0, 1, n) # length at point i in fraction
    xy = np.array(interpolate.splev(length, tck, ext=2)).T # points where we sample in image coords.
    length = length*distance_cube[-1]/image.get_pixelperkpc(z) # length at point i in kpc
    log.info(f"Trace consists of {len(trace)} points. Linear interpolation length: {distance_lin[-1]/image.get_pixelperkpc(z):.2f} kpc,  cubic interpolation length: {distance_cube[-1]/image.get_pixelperkpc(z):.2f} kpc")
    if not np.all(np.isclose([trace[0], trace[-1]], [xy[0], xy[-1]])):  # Assert
        raise ValueError(f'Check points: First diff - {trace[0] - xy[0]}; last diff = {trace[-1] - xy[-1]}')
    return xy, length


def beam_ellipse(ra, dec, image):
    """
    Return pyregion ellpse regon coresponding to the image beam at ra, dec
    Parameters
    ----------
    ra: float, ra in degrees
    dec: float, dec in degrees
    image: obj, lib_fits.Image

    Returns
    -------
    beam_ellipse: obj, pyregion region
    """
    b = image.get_beam()
    ra = Angle(str(ra)+'d', unit=units.deg).to_string(sep = ':', unit=units.hour)
    dec = Angle(dec, unit=units.deg).to_string(sep = ':')
    fact = 1/np.sqrt(4*np.log(2)) # this factor is needed to make the pixels in the beam area match the beam_area from the lib_fits function. Not sure why...
    ell_str = f"ellipse({ra}, {dec}, {b[0]*3600*fact}\", {b[1]*3600*fact}\", {b[2]})"
    return pyregion.parse(ell_str)



def interpolate_path(region, image, n, z, flux_scale_err = 0.0, mode='flux'):
    """
    Interpolate a path defined by ordered ds9 points and calculate 'n' evenly spaced points on this path.
    Slide a psf-sized region along this path and calculate the mean image value along the path.
    Parameters
    ----------
    region: string
        ds9 region filename - MUST be ordered from start to end.
    image: obj, lib_fits.Image
    n: int, number of points to space on the path
    z: float, redshift

    Returns
    -------
    trace_data: array, shape(n,)
        Values of sliding psf region means.
    xy: array, shape (n,2)
        Evenly spaced points on the path
    path_length: array, shape (n,)
        Length of the path from first point to point n in kpc
    flux_scale_err: float, optional. Default = 0.0
        Relative error of the flux scale.
    """

    df = pd.DataFrame()
    xy, l = fit_path_to_regions(region, image, z, n)
    df['l'] = l
    radec = image.get_wcs().all_pix2world(xy,0) #TODO check origin
    df['ra'], df['dec'] = radec.T

    path_data = np.zeros(len(xy))
    path_data_error = np.zeros(len(xy))
    for i,p in enumerate(radec):
        beam = beam_ellipse(p[0], p[1], image).as_imagecoord(image.img_hdr)
        beam_mask = beam.get_mask(hdu=image.img_hdu, header=image.img_hdr, shape=image.img_data.shape)
        npix = np.sum(beam_mask)
        data = image.img_data[beam_mask]
        nndata = data[~np.isnan(data)]
        if mode == 'flux':
            path_data[i] = np.sum(nndata) / image.get_beam_area(unit='pixel')
            path_data_error[i] = image.noise * np.sqrt(npix / image.get_beam_area(unit='pixel'))
            print(f'Flux: {path_data[i]:.6f} +/- {path_data_error[i]:.6f}')
        elif mode == 'mean':
            path_data[i] = np.mean(nndata)
            path_data_error[i] = np.mean(nndata)/np.sqrt(len(nndata))
            print(f'Mean: {path_data[i]:.6f} +/- {path_data_error[i]:.6f}')

    if mode == 'flux':
        df[f'F_{image.freq/1e6:.0f}'] = path_data
        df[f'F_err_{image.freq/1e6:.0f}'] = path_data_error
    elif mode == 'mean':
        df[f'M_{image.freq/1e6:.0f}'] = path_data
        df[f'M_err_{image.freq/1e6:.0f}'] = path_data_error

    return df


parser = argparse.ArgumentParser(description='Trace a path defined by a points in a ds9 region file. \n    trace_path.py <fits image> <ds9 region>')
parser.add_argument('region', help='ds9 point regions defining the path, must be ordered from start to end!')
parser.add_argument('stokesi', nargs='+', default=[], help='List of fits images of Stokes I.')
parser.add_argument('--bg', default=None, type=str, help='Path to ds9 region for background estimation [default: estimate from image].')
parser.add_argument('-z', '--z', default = 0.1259, type=float, help='Source redshift. Defaults to A1033.')
parser.add_argument('--radec', dest='radec', nargs='+', type=float, help='RA/DEC where to center final image in deg (if not given, center on first image)')
parser.add_argument('-b', '--beam', default = None, type=float, help='If specified, convolve all images to a circular beam of this radius (deg). Otherwise, convolve to a circular beam with a radius equal to the largest beam major axis.')
parser.add_argument('-n', '--n', default=100, type=int, help='Number of points to sample.')
parser.add_argument('-o', '--out', default='trace_path', type=str, help='Name of the output image and csv file.')
parser.add_argument('--align', action='store_true', help='Align the images.')
parser.add_argument('--reuse-shift', action='store_true', help='Resue catalogue shifted images if available.')
parser.add_argument('--reuse-regrid', action='store_true', help='Resue regrid images if availalbe.')
parser.add_argument('--reuse-df', action='store_true', help='Resue data frame if availalbe.')
parser.add_argument('-d', '--debug', action='store_true', help='Debug output.')
parser.add_argument('-v', '--verbose', action='store_true', help='Verbosity.')
args = parser.parse_args()

if args.verbose:
    log.root.setLevel(log.INFO)

stokesi = []
all_images = lib_fits.AllImages(args.stokesi)

# find+apply shift w.r.t. first image
if args.align:
    if args.reuse_shift and np.all([os.path.exists(name.replace('.fits', '-shifted.fits')) for name in args.stokesi]):
        log.info('Reuse cat shifted images.')
        all_images = lib_fits.AllImages(
            [name.replace('.fits', '-recenter-convolve-regrid.fits') for name in args.stokesi])
    else:
        log.info('Align images to catalogue matches')
        all_images.align_catalogue()
        if args.debug: all_images.write('shifted')

# convolve images to the same beam (for now force circ beam)
if args.reuse_regrid and np.all([os.path.exists(name.replace('.fits', '-recenter-convolve-regrid.fits')) for name in args.stokesi]):
    log.info('Reuse prepared images.')
    all_images = lib_fits.AllImages([name.replace('.fits', '-recenter-convolve-regrid.fits') for name in args.stokesi])
elif len(all_images) > 1:
    log.info('Recenter, convolve and regrid all images')
    all_images.convolve_to(circbeam=True) # elliptical beam seems buggy in some cases. Also, circ beam is nice to treat covariance matrix of pixels
    if args.debug: all_images.write('recenter-convolve')
    all_images.regrid_common()
    all_images.write('recenter-convolve-regrid')

for image in all_images:
    image.calc_noise(bg_reg=args.bg) # update noise in all images TODO: which is best way? BG region??

if args.reuse_df and os.path.exists(f'{args.out}.csv'):
    df = pd.read_csv(f'{args.out}.csv')
else:
    df_list = []
    for image in all_images:
        df_list.append(interpolate_path(args.region, image, args.n, args.z))

    df = pd.concat(df_list, axis=1)#, join_axes=[df_list[0].index])
    df.reindex(df_list[0].index)
    df = df.loc[:,~df.columns.duplicated()]

    log.info(f'Save DataFrame to {args.out}.csv...')
    df.to_csv(f'{args.out}.csv')

# do the plotting
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
fig, ax = plt.subplots(constrained_layout=True)
ax.set_xlabel('distance [kpc]')

for i, im_obj in enumerate(all_images):
    im = im_obj.imagefile
    # scale data to 54MHz using beta=-0.6
    scale = (144e6/im_obj.freq)**(-0.6)

    ax.plot(df['l'], df[f'F_{im_obj.freq/1e6:.0f}']*scale, color=f'C{i}', label=r'$\nu = $' + f'{im_obj.freq/1e6:.0f} MHz')
    ax.hlines(y=0.0, xmin=df['l'].min(), xmax=df['l'].max(), linewidth=1, color='grey', ls='dotted', alpha=0.5)
    ax.fill_between(df['l'], scale*(df[f'F_{im_obj.freq/1e6:.0f}']-df[f'F_err_{im_obj.freq/1e6:.0f}']), scale*(df[f'F_{im_obj.freq/1e6:.0f}']+df[f'F_err_{im_obj.freq/1e6:.0f}']), color=f'C{i}', alpha=0.3)
    ax.set_ylabel(r'Flux density at 144MHz assuming $\alpha = -0.6$ [Jy]')

ax.legend(loc='best')
ax.set_xlim([0,df['l'].max()])
# ax.set_ylim(bottom = 0)

log.info(f'Save plot to {args.out}.pdf...')
plt.savefig(args.out+'.pdf')
