#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021 - Francesco de Gasperin
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

import os, sys, argparse, logging
import numpy as np
from astropy.io import fits as pyfits
from lib_linearfit import linear_fit, linear_fit_bootstrap
from lib_fits import AllImages
logging.root.setLevel(logging.DEBUG)

def get_args():
    parser = argparse.ArgumentParser(description='Make spectral index maps, e.g. spidxmap.py --region ds9.reg --noise --sigma 5 --save *fits')
    parser.add_argument('images', nargs='+', help='List of images to use for spidx')
    parser.add_argument('--ncpu', dest='ncpu', default=1, type=int, help='Number of cpus to use (default: 1)')
    parser.add_argument('--beam', dest='beam', nargs=3, type=float, help='3 parameters final beam to convolve all images (BMAJ (arcsec), BMIN (arcsec), BPA (deg))')
    parser.add_argument('--region', dest='region', type=str, help='Ds9 region to restrict analysis')
    parser.add_argument('--noiseregion', dest='noiseregion', type=str, help='Ds9 region to calculate rms noise (default: do not use)')
    parser.add_argument('--noisesigma', dest='noisesigma', default=5, type=int, help='Sigma used in the calc_noise function when no region is specified (default: 5)')
    parser.add_argument('--size', dest='size', nargs=2, type=float, help='Size (ra and dec) of final image in degree (example: 3.5 4.0)')
    parser.add_argument('--radec', dest='radec', nargs=2, type=float, help='RA/DEC where to center final image in deg (if not given, center on first image - example: 32.3 30.1)')
    parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before calculating spidx (default: false)')
    parser.add_argument('--noise', dest='noise', action='store_true', help='Calculate noise of each image, necessary for the error map (default: false)')
    parser.add_argument('--fluxerr', dest='fluxerr', type=float, help='Fractional flux density error to be added in quadrature (default: 0 - example: 0.05 for 5 per cent)')
    parser.add_argument('--save', dest='save', action='store_true', help='Save intermediate results (default: false)')
    parser.add_argument('--sigma', dest='sigma', type=float, help='Restrict to pixels above this sigma in all images')
    parser.add_argument('--circbeam', dest='circbeam', action='store_true', help='Force final beam to be circular (default: False, use minimum common beam area)')
    parser.add_argument('--bootstrap', dest='bootstrap', action='store_true', help='Use bootstrap to estimate errors (default: use normal X|Y with errors)')
    parser.add_argument('--curvature', dest='curvature', action='store_true', help='Estimate curvature (default: only spectral index)')
    parser.add_argument('--output', dest='output', default='spidx.fits', type=str, help='Name of output mosaic (default: spidx.fits)')

    args = parser.parse_args()

    # check input
    if len(args.images) < 2:
        logging.error('Requires at lest 2 images.')
        sys.exit(1)
    if len(args.images) < 3 and args.curvature:
        logging.error('Curvature requires at least 3 images.')
        sys.exit(1)
    if args.curvature and args.bootstrap:
        logging.error('Curvature not yet implemented with bootstrap.')
        sys.exit(1)
    if args.curvature and args.ncpu > 1:
        logging.error('Curvature not yet implemented with multi-cpu.')
        sys.exit(1)
    
    return args


if __name__ == '__main__':
    args = get_args()

    ########################################################
    # prepare images and convolve+regrid if needed
    if np.all([os.path.exists(name.replace('.fits', '-conv-regrid.fits')) for name in args.images]):
        logging.info('Found convolved+regridded image... restoring.')
        all_images = AllImages([name.replace('.fits', '-conv-regrid.fits') for name in args.images])
        regrid_hdr = all_images.images[0].img_hdr
    else:
        all_images = AllImages(args.images)
        all_images.convolve_to(beam=args.beam, circbeam=args.circbeam) # no beam convolve to smallest beam
        if args.save: 
            all_images.write('conv', inflate=True) 
        regrid_hdr = all_images.regrid_common(size=args.size, region=args.region, radec=args.radec, action='regrid_header')
        if args.save: 
            all_images.write('conv-regrid', inflate=True)
        
    #####################################################
    # find+apply shift w.r.t. lowest noise image
    if args.shift:
        all_images.align_catalogue()
    
    #########################################################
    # only after regrid+convolve apply mask and find noise
    for image in all_images:
        if args.noise:
            if args.sigma is not None:
                # usually the sigma used for the blanking is rather low, better to increase it for the calc_noise
                image.calc_noise(sigma=args.noisesigma, bg_reg=args.noiseregion, force_recalc=True) # after convolution
                image.blank_noisy(args.sigma)
            else:
                image.calc_noise(force_recalc=True) # after mask?/convolution

        if args.region is not None:
            image.apply_region(args.region, invert=True) # after convolution to minimise bad pixels

    #########################################################
    # do spdix and write output
    xsize = regrid_hdr['NAXIS1']
    ysize = regrid_hdr['NAXIS2']
    frequencies = [ image.get_freq() for image in all_images ]
    if args.noise: 
        rmserr = np.array([ image.noise for image in all_images ])
    else: yerr = None
    spidx_data = np.empty(shape=(ysize,xsize))
    spidx_data[:] = np.nan
    spidx_err_data = np.empty(shape=(ysize,xsize))
    spidx_err_data[:] = np.nan
    spidx_frac_data = np.empty(shape=(ysize,xsize))
    spidx_frac_data[:] = np.nan
    if args.curvature:
        curvature_data = np.empty(shape=(ysize,xsize))
        curvature_data[:] = np.nan
        curvature_err_data = np.empty(shape=(ysize,xsize))
        curvature_err_data[:] = np.nan

    if args.ncpu > 1:
        from lib_multiproc import multiprocManager
        def funct(i,j,frequencies, val4reg, yerr, bootstrap, outQueue=None):
            if bootstrap:
                (a, b, sa, sb) = linear_fit_bootstrap(x=frequencies, y=val4reg, yerr=yerr, tolog=True)
            else:
                (a, b, sa, sb) = linear_fit(x=frequencies, y=val4reg, yerr=yerr, tolog=True)
            
            if outQueue is not None:
                outQueue.put([i,j,a,sa])

        # start processes for multi-thread
        mpm = multiprocManager(args.ncpu, funct)
        for i in range(ysize):
            print('%i/%i' % (i+1,ysize), end='\r')
            sys.stdout.flush()
            for j in range(xsize):
                val4reg = np.array([ image.img_data[i,j] for image in all_images ])
                if np.isnan(val4reg).any() or (np.array(val4reg) < 0).any(): continue
                # add flux error
                if args.fluxerr: 
                    yerr = np.sqrt((args.fluxerr*val4reg)**2+rmserr**2)
                else:
                    yerr = rmserr
                if (np.array(val4reg) <= 0).any(): continue
                mpm.put([i,j,frequencies, val4reg, yerr, args.bootstrap])

        print("Computing...")
        mpm.wait()
        for r in mpm.get():
            i = r[0]; j = r[1]; spidx = r[2]; spidx_err = r[3]
            spidx_data[i,j] = spidx
            spidx_err_data[i,j] = spidx_err
            spidx_frac_data[i,j] = abs(spidx_err/spidx)
    else:
        for i in range(ysize):
            print('%i/%i' % (i+1,ysize), end='\r')
            sys.stdout.flush()
            for j in range(xsize):
                val4reg = np.array([ image.img_data[i,j] for image in all_images ])
                if np.isnan(val4reg).any() or (np.array(val4reg) < 0).any(): continue
                # add flux error
                if args.fluxerr and args.noise: 
                    yerr = np.sqrt((args.fluxerr*val4reg)**2+rmserr**2)
                elif args.noise:
                    yerr = rmserr
                if (np.array(val4reg) <= 0).any(): continue
                if args.bootstrap:
                    (spidx, b, spidx_err, sb) = linear_fit_bootstrap(x=frequencies, y=val4reg, yerr=yerr, tolog=True)
                elif len(frequencies) == 2:
                    spidx = (np.log10(val4reg[1])-np.log10(val4reg[0]))/(np.log10(frequencies[1])-np.log10(frequencies[0]))
                    spidx_err = 1. / (np.log10(frequencies[1])-np.log10(frequencies[0])) * ((yerr[0]/val4reg[0])**2 + (yerr[1]/val4reg[1])**2) ** 0.5
                elif args.curvature:
                    # polynomial fit for curvature
                    x = np.log10(frequencies)
                    y = np.log10(val4reg)
                    # fit y = a + b*x + c*x^2
                    A = np.vstack([np.ones(len(x)), x, x**2]).T
                    if args.noise:
                        W = np.diag(1/(yerr/val4reg)**2)  # weight matrix
                        coeffs, residuals, rank, s = np.linalg.lstsq(A.T @ W @ A, A.T @ W @ y, rcond=None)
                        # error estimation for weighted fit
                        cov = np.linalg.inv(A.T @ W @ A)
                        spidx_err = np.sqrt(cov[1,1])  # error on linear coefficient (spectral index)
                        curvature_err = np.sqrt(cov[2,2])  # error on curvature coefficient
                    else:
                        coeffs, residuals, rank, s = np.linalg.lstsq(A, y, rcond=None)
                        if len(residuals) > 0:
                            mse = residuals[0] / (len(y) - 3)
                            cov = mse * np.linalg.inv(A.T @ A)
                            spidx_err = np.sqrt(cov[1,1])
                            curvature_err = np.sqrt(cov[2,2])
                        else:
                            spidx_err = np.nan
                            curvature_err = np.nan
                    spidx = coeffs[1]  # linear coefficient (spectral index)
                    curvature = coeffs[2]  # curvature coefficient
                else:
                    (spidx, b, spidx_err, sb) = linear_fit(x=frequencies, y=val4reg, yerr=yerr, tolog=True)
                spidx_data[i,j] = spidx
                spidx_err_data[i,j] = spidx_err
                spidx_frac_data[i,j] = abs(spidx_err/spidx)
                if args.curvature:
                    curvature_data[i,j] = curvature
                    curvature_err_data[i,j] = curvature_err

    if 'FREQ' in regrid_hdr.keys():
        del regrid_hdr['FREQ']
    if 'RESTFREQ' in regrid_hdr.keys():
        del regrid_hdr['RESTFREQ']
    regrid_hdr['BTYPE'] = 'SPIDX'

    if args.output[-5:] == '.fits':
        filename_out = args.output
    else:
        filename_out = args.output+'.fits'
    logging.info('Save %s (and errors)' % filename_out)
    regrid_hdr.add_history(' '.join(sys.argv))
    regrid_hdr.add_history(' - Frequencies: '+ ','.join(['%s' % (f/1e6) for f in frequencies]) + ' MHz.')
    # all images have the same beam after convolution, so we can just take the beam info from the first one
    bamj, bmin, bpa = all_images[0].get_beam()
    regrid_hdr['BMAJ'] = bamj; regrid_hdr['BMIN'] = bmin; regrid_hdr['BPA'] = bpa
    regrid_hdr['CTYPE3'] = 'ALPHA'
    pyfits.writeto(filename_out, spidx_data, regrid_hdr, overwrite=True, output_verify='fix')
    regrid_hdr['CTYPE3'] = 'ALPHAERR'
    pyfits.writeto(filename_out.replace('.fits','-err.fits'), spidx_err_data, regrid_hdr, overwrite=True, output_verify='fix')
    regrid_hdr['CTYPE3'] = 'ALPHAFRAC'
    pyfits.writeto(filename_out.replace('.fits','-err-frac.fits'), spidx_frac_data, regrid_hdr, overwrite=True, output_verify='fix')
    if args.curvature:
        regrid_hdr['CTYPE3'] = 'CURVATURE'
        pyfits.writeto(filename_out.replace('.fits','-curvature.fits'), curvature_data, regrid_hdr, overwrite=True, output_verify='fix')
        regrid_hdr['CTYPE3'] = 'CURVATUREERR'
        pyfits.writeto(filename_out.replace('.fits','-curvature-err.fits'), curvature_err_data, regrid_hdr, overwrite=True, output_verify='fix')