#!/usr/bin/python3

# MeerKAT Cosine beam shape function
import os,sys
from astropy.io import fits
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Save/Correct MeerKAT primary beam.')
parser.add_argument('--savecorr', dest='savecorr', action='store_true', help='Save a beam corrected file (defailt: False)')
parser.add_argument('--savebeam', dest='savebeam', action='store_true', help='Save the beam file (default: False)')
parser.add_argument('--outputcorr', dest='outputcorr', default='pbcor', help='Name of the suffix for the corrected fits file (default: {input}_pbcor.fits)')
parser.add_argument('--outputbeam', dest='outputbeam', default='beam', help='Name of ths suffix for the output beam fits file (default: {input}_beam.fits)')
parser.add_argument('images', nargs='+', help='List of input images.')

args = parser.parse_args()

if args.savecorr == False and args.savebeam == False:
    print('Select either --savecorr or --savebeam or both.')
    sys.exit()

if len(args.images) == 0:
    print('Missing input images.')
    sys.exit()

def MKCosBeam(rho, nu):
    """
    Calculate cosine beam shape (Condon & Ransom, Essential Radio Astronomy eq 3.95)

    Return power gain of circularly symmetric beam
    * rho   = offset from center (degrees)
    * nu    = Frequency (Hz)
    * D     = Antenna diameter (m)
    * beam  = beam FWHM (amin)
    """
    ################################################################
    #theta_b = radians(57.5/60) * (1.5e9/nu)
    theta_b = 0.0167261 * (1.5e9/nu)
    rhor = 1.18896*(rho*np.pi/180)/theta_b
    div = (1.-4.*(rhor**2))
    div[abs(div)<1e-5] = 1.0e-5
    gain = (np.cos(np.pi*rhor)/div)**2
    return gain

def get_beam(data, header):
    """
    apply the beam to the the 2d array "data" using MKCosBeam
    """
    # freq
    assert 'FREQ' in header['CTYPE3']
    nu = header['CRVAL3']
    # find distance in deg from image center to each pixel
    
    pix2deg = abs(header['CDELT1']) # in deg
    pixPhaseCentre = [header['CRPIX1'], header['CRPIX2']]

    def beam_creator(i,j):
        # get distance from phase centre pixel in deg
        rho = np.sqrt( (pixPhaseCentre[0] - i)**2 + (pixPhaseCentre[1] - j)**2 ) * pix2deg
        return MKCosBeam(rho, nu)

    beam = np.fromfunction(beam_creator, [data.shape[2],data.shape[3]])
    return np.array([[beam]])

for fits_file in args.images:
    # read
    print('Work on: %s' % fits_file)
    with fits.open(fits_file, readonly=True) as f:
        # apply beam
        data = f[0].data
        header = f[0].header
        beam = get_beam(data, header)
        if args.savecorr:
            data = data/beam
            # write
            fits_output = fits_file.replace('.fits','_'+args.outputcorr+'.fits')
            print('Save: %s' % fits_output)
            fits.writeto(fits_output, data.astype(np.float32), header, overwrite=True)
        if args.savebeam:
            # write
            fits_output = fits_file.replace('.fits','_'+args.outputbeam+'.fits')
            print('Save: %s' % fits_output)
            fits.writeto(fits_output, beam.astype(np.float32), header, overwrite=True)

print("Done.")
