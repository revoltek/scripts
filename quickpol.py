#!/usr/bin/env python3
import argparse
import numpy as np
from astropy.io import fits
import os
def compute_polarisation(Q, U):
    """Compute the polarisation as sqrt(Q^2 + U^2)"""
    P = np.sqrt(Q**2 + U**2)
    return P
def compute_polarisation_fraction(I, Q, U):
    """Compute the polarisation fraction as sqrt(Q^2 + U^2) / I"""
    P_lin = np.sqrt(Q**2 + U**2)
    with np.errstate(invalid='ignore', divide='ignore'):
        frac_pol = np.where(I != 0, P_lin / I, 0)
    return frac_pol
def compute_polarisation_angle(Q, U):
    """Compute the polarisation angle in degrees as 0.5 * arctan(U / Q)"""
    pol_angle = 0.5 * np.arctan2(U, Q)  # Arctan2 handles the quadrant correctly
    pol_angle_deg = np.degrees(pol_angle)
    return pol_angle_deg
def load_fits_data(file_path):
    """Load FITS file data."""
    with fits.open(file_path) as hdulist:
        data = hdulist[0].data
        header = hdulist[0].header
    return data, header
def save_fits(data, header, output_file):
    """Save the data as a FITS file."""
    hdu = fits.PrimaryHDU(data, header=header)
    hdu.writeto(output_file, overwrite=True)
    print(f"Saved {output_file}")
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Compute polarisation fraction and angle from Stokes I, Q, U fits files")
    parser.add_argument("--I", required=True, help="FITS file for Stokes I")
    parser.add_argument("--Q", required=True, help="FITS file for Stokes Q")
    parser.add_argument("--U", required=True, help="FITS file for Stokes U")
    parser.add_argument("--output_prefix", "-o", required=True, help="Output prefix for the generated FITS files")
    args = parser.parse_args()
    # Load data from FITS files
    I_data, I_header = load_fits_data(args.I)
    Q_data, Q_header = load_fits_data(args.Q)
    U_data, U_header = load_fits_data(args.U)
    
    # Compute polarisation fraction and angle
    pol = compute_polarisation(Q_data, U_data)
    pol_fraction = compute_polarisation_fraction(I_data, Q_data, U_data)
    pol_angle = compute_polarisation_angle(Q_data, U_data)
    
    # Save the polarisation fraction and angle as FITS files
    save_fits(pol, I_header, f"{args.output_prefix}_polarisation.fits")
    save_fits(pol_fraction, I_header, f"{args.output_prefix}_polarisation_fraction.fits")
    save_fits(pol_angle, I_header, f"{args.output_prefix}_polarisation_angle.fits")
if __name__ == "__main__":
    main()
