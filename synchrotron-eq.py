#!/usr/bin/env python3

import numpy as np
from astropy import constants as const
import astropy.units as u

def eV_from_gamma( gamma ):
    return const.m_e * const.c**2 * gamma
def gamma_from_eV( eV ):
    return eV.to(u.J) / (const.m_e * const.c**2)

print( eV_from_gamma(1000).to(u.eV) )
print( gamma_from_eV(1e6*u.eV) )

def freq_from_gamma(gamma, B):
    """
    get the radio freq from a gamma and B
    """
    return gamma**2 * 1.6021766208e-19*u.C * B / (const.m_e * 2 * np.pi)

# e.g. get the MHz of electrons at gamma=1000
print (freq_from_gamma(1000, 20e-6*u.Gauss).to(u.MHz))


def gamma_from_freq(freq, B):
    """
    get the gamma from observing freq and a B
    """
    return np.sqrt((freq * const.m_e * 2 * np.pi) / (1.6021766208e-19*u.C * B)).to(u.dimensionless_unscaled)

# e.g. get the gamma of electrons emitting at 14 MHz in a 5 uG field
print ( gamma_from_freq( 14e6*u.Hz, 5e-6*u.Gauss).to(u.dimensionless_unscaled) )


def B_cmb (z):
    return (1e-6*3.25*(1+z)**2)*u.G
def B_min (z):
    return B_cmb(z)/np.sqrt(3)

print ("B_cmb:", B_cmb(0.225))

def age(gamma, B, z):
    B_ic = B_cmb(z)
    return 2.5e13*u.yr / (( (B/(1e-6*u.G))**2+(B_ic/(1e-6*u.G))**2) * gamma)

B_min = B_min(0.225)/2
print ('B_min:', B_min)
gamma = gamma_from_freq(58e6*u.Hz, B_min)
print ('gamma:', gamma)
age = age(gamma, B_min, 0.225).to(u.yr)
print('age:', age)



