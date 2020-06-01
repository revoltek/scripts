#!/usr/bin/env python

import numpy as np
from astropy import constants as const
import astropy.units as u

def radio_freq(gamma, B):
    return gamma**2 * const.e * B / (const.m_e * 2 * np.pi)

radio_freq(1000, 1e6*u.Gauss)

