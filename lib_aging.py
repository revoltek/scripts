#!/usr/bin/env python

# Simple lib to estimate plasma aging - Henrik Edler, 2022

import sys
import multiprocessing as mp
import astropy.units as u
from scipy import special, integrate, interpolate
import logging as log
import numpy as np

### Define constants
### Throughout this script, SI units are used!!
c = 299792458.0 # m/s
m_e = 9.1093837015e-31 # kg
sigma_T = 6.6524587321e-29 # m2
mu0 = 1.25663706212e-06 # N A-2
e = 1.602176634e-19 # Coulomb
eps0 = 8.8541878128e-12 # F/m
U_cmb = 4.19e-14  # J m-3 CMB energy density at z=0
B_cmb = 3.25e-10 # T, CMB B field at z=0 in


def nu_c(E, B, alpha):
    """
    Critical frequency for synchrotron
    Parameters
    ----------
    E: float, energy in SI
    B
    alpha

    Returns
    -------

    """
    # Critical frequency
    return (3*(E / (c ** 2 * m_e)) ** 2 * e * B * np.sin(alpha) / (4 * np.pi * m_e)) # Hardcastle and Longair 3/4, Hardwood 6/4...

def F_accurate(x):
    # numerical integral defining F(x). Use this to calculate lookup-table.
    # TODO: epsrel
    return x*integrate.quad(lambda z: special.kv(5./3., z), x, np.inf, limit = 2000)[0]

def create_F_lookup():
    xvals = np.logspace(-4,1.4,1000)
    with mp.Pool() as p:
        results = p.map(F_accurate, xvals)
    np.save(__file__.replace('lib_aging.py','')+'lib_aging_data/F(X)_lookup.npy', np.array([xvals, results]))

def n_e(E, iidx, B, t, z):
    """
    Electron density taking into account Jaffe-Perola + IC losses
    Parameters
    ----------
    E: float, energy in J
    iidx: float, spectral index
    B: float, magnetic field in Tesla
    t: float, age in s
    z: float, redshift z

    Returns
    -------
    electron density, float
    """
    # C = B**2*(4*sigma_T/(6*m_e**2*c**3*mu0))
    beta = E*t*(B**2/(2*mu0) + U_cmb*(1+z)**4)*(4*sigma_T/(3*m_e**2*c**3)) # Harwood 2013 has nu_c**3 instead of c. But this can't be right?
    # slightly confusing with exponent indices here since different sign conventions exist...
    if beta >= 1:
        return 0.
    else:
        return E ** (-2 * iidx - 1) * ((1 - beta) ** ((2 * iidx + 1) - 2))


def n_e_fermi2_steady_state(E, B, q, D_0, z, C_i=0):
    """
    Electron density taking into account Synch. + IC losses and 2nd order fermit accelleration.
    Underlying assumption is a power-law energy dependency of the acceleration time -> t_acc = t_0 * p^k
    Steady state solution.
    https://www.science.org/action/downloadSupplement?doi=10.1126%2Fsciadv.1701634&file=1701634_sm.pdf
    Parameters
    ----------
    E: float, energy in J
    B: float, magnetic field in Tesla
    t_0: float, acceleration timescale normalization
    k: float, exponent of acceleration time
    z: float, redshift z

    Returns
    -------
    electron density, float
    """
    p = c**-1 * (E**2 - m_e**2 * c**4)**0.5 # probably E \approx pc would be fine
    C = (4*sigma_T * (3.2e-10)**2 )/ (3 * m_e**2 * c**3 * mu0)
    print(C)
    F = 4.8e-4 * ((B / (3.2e-10))**2 + (1+z)**4)
    C_i = 0 # integration constant

    k = 3 - q
    mp.dps = 15
    #print((k-4)/k, -(F*p**k)/(D_0*k), np.complex(mpmath.gammainc((k-4)/k, -(F*p**k)/(D_0*k))))
    sum1 = (1/F)*(-F/(D_0*k))**(4/k)*np.complex(mpmath.gammainc((k-4)/k, -(F*p**k)/(D_0*k))).real
    sum2 = C_i/D_0
    mul = - p**2 * np.exp((-F*p**k)/(D_0*k))
    #print(sum1, sum2, mul, -F*p**k/(D_0*k))
    res = mul * (sum1 + sum2)
    if res < 0:
        return 0
    else:
        return res

def n_e_fermi2_steady_state2(E, B, q, D_0, z):
    """
    Electron density taking into account Synch. + IC losses and 2nd order fermit accelleration.
    Underlying assumption is a power-law energy dependency of the acceleration time -> t_acc = t_0 * p^k
    Steady state solution ignoring the integration constant in the ODE.
    See also Stawarz et al. 2008.

    Parameters
    ----------
    E: float, energy in J
    B: float, magnetic field in Tesla
    q: float, spectral index of turbulence. 5/3 is Kolmogorov case, 2 is hard-sphere
    D_0: float, reacceleration coefficient
    z: float, redshift z

    Returns
    -------
    electron density, float
    """
    # p = c**-1 * (E**2 - m_e**2 * c**4)**0.5 # probably E \approx pc would be fine
    p = E
    C = (4*sigma_T * (3.2e-10)**2 )/ (3 * m_e**2 * c**3 * mu0)
    F = 4.8e-4 * ((B / (3.2e-10))**2 + (1+z)**4)
    F = C * ((B / (3.2e-10))**2 + (1+z)**4)
    # print(-F*p**(3-q)/(D_0*(3-q)))
    return p**2 * np.exp((-F*p**(3-q))/(D_0*(3-q)))


class S_model():

    def __init__(self, epsrel=1.e-3):
        # use lookup table for F(x)
        with open(__file__.replace('lib_aging.py', '') + 'lib_aging_data/F(X)_lookup.npy', 'rb') as f:
            xF_dat = np.load(f)
        self.F_interp = interpolate.InterpolatedUnivariateSpline(xF_dat[0], xF_dat[1], ext=0, check_finite=False)
        self.epsrel = epsrel

    def evaluate(self, nu, B, iidx, t, z, N0=1.):
        """
        Synchrotron flux density including basic Jaffe-Perola aging.
        Parameters
        ----------
        nu: float or numpy array
            Frequency in Hertz
        B: float
            B in Tesla
        iidx: float
            Injection index, positive definition
        t: float, or numpy array
            Age in Myrs
        N0: float,
            Normaliation factor, optional.

        Returns
        -------
        flux density: float or numpy array
            Arbitrary units
        """
        t *= 1e6*3.154e7 # Myrs to seconds
        nu = nu*(1+z) # redshift frequency
        C0 = (z+1)**-2*N0*3**0.5*e**3*B/(8*np.pi*eps0*c*m_e)
        E_min, E_max = 0.5e6*1.60218e-19, 1.e11*1.60218e-19 # eV, TODO: units...

        def integrand(E, alpha):
            """
            Integrand to call
            Parameters
            ----------
            E: float, energy in J
            alpha: float, impact angle
            """
            try:
                E = E[:,np.newaxis]
                alpha = alpha[np.newaxis]
                result = self.F(nu/nu_c(E,B,alpha))*0.5*np.sin(alpha)**2*n_e(E, iidx, B, t, z)
                result[np.isnan(result)] = 0.0  # case zero times infinity
                return result
            except (IndexError, TypeError) as e:
                if alpha == 0.:
                    return 0.
                else:
                    return self.F(nu/nu_c(E,B,alpha))*0.5*np.sin(alpha)**2*n_e(E, iidx, B, t, z)

        res_quad =  C0 * integrate.dblquad(integrand, 1e-4, np.pi, E_min, E_max, epsrel=self.epsrel)[0] # rough integration
        return res_quad

    def evaluate_fermi2_steady_state(self, nu, B, z, q, D_0, N0=1.):
        """
        Synchrotron flux density including Jaffe-Perola model, synch. and IC aging as well as 2nd order fermi acceleration.
        Parameters
        ----------
        nu: float or numpy array
            Frequency in Hertz
        B: float
            B in Tesla
        z: float, redshift
        t_acc: float
           acceleration timescale in Myrs
        N0: float,
            Normaliation factor, optional.

        Returns
        -------
        flux density: float or numpy array
            Arbitrary units
        """
        # t_acc *= 1e6*3.154e7 # Myrs to seconds
        nu = nu*(1+z) # redshift frequency
        C0 = (z+1)**-2*N0*3**0.5*e**3*B/(8*np.pi*eps0*c*m_e)
        E_min, E_max = 0.5e6*1.60218e-19, 1.e10*1.60218e-19 # eV, TODO: units...

        def integrand(E, alpha):
            """
            Integrand to call
            Parameters
            ----------
            E: float, energy in J
            alpha: float, impact angle
            """
            return self.F(nu/nu_c(E,B,alpha))*0.5*np.sin(alpha)**2*n_e_fermi2_steady_state2(E, B, q, D_0, z)

        return C0*integrate.dblquad(integrand, 0, np.pi, E_min, E_max, epsrel=self.epsrel)[0] # rough integration

    def F(self, x):
        return np.vectorize(self._F)(x)

    def _F(self, x):
        # F(x): Use asymptotes below and above, in between interpolate lookup table
        # Mourad Fouka1and Saad Ouichaoui, 2013
        if x > 25:
            return np.sqrt(np.pi*x/2)*np.exp(-x)
        elif x < 1e-4:
            return np.pi*2**(5/3)/(special.gamma(1/3)*np.sqrt(3))*x**(1/3)
        else:
            return self.F_interp(x)

def get_si(nu1, nu2, S1, S2): # TODO use lib_linearfit
    return np.log(S1 / S2) / np.log(nu1 / nu2)

def get_aging_si(nu1, nu2, B, injection_index, times, z, model=None):
    """
    Return the Jaffe-Perola aging path in a color-color plot.
    Parameters
    ----------
    nu1: float
        lower spectral index HERTZ.
    nu2: list of two floats
        upper spectral index HERTZ.
    B: float
        Magnetic field in Tesla
    injection_index: float
        injection spectral index (positive definition)
    times: array of floats
        times at which to evaluate the SI in Myr
    z: float
        Redshift
    Returns
    -------
    si: array, sequence of the  spectral indices at different times
    """
    try:
        times[0]
    except IndexError:
        times = [times]
    S_array = np.zeros((len(times), 2))
    if model is None:
        model = S_model()
    for i, t in enumerate(times):
        S_array[i,0] = model.evaluate(nu1, B, injection_index, t, z)
        S_array[i,1] = model.evaluate(nu2, B, injection_index, t, z)
    si = get_si(nu1, nu2, S_array[:,0], S_array[:,1])
    return si

#!/usr/bin/env python

# Simple lib to estimate plasma aging - Henrik Edler

import sys
import multiprocessing as mp
import astropy.units as u
from scipy import special, integrate, interpolate
from collections.abc import Iterable
import mpmath
import numba
import quadpy
import logging as log
import numpy as np

### Define constants
### Throughout this script, SI units are used!!
c = 299792458.0 # m/s
m_e = 9.1093837015e-31 # kg
sigma_T = 6.6524587321e-29 # m2
mu0 = 1.25663706212e-06 # N A-2
e = 1.602176634e-19 # Coulomb
eps0 = 8.8541878128e-12 # F/m
U_cmb = 4.19e-14  # J m-3 CMB energy density at z=0
B_cmb = 3.25e-10 # T, CMB B field at z=0 in

@numba.njit()
def nu_c(E, B, alpha):
    """
    Critical frequency for synchrotron
    Parameters
    ----------
    E: float, energy in SI
    B
    alpha

    Returns
    -------

    """
    # Critical frequency
    return (3*(E / (c ** 2 * m_e)) ** 2 * e * B * np.sin(alpha) / (4 * np.pi * m_e)) # Hardcastle and Longair 3/4, Hardwood 6/4...

def F_accurate(x):
    # numerical integral defining F(x). Use this to calculate lookup-table.
    return x*integrate.quad(lambda z: special.kv(5./3., z), x, np.inf, epsrel=1.e-5, limit=2000)[0]

def create_F_lookup():
    xvals = np.logspace(-4,1.4,1000)
    with mp.Pool() as p:
        results = p.map(F_accurate, xvals)
    np.save(__file__.replace('lib_aging.py','')+'lib_aging_data/F(X)_lookup.npy', np.array([xvals, results]))

# @numba.njit()
def n_e(E, iidx, B, t, z):
    """
    Electron density taking into account Jaffe-Perola + IC losses
    Parameters
    ----------
    E: float, energy in J
    iidx: float, spectral index
    B: float, magnetic field in Tesla
    t: float, age in s
    z: float, redshift z

    Returns
    -------
    electron density, float
    """
    # C = B**2*(4*sigma_T/(6*m_e**2*c**3*mu0))
    beta = E*t*(B**2/(2*mu0) + U_cmb*(1+z)**4)*(4*sigma_T/(3*m_e**2*c**3)) # Harwood 2013 has nu_c**3 instead of c. But this can't be right?
    # slightly confusing with exponent indices here since different sign conventions exist...
    if beta >= 1:
        return 0.
    else:
        return E ** (-2 * iidx - 1) * ((1 - beta) ** ((2 * iidx + 1) - 2))


def n_e_fermi2_steady_state(E, B, q, D_0, z, C_i=0):
    """
    Electron density taking into account Synch. + IC losses and 2nd order fermit accelleration.
    Underlying assumption is a power-law energy dependency of the acceleration time -> t_acc = t_0 * p^k
    Steady state solution.
    https://www.science.org/action/downloadSupplement?doi=10.1126%2Fsciadv.1701634&file=1701634_sm.pdf
    Parameters
    ----------
    E: float, energy in J
    B: float, magnetic field in Tesla
    t_0: float, acceleration timescale normalization
    k: float, exponent of acceleration time
    z: float, redshift z

    Returns
    -------
    electron density, float
    """
    p = c**-1 * (E**2 - m_e**2 * c**4)**0.5 # probably E \approx pc would be fine
    C = (4*sigma_T * (3.2e-10)**2 )/ (3 * m_e**2 * c**3 * mu0)
    print(C)
    F = 4.8e-4 * ((B / (3.2e-10))**2 + (1+z)**4)
    C_i = 0 # integration constant

    k = 3 - q
    mp.dps = 15
    #print((k-4)/k, -(F*p**k)/(D_0*k), np.complex(mpmath.gammainc((k-4)/k, -(F*p**k)/(D_0*k))))
    sum1 = (1/F)*(-F/(D_0*k))**(4/k)*np.complex(mpmath.gammainc((k-4)/k, -(F*p**k)/(D_0*k))).real
    sum2 = C_i/D_0
    mul = - p**2 * np.exp((-F*p**k)/(D_0*k))
    #print(sum1, sum2, mul, -F*p**k/(D_0*k))
    res = mul * (sum1 + sum2)
    if res < 0:
        return 0
    else:
        return res

def n_e_fermi2_steady_state2(E, B, q, D_0, z):
    """
    Electron density taking into account Synch. + IC losses and 2nd order fermit accelleration.
    Underlying assumption is a power-law energy dependency of the acceleration time -> t_acc = t_0 * p^k
    Steady state solution ignoring the integration constant in the ODE.
    See also Stawarz et al. 2008.

    Parameters
    ----------
    E: float, energy in J
    B: float, magnetic field in Tesla
    q: float, spectral index of turbulence. 5/3 is Kolmogorov case, 2 is hard-sphere
    D_0: float, reacceleration coefficient
    z: float, redshift z

    Returns
    -------
    electron density, float
    """
    # p = c**-1 * (E**2 - m_e**2 * c**4)**0.5 # probably E \approx pc would be fine
    p = E
    C = (4*sigma_T * (3.2e-10)**2 )/ (3 * m_e**2 * c**3 * mu0)
    F = 4.8e-4 * ((B / (3.2e-10))**2 + (1+z)**4)
    F = C * ((B / (3.2e-10))**2 + (1+z)**4)
    # print(-F*p**(3-q)/(D_0*(3-q)))
    return p**2 * np.exp((-F*p**(3-q))/(D_0*(3-q)))


class S_model():

    def __init__(self, epsrel=1.e-3):
        # use lookup table for F(x)
        with open(__file__.replace('lib_aging.py', '') + 'lib_aging_data/F(X)_lookup.npy', 'rb') as f:
            xF_dat = np.load(f)
        self.F_interp = interpolate.InterpolatedUnivariateSpline(xF_dat[0], xF_dat[1], ext=0, check_finite=False)
        self.epsrel = epsrel

    def evaluate(self, nu, B, iidx, t, z, N0=1.):
        """
        Synchrotron flux density including basic Jaffe-Perola aging.
        Parameters
        ----------
        nu: float or numpy array
            Frequency in Hertz
        B: float
            B in Tesla
        iidx: float
            Injection index, positive definition
        t: float, or numpy array
            Age in Myrs
        N0: float,
            Normaliation factor, optional.

        Returns
        -------
        flux density: float or numpy array
            Arbitrary units
        """
        t *= 1e6*3.154e7 # Myrs to seconds
        nu = nu*(1+z) # redshift frequency
        C0 = (z+1)**-2*N0*3**0.5*e**3*B/(8*np.pi*eps0*c*m_e)
        E_min, E_max = 0.5e6*1.60218e-19, 1.e11*1.60218e-19 # eV, TODO: units...

        def integrand(E, alpha):
            """
            Integrand to call
            Parameters
            ----------
            E: float, energy in J
            alpha: float, impact angle
            """
            try:
                E = E[:,np.newaxis]
                alpha = alpha[np.newaxis]
                result = self.F(nu/nu_c(E,B,alpha))*0.5*np.sin(alpha)**2*n_e(E, iidx, B, t, z)
                result[np.isnan(result)] = 0.0  # case zero times infinity
                return result
            except (IndexError, TypeError) as e:
                if alpha == 0.:
                    return 0.
                else:
                    return self.F(nu/nu_c(E,B,alpha))*0.5*np.sin(alpha)**2*n_e(E, iidx, B, t, z)

        res_quad =  C0 * integrate.dblquad(integrand, 1e-4, np.pi, E_min, E_max, epsrel=self.epsrel)[0] # rough integration
        return res_quad

    def evaluate_fermi2_steady_state(self, nu, B, z, q, D_0, N0=1.):
        """
        Synchrotron flux density including Jaffe-Perola model, synch. and IC aging as well as 2nd order fermi acceleration.
        Parameters
        ----------
        nu: float or numpy array
            Frequency in Hertz
        B: float
            B in Tesla
        z: float, redshift
        t_acc: float
           acceleration timescale in Myrs
        N0: float,
            Normaliation factor, optional.

        Returns
        -------
        flux density: float or numpy array
            Arbitrary units
        """
        # t_acc *= 1e6*3.154e7 # Myrs to seconds
        nu = nu*(1+z) # redshift frequency
        C0 = (z+1)**-2*N0*3**0.5*e**3*B/(8*np.pi*eps0*c*m_e)
        E_min, E_max = 0.5e6*1.60218e-19, 1.e10*1.60218e-19 # eV, TODO: units...

        def integrand(E, alpha):
            """
            Integrand to call
            Parameters
            ----------
            E: float, energy in J
            alpha: float, impact angle
            """
            return self.F(nu/nu_c(E,B,alpha))*0.5*np.sin(alpha)**2*n_e_fermi2_steady_state2(E, B, q, D_0, z)

        return C0*integrate.dblquad(integrand, 0, np.pi, E_min, E_max, epsrel=self.epsrel)[0] # rough integration

    def F(self, x):
        return np.vectorize(self._F)(x)

    def _F(self, x):
        # F(x): Use asymptotes below and above, in between interpolate lookup table
        # Mourad Fouka1and Saad Ouichaoui, 2013
        if x > 25:
            return np.sqrt(np.pi*x/2)*np.exp(-x)
        elif x < 1e-4:
            return np.pi*2**(5/3)/(special.gamma(1/3)*np.sqrt(3))*x**(1/3)
        else:
            return self.F_interp(x)

def get_si(nu1, nu2, S1, S2):
    return np.log(S1 / S2) / np.log(nu1 / nu2)

def get_aging_si(nu1, nu2, B, injection_index, times, z, model=None):
    """
    Return the Jaffe-Perola aging path in a color-color plot.
    Parameters
    ----------
    nu1: float
        lower spectral index HERTZ.
    nu2: list of two floats
        upper spectral index HERTZ.
    B: float
        Magnetic field in Tesla
    injection_index: float
        injection spectral index (positive definition)
    times: array of floats
        times at which to evaluate the SI in Myr
    z: float
        Redshift
    Returns
    -------
    si: array, sequence of the  spectral indices at different times
    """
    try:
        times[0]
    except IndexError:
        times = [times]
    S_array = np.zeros((len(times), 2))
    if model is None:
        model = S_model()
    for i, t in enumerate(times):
        S_array[i,0] = model.evaluate(nu1, B, injection_index, t, z)
        S_array[i,1] = model.evaluate(nu2, B, injection_index, t, z)
    si = get_si(nu1, nu2, S_array[:,0], S_array[:,1])
    return si

def get_aging_si_steady_state(nu1, nu2, B, q, D_0, z, model=None):
    """
    Return the Jaffe-Perola aging path in a color-color plot.
    Parameters
    ----------
    nu1: float
        lower spectral index HERTZ.
    nu2: list of two floats
        upper spectral index HERTZ.
    B: float
        Magnetic field in Tesla
    q: float
        spectral index of momentum space diffusion coefficient
    D_0: float
        reacceleration factor
    z: float
        Redshift
    Returns
    -------
    si: array,
    """
    try:
        D_0[0]
    except IndexError:
        D_0 = [D_0]
    S_array = np.zeros((len(D_0), 2))
    if model is None:
        model = S_model()
    for i, t in enumerate(D_0):
        S_array[i,0] = model.evaluate_fermi2_steady_state(nu1,B,z,q,np.array(D_0))
        S_array[i,1] = model.evaluate_fermi2_steady_state(nu2,B,z,q,np.array(D_0))
    si = get_si(nu1, nu2, S_array[:,0], S_array[:,1])
    return si


def get_model_si_vs_B(nu1, nu2, B_range, injection_index, z, t):
    """
    Return the Jaffe-Perola aging path in a color-color plot.
    Parameters
    ----------
    nu1: float
        lower spectral index HERTZ.
    nu2: list of two floats
        upper spectral index HERTZ.
    B_range: list of len 2, [lower, upper]
        Magnetic field in Tesla
    injection_index: float
        injection spectral index (positive definition)
    z: float
        Redshift
    t:  of floats
        times at which to evaluate the SI in Myr
    Returns
    -------
    si: array, sequence of the  spectral indices at different times
    """
    S_array = np.zeros((len(B_range), 2))
    for i, B in enumerate(B_range):
        S_array[i,0] = S_model(nu1, B, injection_index, t, z)
        S_array[i,1] = S_model(nu2, B, injection_index, t, z)
    si = get_si(nu1, nu2, S_array[:,0], S_array[:,1])
    return si

def characteristic_lifetime(nu, B, z):
    """
    Characteristic lifetime of electrons observed at frequency nu
    Taken from van Weeren 2018 review paper

    Parameters
    ----------
    nu: float, freq in Hz
    B: float, B in Tesla
    z: float, redshift

    Returns
    -------
    t_age: float, characteristic lifetime in Myr
    """
    log.error('not fully implemented check RAiSE III: 3C radio AGN energetics and composition to fix ')
    sys.exit()

    return 3.2e10*(B**0.5 / (B**2+(B_cmb*(1+z)**2)**2))*((1+z)*nu)**-0.5*1e-6

def plot_S_model():
    # Debug plotting
    B = 5e-10
    print('huhu')
    nu_range = np.logspace(np.log10(30e6), 9, 6)
    age_range = np.linspace(0, 200, 5)

    from agnpy.emission_regions import Blob
    from agnpy.synchrotron import Synchrotron
    blob = Blob(z=0.001, B = 10**4*B*u.gauss, spectrum_dict = {"type": "PowerLaw", "parameters": {"p": 2.3,"gamma_min": 2,"gamma_max": 1e7}})
    synch = Synchrotron(blob)
    sed = synch.sed_flux(nu_range*u.Hz)
    sed = sed.value / nu_range

    results = np.zeros((len(age_range), len(nu_range)))
    for i, age in enumerate(age_range):
        with mp.Pool() as p:
            results[i] = p.starmap(S_model, [[nu, B, 0.65, 1000, age] for nu in nu_range])
        print(results[i], nu_range)

    PL = (nu_range**-0.65)
    PL /= (PL[0]/sed[0])
    print((results[0,0]/sed[0]))
    results /= (results[0,0]/sed[0])

    import matplotlib.pyplot as plt
    plt.close()
    print(nu_range, sed)
    plt.plot(nu_range, sed, c='k', label=f'AGNPY for 0 Myr; B = {B}T')
    plt.plot(nu_range, PL, label=f'PL alpha = 0.65', c='k', ls='dotted')
    for age, res in zip(age_range, results):
        plt.plot(nu_range, res, label=f'{age}Myr; B = {B}T')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('frequency [Hz]')
    plt.ylabel('S')
    # plt.xlim([np.min(nu_range), np.max(nu_range)])
    # plt.ylim([np.min(res), 1.05*np.max(res)])
    plt.legend()
    plt.savefig(__file__.replace('lib_aging.py','')+'lib_aging_data/synch_vs_nu.png')

