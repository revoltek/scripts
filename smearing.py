#!/usr/bin/python3

from math import *
import numpy as np
import scipy.special

# USER INPUT
freq = 54.0e6 # Hz. Observing freqency, from msoverview
nchan= 16 # Channels per subband
bw=195.e3/nchan # Hz. 195kHz is the bandwidt of one subband
delta_T = 2. # seconds. Averaging
b = 1292.0e3 # meters. IS longest baseline in meters. 
#b = 100.0e3 # meters. RS longest baseline in meters.
delta_Theta = 2.0 # deg. Distance in image from phase centre, i.e. half the FWHM.

c = 2.998e8 # m/s. Speed of light
resolution = (c/freq)/b*180./np.pi # deg. Resolution
print("Resolution is ~",resolution*3600, "arcsec")

def time_smearing(resolution,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # output and input units are equal but print is in deg
    # Condition that delta_Theta * delta_T << 1.37E4 * Resolution
    # where delta_Theta is the angular radius of the image and delta_T is the time smoothing.

    # Given resolution and field of view (delta_theta) this returns a condition on the time smearing at the correlator
    delta_T = 1.37E4*resolution/delta_Theta
    print('Time averaging should be less than %s s.' % delta_T)

    time_smearing2(delta_T,delta_Theta,resolution)
    
    return delta_T

def time_smearing2(delta_T,delta_Theta,resolution):
	#Same as above but provides the flux loss for a given time averaging
	Reduction = 1-1.22E-9*(delta_Theta/resolution)**2.0 * delta_T**2.0
	print('At radius %s deg and resolution %s arcsec the source will have %s percent of its flux if data smoothed to %s s.' % (delta_Theta,resolution*3600,Reduction,delta_T))

	return
    
def bandwidth_smearing(resolution,freq,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # Output and input units are equal 
    # Condition that delta_Theta * delta_freq << Freq * resolution
    # where delta_Theta is the offset from the pointing centre and delta_freq is the bandwidth smoothing.

    # Given resolution, freq and offset this gives the condition for the delta_freq

    delta_freq = freq*resolution/delta_Theta

    print('Bandwidth averaging should be much less than %s chan.' % (195.e3/delta_freq))

    bandwidth_smearing2(resolution,freq,delta_Theta,delta_freq)

    return delta_freq

def bandwidth_smearing2(resolution,freq,delta_Theta,delta_freq):
	# Same as above but gives the flux loss for a given frequency averaging.

    # Bandwidth smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B
    beta = (delta_freq/freq) * (delta_Theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    Reduction = ((np.pi**0.5)/(gamma * beta)) * (scipy.special.erf(beta*gamma/2.0))

    print('At radius %s deg and resolution %s arcsec and at frequency of %s Hz a source will have %s percent of its flux if data smoothed in freq to %s chan.' % (delta_Theta,resolution*3600,freq,Reduction,195.e3/delta_freq))

    return

time_smearing(resolution, delta_Theta)
time_smearing2(delta_T,delta_Theta, resolution)
bandwidth_smearing(resolution,freq,delta_Theta)
bandwidth_smearing2(resolution,freq,delta_Theta,bw)
