#!/usr/bin/python

# Python script to calculate smearing in time, frequency and also residual delays and rates.
# Written by Eskil Varenius, spring 2014. The frequency averaging is exact, while the time
# averaging is a simple upper limit based on a point source at the north pole.
# Note that the calcualtions use the longest baseline, the baselines will infact be shorter
# at all times since there is always some projection angle to the source. Still,
# the numbers are useful as an upper bound.

# The total coherence loss should not be more than 10%, else the signal to noise might be too low even if correcting for
# the amplitude loss on longest baselines. (Without correcting for the amplitude loss, one will also get larger sources
# far from the phase center, which might affect the science - so beware of averaging too much!)

# Note that this script does NOT include any primary beam effects, make sure you are not limited by these.

# A longer explanation of the formulas is in Chapter 3 of the licentiate thesis by Eskil.

import numpy as np

# USER INPUT
nu0 = 60.0e6 # Hz, Observing freqency, from msoverview
nchan= 8. # Channels per subband
bw=195.0e3/nchan # in Hz. 195kHz is the bandwidt of one subband
time = 1. # seconds. 
#b = 1158.0e3 # meters. IS longest baseline in meters. 
b = 121.0e3 # meters. RS longest baseline in meters.
imoffset = 2.*60.*60.0 # arcseconds. Distance in image from phase centre, i.e. half the FoV.
maxrate = 3e-3 # Hz, typical residual rates on long baselines, value for M82 data
maxdelay = 300e-9 # seconds. typical residual delays on long baselines, value for M82 data

# If true, add also the delay to the image offset. In principle, one should use this
# to get a true measure of the averaging. If true, the delay
# given as max delay above will be converted to a corresponding offset from the phase center,
# and this offset will be added to the offset given. 
add_delayloss = False

# Non user input below
imoffset_rad = (imoffset/3600.0)*(np.pi/180.0) # auto-conversion to radians
c = 2.998e8 # speed of light in m/s
lam = c/nu0 # Wavelength of observing frequency
omega = 2*np.pi/(23*3600.0+56*60) # One rotation in one sidereal day. Rad/sec.

if add_delayloss:
    # Calculate position offset from delay and add this to FoV stated in imoffset
    delayoff_as = 3600.0*(180.0/np.pi)*np.arcsin(maxdelay*c/b)
    print 'Adding delay error of ' + str(round(delayoff_as)) + ' to desired FoV radius of ' + str(imoffset) + '.' 
    imoffset = imoffset + delayoff_as
    print 'New FoV radius required is ' + str(round(imoffset)) + ' arcsec.'
    imoffset_rad = (imoffset/3600.0)*(np.pi/180.0)

# Calculate the loss due to averaging in fourier space, over a phase change dtheta.
def loss(dtheta):
    return (0.5*dtheta)**2/6.0

# Calculate loss due to frequency
def print_freqloss():
    path = bw * b/c
    avgangle = path*2*np.pi*np.sin(imoffset_rad)
    if avgangle > np.pi:
        print 'Warning: Too much averaging in frequency!'
    print 'Loss due to frequency averaging is <' + str(np.round(100*loss(avgangle),2)) + '%.'

# Calculate loss due to time
def print_timeloss():
    path = (b/lam) * time * omega
    avgangle = path * 2*np.pi * np.sin(imoffset_rad)
    if avgangle > np.pi:
        print 'Warning: Too much averaging in time!'
    print 'Loss due to time averaging is <' + str(np.round(100*loss(avgangle),2)) + '%.'

# Calculate loss due to residual rates
def print_rateloss():
    avgangle = maxrate * time * 2*np.pi
    if avgangle > np.pi:
        print 'Warning: Too high rate!!'
    print 'Loss due to rate is  <' + str(np.round(100*loss(avgangle),2)) + '%.'

print_freqloss()
print_timeloss()
print_rateloss()
