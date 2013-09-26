#!/usr/bin/python
# Create a fame skymodel, useful for tests
import numpy as np

# number of point sourcer
Nsources = 10
# rmax of the spiral in deg
Rmax = 5.
# number of angles
Nang = 8.

# Center of the spiral
#virgo
#12:30:49.420000
rah=12.
ram=30.
ras=49.42
#+012.23.28.000000
decd=12.
decm=23.
decs=28.

# print necessary lines
print "# (Name, Type, Ra, Dec, I, Q, U, V, ReferenceFrequency='74e6', SpectralIndex='[]') = format"
print "C0, POINT, 12:30:49.42, +12.23.28.00, 10.0, 0.0, 0.0, 0.0, 136.0e+06, [0]"
print "center, POINT, 12:30:49.42, +12.23.28.00, 0.0, 0.0, 0.0, 0.0, 136.0e+06, [0]"

def hmstora(rah,ram,ras):
    """Convert RA in hours, minutes, seconds format to decimal
    degrees format.

    Keyword arguments:
    rah,ram,ras -- RA values (h,m,s)

    Return value:
    radegs -- RA in decimal degrees

    """
    hrs = (float(rah)+(float(ram)/60)+(float(ras)/3600.0)) % 24

    return 15*hrs
# Converts an hms format RA to decimal degrees

def dmstodec(decd,decm,decs):
    """Convert Dec in degrees, minutes, seconds format to decimal
    degrees format.

    Keyword arguments:
    decd,decm,decs -- list of Dec values (d,m,s)

    Return value:
    decdegs -- Dec in decimal degrees

    """
    if decd < 0:
        decm = -1*decm
        decs = -1*decs

    decdegs = float(decd)+(float(decm)/60)+(float(decs)/3600.0)

    if abs(decdegs) > 90:
        raise ValueError

    return decdegs
# Converts a dms format Dec to decimal degrees  

def sec2hms(seconds):
    """Seconds to hours, minutes, seconds"""
    hours, seconds = divmod(seconds, 60**2)
    minutes, seconds = divmod(seconds, 60)
    return (int(hours), int(minutes), seconds)


def ratohms(radegs):
    """Convert RA in decimal degrees format to hours, minutes,
    seconds format.

    Keyword arguments:
    radegs -- RA in degrees format

    Return value:
    ra -- tuple of 3 values, [hours,minutes,seconds]

    """

    radegs %= 360
    raseconds = radegs * 3600 / 15.0
    return sec2hms(raseconds)

def dectodms(decdegs):
    """Convert Declination in decimal degrees format to hours, minutes,
    seconds format.

    Keyword arguments:
    decdegs -- Dec. in degrees format

    Return value:
    dec -- list of 3 values, [degrees,minutes,seconds]

    """

    if abs(decdegs) > 90:
        raise ValueError
    decd = int(decdegs)
    decm = int((decdegs-decd)*60)
    decs = (((decdegs-decd)*60)-decm)*60
    if decd < 0:
        decm = -1*decm
        decs = -1*decs
    dec = (decd,decm,decs)
    return dec

radeg = hmstora(rah,ram,ras)
decdeg = dmstodec(decd,decm,decs)

Rinc = Rmax / float(Nsources)
for i in xrange(Nsources):
	# compute radii values
	R = (i+1)*Rinc
	# compute angles value
	ang = i*(2*np.pi)/Nang
	# convert from polar to cartesian
	X=R*np.cos(ang)
	Y=R*np.sin(ang)
	# convert to RA,DEC
	(rah,ram,ras) = ratohms(radeg+X)
	(decd,decm,decs) = dectodms(decdeg+Y)
    	ra = str(rah)+":"+str(ram)+":"+str(ras)
    	dec = "+"+str(decd)+"."+str(decm)+"."+str(decs)
    	print "C"+str(i+1)+", POINT, "+str(ra)+", "+str(dec)+", 1.0, 0.0, 0.0, 0.0, 136.0e+06, [0]"
