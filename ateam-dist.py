#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
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

# Usage: ./ateam-dist.py RA (00h00m00s) DEC (+/-00d00m00s)
# Example: ./ateam-dist.py 12h30m49.4s +12d23m28s
# calculate the angular distance between an object and ateams/calibrators


import sys, os, re
import numpy as np
from lib_coordinates_mode import *

# A-team
ateam={'CygA':{'ra':299.8679167,'dec':40.7338889},\
'CasA':{'ra':350.8583333,'dec':58.8000000},\
'TauA':{'ra':83.6333333,'dec':22.0144444},\
'VirA':{'ra':187.7058333,'dec':12.3911111}\
}
# Calibrators
cal={'3c48':{'ra':24.4220808,'dec':33.1597594},\
'3c147':{'ra':85.6505746,'dec':49.8520094},\
'3c196':{'ra':123.4001379,'dec':48.2173778},\
'3c286':{'ra':202.7845329,'dec':30.5091550},\
'3c295':{'ra':212.835495,'dec':52.202770},\
'3c380':{'ra':277.3824204,'dec':48.7461556}\
}

def distance(ra1, dec1, ra2, dec2):
    """Return the angular distance between two sources

    Arguments:
    ra1, dec1: decimal degrees values of ra and dec for the first sources
    ra2, dec2: decimal degrees values of ra and dec for the second sources

    Return value:
    angular distance in deg
    """
    ra1=ra1*np.pi/180
    dec1=dec1*np.pi/180
    ra2=ra2*np.pi/180
    dec2=dec2*np.pi/180

    cosdist = np.cos(np.pi/2 - dec1)*np.cos(np.pi/2 - dec2)+np.sin(np.pi/2 - dec1)*np.sin(np.pi/2-dec2)*np.cos(ra1-ra2)
    return np.arccos(cosdist)*180/np.pi

def printdist(name,d):
	print(name, str(d), "deg", end=' ')
	if d<25: print('*CLOSE*')
	else: print('')


rah, ram, ras = re.sub(r'[h|m|s]', ' ', sys.argv[1]).split()
decd, decm, decs = re.sub(r'[d|m|s]', ' ', sys.argv[2]).split()
objra = hmstora(rah, ram, ras)
objdec = dmstodec(decd, decm, decs)

print("Distance form A-team:")
for name in ateam:
	d=distance(objra,objdec,ateam[name]['ra'],ateam[name]['dec'])
	printdist(name, d)

print("Distance form Calibrators:")
for name in cal:
	d=distance(objra,objdec,cal[name]['ra'],cal[name]['dec'])
	printdist(name, d)
