#!/usr/bin/python

from datetime import datetime as dt
from astrotime import *
import sys

vlalong = 107.6184

if len(sys.argv) == 0:
    print("Usare: vla_utc2lst.py 4-11-2013 22:07:05")
    print("date: dd-mm-yyyy -- time: hh:mm:ss")

day, month, yr = sys.argv[1].split('-')
h, m, s = sys.argv[2].split(':')

t = dt(int(yr),int(month),int(day),int(h),int(m),int(s))

gmst = utcDatetime2gmst(t)
print(gmst2lst(vlalong, hour=gmst.hour, minute=gmst.minute, second=gmst.second))
