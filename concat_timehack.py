#!/usr/bin/env python

import pyrap.tables as pt
import glob, optparse, sys, subprocess, os
import numpy as np
from pyrap.quanta import quantity
from datetime import datetime

usage = "usage: python %prog [options] MStoConcat1 MS2 .."
description="Temporary script to concat data in such a way to eliminate time gaps between the data to enable imaging with AWimager."
vers="1.0"

parser = optparse.OptionParser(usage=usage, version="%prog v{0}".format(vers), description=description)
parser.add_option("-f", "--force", action="store_true", dest="force", default=False ,help="Force concat despite different time intervals [default: %default]")
parser.add_option("-g", "--gap", action="store", type='float', dest="gap", default=120. ,help="Define time gap, in seconds, to include between observations [default: %default]")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default='TIME_HACK.MS',help="Define the name of the output measurement set [default: %default]")
parser.add_option("-w", "--overwrite", action="store_true", dest="ow", default=False,help="Overwrite output MS if it already exists [default: %default]")
(options, args) = parser.parse_args()

def checkintervals(mss):
	"""Checks the intervals of each measurement set"""
	intervals=[]
	print "Measurement Set\t\t\tTime Interval\t\tStart\t\tEnd"
	print "-----------------------------------------------------------------------------------------------"
	for i in mss:
		temp=pt.table(i, ack=False)
		t=temp.getcol('INTERVAL')[0]
		intervals.append(t)
		temp.close()
		#Open table again to get start and end tomes
		temp=pt.table(i+'/OBSERVATION', ack=False)
		tempst=temp.getcell("LOFAR_OBSERVATION_START", 0)
		tempend=temp.getcell("LOFAR_OBSERVATION_END", 0)
		print "{0}\t{1}s\t{2}\t{3}".format(i, t, datetime.utcfromtimestamp(quantity('{0}s'.format(tempst)).to_unix_time()), 
		datetime.utcfromtimestamp(quantity('{0}s'.format(tempend)).to_unix_time()))
	intervals=np.array(intervals)
	uniq=np.unique(intervals)
	if len(uniq) > 1:
		return False, uniq
	else:
		return True, uniq

def concat(sets, outname):
	"""Just the function from concat.py"""
	sys.stdout.write("Concatenating Measurement Sets...")
	sys.stdout.flush()
	try:
		newtable = pt.table(sets, ack=False)
		newtable.sort('TIME').copy(outname, deep = True)
		print "Done!"
		return True
	except:
		print "Error!"
		print "Concat failed most likely channel number differences?"
		return False

def newtimearray(times, cen_times, interval, g):
	newtimes={} #to hold the new times for each unique old time
	newtimes_cen={}
	diff_times=np.unique(times)
	diff_cen_times=np.unique(times)
	newtimes[diff_times[0]]=diff_times[0]
	newtimes_cen[diff_cen_times[0]]=diff_cen_times[0]
	boost=False
	for i in range(1,len(diff_times)):
		if diff_times[i] not in newtimes.keys():
			if diff_times[i]-diff_times[i-1]>interval*1.01:
				boost=True
				newtimes[diff_times[i]]=newtimes[diff_times[i-1]]+interval+g
			else:
				newtimes[diff_times[i]]=newtimes[diff_times[i-1]]+interval
		if diff_cen_times[i] not in newtimes_cen.keys():
			if boost:
				newtimes_cen[diff_cen_times[i]]=newtimes_cen[diff_cen_times[i-1]]+interval+g
			else:
				newtimes_cen[diff_cen_times[i]]=newtimes_cen[diff_cen_times[i-1]]+interval
		boost=False
	# for j in diff_times:
		# print "{0} --> {1}".format(datetime.utcfromtimestamp(quantity('{0}s'.format(j)).to_unix_time()), datetime.utcfromtimestamp(quantity('{0}s'.format(newtimes[j])).to_unix_time()))
	newtimes_list=[newtimes[k] for k in times]
	newcentimes_list=[newtimes_cen[l] for l in cen_times]
	return np.array(newtimes_list), np.array(newcentimes_list), newtimes_list[0], newtimes_list[-1]

oname=options.output
gap=options.gap

#check for measurement sets arguments
if len(args)<1:
	print "You must enter some datasets! - 'python hack.py MS1 MS2 ... MSX'"
	sys.exit()

#check if output already exists
if os.path.isdir(oname):
	if not options.ow:
		print "{0} already exists! Delete or run with '-w' option".format(oname)
		sys.exit()
	else:
		subprocess.call(["rm", "-r", oname])

#check intervals
print "\nChecking intervals of sets {0}...\n".format(", ".join(args))
check, interv=checkintervals(args)
if check:
	inter=interv[0]
else:
	if options.force:
		print "\nForcing concat with different time intervals. Will take first interval as interval to use"
		inter=interv[0]
	else:
		print "\nMeasurement sets intervals are not the same!"
		print "Intervals detected: {0}".format(interval)
		sys.exit()

print "\nInterval = {0}s".format(inter)	

if concat(args, oname):
# if True:
	print "Changing Times on Set..."
	mstochange=pt.table(oname, ack=False, readonly=False)
	amps_time = mstochange.getcol('TIME')
	amps_time_cen = mstochange.getcol('TIME_CENTROID')
	new_time, newtime_cen, start, end=newtimearray(amps_time, amps_time_cen, inter, gap)
	mstochange.putcol('TIME',new_time)
	mstochange.putcol('TIME_CENTROID',newtime_cen)
	mstochange.close()
	mstochange=pt.table(oname+'/OBSERVATION', ack=False, readonly=False)
	mstochange.putcell("LOFAR_OBSERVATION_START", 0, start)
	mstochange.putcell("LOFAR_OBSERVATION_END", 0, start)
	newrange=np.array([start, end])
	mstochange.putcell("TIME_RANGE", 0, newrange)
	mstochange.close()
	inttime=end-start
	print "Start Time: {0}".format(datetime.utcfromtimestamp(quantity('{0}s'.format(start)).to_unix_time()))
	print "NEW End Time: {0}".format(datetime.utcfromtimestamp(quantity('{0}s'.format(end)).to_unix_time()))
	print "Integration Time = {0:.2f}s ({1:.2f} hours)".format(inttime, (inttime)/3600.)
else:
	sys.exit()
