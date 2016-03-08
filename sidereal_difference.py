#! /usr/bin/env python

import sys
import pydal as dal
from pylab import *
from numpy import *
from scipy import *

# check usage
if len(sys.argv) < 4 or len(sys.argv) > 12:
        print "Usage:"
        print "\tsidereal_difference.py <file> x 3 <antenna1> <antenna2> " + \
	      "[sub-band [channel begin:channel end [start time " + \
              "[ time adv [polarization]]]]]"
        print "\t<> required"
        print ""
        sys.exit(1)

# other customizations
data_name = "DATA"
timecut = (23+56/60.)/24    # 1 sidereal day, 23h56m

# set default values
print "File 1: ", sys.argv[1]
print "File 2: ", sys.argv[2]
print "File 3: ", sys.argv[3]
print "Antenna 1: ", sys.argv[4]
print "Antenna 2: ", sys.argv[5]
if len(sys.argv) < 7:  subband = str(0)
else:
	subband = sys.argv[6]
	print "Sub-band: ", sys.argv[6]
if len(sys.argv) < 8:  range_plot = [0]
elif len((sys.argv[7]).split(':')) > 1:
        range_plot = range(int((sys.argv[7]).split(':')[0]),int((sys.argv[7]).split(':')[1])+1)
	print "Channels: ", sys.argv[7]
else:
	range_plot = [int(sys.argv[7])]
	print "Channels: ", sys.argv[7]
if len(sys.argv) < 9:  start_time = 0
else:
	start_time = int(sys.argv[8])
	print "Start time: ", sys.argv[8]
if len(sys.argv) < 10:  time_adv = 0
else:	
	time_adv = int(sys.argv[9])
	print "Time adv: ", sys.argv[9]
if len(sys.argv) < 11:
	pol = 0
else:	
	pol = int(sys.argv[10])
	print "Polarisation: ", sys.argv[10]

# open files
msds1 = dal.dalDataset()
if ( True != msds1.open(sys.argv[1]) ):
        print "ERROR: Could not open file: " + sys.argv[1]
        print "       Please check the file and try again."
        sys.exit(1)
msds2 = dal.dalDataset()
if ( True != msds2.open(sys.argv[2]) ):
        print "ERROR: Could not open file: " + sys.argv[2]
        print "       Please check the file and try again."
        sys.exit(1)
msds3 = dal.dalDataset()
if ( True != msds3.open(sys.argv[3]) ):
        print "ERROR: Could not open file: " + sys.argv[3]
        print "       Please check the file and try again."
        sys.exit(1)


# open tables
tablename = "MAIN";
msds1.setFilter( "TIME," + data_name, \
        "ANTENNA1 = " + sys.argv[4] + " AND ANTENNA2 = " + sys.argv[5] + \
        " AND DATA_DESC_ID = " + subband )
maintable1 = msds1.openTable( tablename );
msds2.setFilter( "TIME," + data_name, \
        "ANTENNA1 = " + sys.argv[4] + " AND ANTENNA2 = " + sys.argv[5] + \
        " AND DATA_DESC_ID = " + subband )
maintable2 = msds2.openTable( tablename );
msds3.setFilter( "TIME," + data_name, \
        "ANTENNA1 = " + sys.argv[4] + " AND ANTENNA2 = " + sys.argv[5] + \
        " AND DATA_DESC_ID = " + subband )
maintable3 = msds3.openTable( tablename );

# get times
print "Getting time..."
time_col = maintable1.getColumn("TIME")
time = time_col.data()
time = array(time[start_time:]/(24*3600))    # convert from MJD in seconds to days

# get datasets
print "Getting data 1..."
data_col = maintable1.getColumn(data_name)
data1 = array(data_col.data()[start_time:,range_plot,:])
data1_reduce = add.reduce(array=data1[:,:,:],axis=1)/len(range_plot)
print "Getting data 2..."
data_col = maintable2.getColumn(data_name)
data2 = array(data_col.data()[start_time:,range_plot,:])
data2_reduce = add.reduce(array=data2[:,:,:],axis=1)/len(range_plot)
print "Getting data 3..."
data_col = maintable3.getColumn(data_name)
data3 = array(data_col.data()[start_time:,range_plot,:])
data3_reduce = add.reduce(array=data3[:,:,:],axis=1)/len(range_plot)

# get filtered data and time (new way;  assumes all time bins are same size)
indices_day1 = where((time-time[0]) <= timecut)[0]
time_day1 = time[indices_day1]
data1_day1 = data1_reduce[indices_day1,:]
data2_day1 = data2_reduce[indices_day1,:]
data3_day1 = data3_reduce[indices_day1,:]
time_day2 = time[indices_day1+len(indices_day1)]
data1_day2 = data1_reduce[indices_day1+len(indices_day1),:]
data2_day2 = data2_reduce[indices_day1+len(indices_day1),:]
data3_day2 = data3_reduce[indices_day1+len(indices_day1),:]

print 'Start of data (MJD from 1Jan2000): ' + str(time[0]-51544)
print 'Start of data (MJD from 17Oct2007): ' + str(time[0]-54390) + '.  Fractional day is ' + str(24*(time[0]-int(time[0]))) + ' hrs'
print 'Total length of data (integrations): ' + str(len(time)+start_time) + ', Day 1 (hrs): ' + str((time_day1[-1]-time_day1[0])*24) + ', Day 2 (hrs): ' + str((time_day2[-1]-time_day2[0])*24)

# time adveraging
if time_adv != 0:
	print "Time adveraging..."
	num = len(time_day1)/time_adv
	for i in range(num):
		data1_day1[i,:] = mean(data1_day1[i*time_adv:(i+1)*time_adv],axis=0)
		data2_day1[i,:] = mean(data2_day1[i*time_adv:(i+1)*time_adv],axis=0)
		data3_day1[i,:] = mean(data3_day1[i*time_adv:(i+1)*time_adv],axis=0)
		data1_day2[i,:] = mean(data1_day2[i*time_adv:(i+1)*time_adv],axis=0)
		data2_day2[i,:] = mean(data2_day2[i*time_adv:(i+1)*time_adv],axis=0)
		data3_day2[i,:] = mean(data3_day2[i*time_adv:(i+1)*time_adv],axis=0)
		time_day1[i] = time_day1[i*time_adv+time_adv/2]
		time_day2[i] = time_day2[i*time_adv+time_adv/2]
	data1_day1 = data1_day1[:num,:]
	data2_day1 = data2_day1[:num,:]
	data3_day1 = data3_day1[:num,:]
	data1_day2 = data1_day2[:num,:]
	data2_day2 = data2_day2[:num,:]
	data3_day2 = data3_day2[:num,:]
	time_day1 = time_day1[:num]
	time_day2 = time_day2[:num]

print "Plotting..."
# set plot label and title
titlestring = "Amplitude and Division, 118 MHz - 149 MHz - 180 MHz - File1:"+sys.argv[1]
ylabelstring = "Intensity"

current_value1_day1 = hypot((data1_day1[:,pol]).real,(data1_day1[:,pol]).imag)
current_value2_day1 = hypot((data2_day1[:,pol]).real,(data2_day1[:,pol]).imag)
current_value3_day1 = hypot((data3_day1[:,pol]).real,(data3_day1[:,pol]).imag)
current_value1_day2 = hypot((data1_day2[:,pol]).real,(data1_day2[:,pol]).imag)
current_value2_day2 = hypot((data2_day2[:,pol]).real,(data2_day2[:,pol]).imag)
current_value3_day2 = hypot((data3_day2[:,pol]).real,(data3_day2[:,pol]).imag)

params = {
	'axes.labelsize': 10,
	'text.fontsize': 10,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'figure.figsize': [12,8]}
rcParams.update(params)

figure(1)

subplot(121)
plot( time_day1, current_value1_day1 + 0.01, ",", label="Day 1 + 0.010 @ 118 MHz")
plot( time_day2-timecut, current_value1_day2 + 0.008, ",", label="Day 2 + 0.008 @ 118 MHz" )
plot( time_day1, current_value2_day1 + 0.006, ",", label="Day 1 + 0.006 @ 149 MHz")
plot( time_day2-timecut, current_value2_day2 + 0.004, ",", label="Day 2  + 0.004 @ 149 MHz" )
plot( time_day1, current_value3_day1 + 0.002, ",", label="Day 1 + 0.002 @ 180 MHz")
plot( time_day2-timecut, current_value3_day2, ",", label="Day 2 @ 180 MHz" )
title(titlestring)
ylabel(ylabelstring)
axis([time_day1[0], time_day1[-1], 0, 0.016])
leg = legend(numpoints=10)
frame = leg.get_frame()
frame.set_alpha(0.7) # make it semi-transparent

subplot(122)
plot( time_day1, (current_value1_day2-current_value1_day1)/(0.5*(current_value1_day2+current_value1_day1)) + 0.1, "-", label="+0.2 @ 118 MHz")
plot( time_day1, (current_value2_day2-current_value2_day1)/(0.5*(current_value2_day2+current_value2_day1)), "-", label="@ 149 MHz" )
plot( time_day1, (current_value3_day2-current_value3_day1)/(0.5*(current_value3_day2+current_value3_day1)) - 0.1, "-", label="-0.2 @ 180 MHz" )
ylabel("Day1-Day2/mean")
xlabel("Time (MJD)")
axis([time_day1[0], time_day1[-1], -0.2, 0.2])
leg = legend(numpoints=10)
frame = leg.get_frame()
frame.set_alpha(0.7) # make it semi-transparent

outputfile = "out.png"
print "\nSaving output file as: ", outputfile
savefig( outputfile )
show()

#print "Computing displacement..."
#displacement = (current_value_day1[:] - current_value_day2[:]) / ((current_value_day1[:] + current_value_day2[:])/2)
#print "Mean displacement: ", mean(displacement), "%"
