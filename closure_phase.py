#! /usr/bin/env python

import sys
import pydal as dal
from pylab import *
from scipy import *

# check usage
if len(sys.argv) < 5 or len(sys.argv) > 9:
        print("Usage:")
        print("\tclosure_phase.py <file> <antenna1> <antenna2> <antenna3> " + \
              "[sub-band] ['channel' or 'time'] [channel/time range] " + \
              "[polarization (0-6)]")
        print("\t<> required")
        print("\tE.g., closure_phase.py /lifs003/SB3.MS 1 5 10 0 channel 100:110 6")
        print("\tSingle or range of channels and times for averaging (in relative bin coordinates; ranges colon-delimited);")
        print("\tPolarization 0-3 plots them individually, 4 plots xx and xy, 5 plots xx and yx, 6 plots xx and yy")
        print("\t-1 means plot all.  Default plots the closure phase of the xx for all channels as function of time for subband 0.")
        print("\tNote that antenna numbers are zero-based and must be given in order from lowest to highest.")
        print("")
        sys.exit(1)

# other customizations
data_name = "DATA"

# set default values
if len(sys.argv) < 6:  subband = str(0)
else: subband = sys.argv[5]
if len(sys.argv) < 7:  quantity_plot = 'channel'
else: quantity_plot = sys.argv[6]
if len(sys.argv) < 8:  range_plot = [0]
elif len((sys.argv[7]).split(':')) > 1:
        range_plot = list(range(int((sys.argv[7]).split(':')[0]),int((sys.argv[7]).split(':')[1])))
else: range_plot = [int(sys.argv[7])]
if len(sys.argv) < 9:  pol = 0
else:	pol = int(sys.argv[8])
if pol > 3:
        pol2 = pol-3
        pol = 0
        print("plotting pols %i and %i" % (pol, pol2))
else:
        pol2 = 0
        print("plotting pol %i" % (pol))

# open file
msds= dal.dalDataset()
if ( True != msds.open(sys.argv[1]) ):
        print("ERROR: Could not open file: " + sys.argv[1])
        print("       Please check the file and try again.")
        sys.exit(1)

# open tables
tablename = "MAIN";
msds.setFilter( "TIME," + data_name, \
        "ANTENNA1 = " + sys.argv[2] + " AND ANTENNA2 = " + sys.argv[3] + \
        " AND DATA_DESC_ID = " + subband )
table12 = msds.openTable( tablename );

msds.setFilter( "TIME," + data_name, \
        "ANTENNA1 = " + sys.argv[3] + " AND ANTENNA2 = " + sys.argv[4] + \
        " AND DATA_DESC_ID = " + subband )
table23 = msds.openTable( tablename );

msds.setFilter( "TIME," + data_name, \
        "ANTENNA1 = " + sys.argv[2] + " AND ANTENNA2 = " + sys.argv[4] + \
        " AND DATA_DESC_ID = " + subband )
table13 = msds.openTable( tablename );

# get times
time_col = table12.getColumn("TIME")
time = time_col.data()
time = time/(24*3600)    # convert from MJD in seconds to days

# get data
data_col12 = table12.getColumn(data_name)
data12 = data_col12.data()
data_col23 = table23.getColumn(data_name)
data23 = data_col23.data()
data_col13 = table13.getColumn(data_name)
data13 = data_col13.data()
nchannels = data12.shape[1] # second element of the data shape i.e. (nrows,256,4)

# calculate phases
if quantity_plot == 'channel':       
        data12_reduce = add.reduce(array=data12[:,range_plot,:],axis=1)/len(range_plot)
        data23_reduce = add.reduce(array=data23[:,range_plot,:],axis=1)/len(range_plot)
        data13_reduce = add.reduce(array=data13[:,range_plot,:],axis=1)/len(range_plot)
elif quantity_plot == 'time':
        data12_reduce = add.reduce(array=data12[range_plot,:,:],axis=0)/len(range_plot)
        data23_reduce = add.reduce(array=data23[range_plot,:,:],axis=0)/len(range_plot)
        data13_reduce = add.reduce(array=data13[range_plot,:,:],axis=0)/len(range_plot)

phase12 = arctan2((data12_reduce[:,pol]).imag,(data12_reduce[:,pol]).real)
phase23 = arctan2((data23_reduce[:,pol]).imag,(data23_reduce[:,pol]).real)
phase13 = arctan2((data13_reduce[:,pol]).imag,(data13_reduce[:,pol]).real)
if pol2:
        phase12_2 = arctan2((data12_reduce[:,pol2]).imag,(data12_reduce[:,pol2]).real)
        phase23_2 = arctan2((data23_reduce[:,pol2]).imag,(data23_reduce[:,pol2]).real)
        phase13_2 = arctan2((data13_reduce[:,pol2]).imag,(data13_reduce[:,pol2]).real)

closure = (phase12 + phase23 - phase13)%(2*pi)
if pol2: closure_2 = (phase12_2 + phase23_2 - phase13_2)%(2*pi)
closure[numpy.where(closure > pi)] = closure[numpy.where(closure > pi)] - 2*pi
if pol2: closure_2[where(closure_2 > pi)] = closure_2[where(closure_2 > pi)] - 2*pi
        
# if the optional channel argument is present
#  plot for this channel/time
if (range_plot != -1):
        if quantity_plot == 'channel':
		subplot(411)
                title("Time vs. Closure Phase, Antennas " + \
                      sys.argv[2] + '-' + sys.argv[3] + '-' + sys.argv[4] + ", Sub-band(" + subband + ") " + sys.argv[1] )	
		plot( time, phase12, "," )
                if pol2: plot( time, phase12_2, "," )
		subplot(412)
		plot( time, phase23, "," )
                if pol2: plot( time, phase23_2, "," )
		subplot(413)
		plot( time, phase13, "," )
                if pol2: plot( time, phase13_2, "," )
                # plot data of given data vs. time
		subplot(414)
                plot( time, closure, "," )
                if pol2: plot( time, closure_2, "," )	
                xlabel("Time (MJD)")

        elif quantity_plot == 'time':
                # plot intensity of given data vs. channel
                plot( closure, "," )
                if pol2: plot( closure_2, "," )
                title("Channel vs. Closure Phase, Antennas " + \
                      sys.argv[2] + '-' + sys.argv[3] + '-' + sys.argv[4] + ", Sub-band(" + subband +
                      ") " + " Time(" + str(time[range_plot[0]]) + " MJD)\n dtime(" + str(len(range_plot)) + ") " + sys.argv[1] )
                xlabel("Channel")

## This is broken after adding the data averaging ##
#
# otherwise, plot all channels/times
#else:
#	if quantity_plot == 'channel':
#		# plot intensity of each channel vs. time
#		for channel in range( nchannels ):
#			plot( closure[channel,:], "," )
#			if pol2:  plot( closure_2[channel,:], "," )
#			title("Channel vs. Closure Phase, Antennas " + \
#			      sys.argv[2] + '-' + sys.argv[3] + '-' + sys.argv[4] + ", Sub-band(" + subband +
#			      ") " + " Channel(" + str(range_plot) + ")\n" + sys.argv[1] )
#		xlabel("Time (MJD)")
#	if quantity_plot == 'time':
#		# plot intensity at each time vs. channel
#		for t in range( len(time) ):
#			plot( closure[:,t], "," )
#			if pol2:  plot( closure_2[:,t], "," )
#			title("Time vs. Closure Phase, Antennas " + \
#			      sys.argv[2] + '-' + sys.argv[3] + '-' + sys.argv[4] + ", Sub-band(" + subband +
#			      ") " + str(len(time)) + " times" + '\n' + sys.argv[1] )
#		xlabel("Channel")

ylabel("Closure phase (rad)")


outputfile = "out.png"
savefig( outputfile )
show()
