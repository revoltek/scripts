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

import sys
import pydal as dal
from pylab import *
from numpy import *
from scipy import *

# check usage
if len(sys.argv) < 4 or len(sys.argv) > 10:
        print("Usage:")
        print("\tbaseline.py <file> <antenna1> <antenna2> " + \
              "[sub-band] ['channel' or 'time'] [channel/time range] " + \
              "['p'hase or 'a'mplitude] [polarization (0-3,4-6)] [plot range xmin,xmax,ymin,ymax (no spaces)]")
        print("\t<> required")
        print("\t E.g., baseline.py /lifs006/SB3.MS 1 10 0 channel 3:10 a 6 4.69e9,4.7e9,0,0.001")
        print("\tSingle or range of channels and times for averaging (in relative bin coordinates; ranges colon-delimited);")
        print("\tPolarization 0-3 plots individually, 4 plots xx and xy, 5 plots xx and yx, 6 plots xx and yy.")
        print("\t-1 means plot all.  Default plots the closure phase of the xx for all channels as function of time.")
        print("Note that antenna numbers are zero-based and must be given in order from lowest to highest.")
        print("")
        sys.exit(1)

# other customizations
data_name = "DATA"

# set default values
if len(sys.argv) < 5:  subband = str(0)
else: subband = sys.argv[4]
if len(sys.argv) < 6:  quantity_plot = 'channel'
else: quantity_plot = sys.argv[5]
if len(sys.argv) < 7:  range_plot = -1
elif len((sys.argv[6]).split(':')) > 1:
        range_plot = list(range(int((sys.argv[6]).split(':')[0]),int((sys.argv[6]).split(':')[1])))
else: range_plot = [int(sys.argv[6])]
if len(sys.argv) < 8:  data_plot = 'a'
else:	data_plot = sys.argv[7]
if len(sys.argv) < 9:  pol = 0
else:	pol = int(sys.argv[8])
if len(sys.argv) < 10: axis_range = 0
else:
        axis_range = []
        for i in [0,1,2,3]:  axis_range.append(float(sys.argv[9].split(',')[i]))

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

# open table
tablename = "MAIN";
msds.setFilter( "TIME," + data_name, \
        "ANTENNA1 = " + sys.argv[2] + " AND ANTENNA2 = " + sys.argv[3] + \
        " AND DATA_DESC_ID = " + subband )
maintable = msds.openTable( tablename );

# get times
time_col = maintable.getColumn("TIME")
time = time_col.data()
time = time/(24*3600)    # convert from MJD in seconds to days

print('Start of data (MJD from 1Oct2007): ' + str(time[0]-54374))
print('Fractional day is ' + str(24*(time[0]-int(time[0]))) + ' hrs')
print('Total length of data (integrations): ' + str(len(time)) + ', (hrs): ' + str((time[-1]-time[0])*24))
#event_time = 54379+(06+(41/60.))/24.      # 071006, obs L3907, el=-50
#event_time = 54381+(21+(55/60.))/24.      # 071008, obs L3914, el=+8       ## look at SB1, ch34 for best rfi (38.6009216 MHz)
#event_time = 54383+(03+(41/60.))/24.      # 071010A, obs L3917, el=-60
#event_time = 54383+(20+(45/60.))/24.      # 071010B, obs L3940, el=+8, stronger
#event_time = 54383+(22+(20/60.))/24.      # 071010C, obs L3940, el=+72
#event_time = 54386+(12+(9/60.))/24.       # 071013, obs L3980, el=+37
event_time = 0

# get data
data_col = maintable.getColumn(data_name)
data = data_col.data()
nchannels = data.shape[1] # second element of the data shape i.e. (nrows,256,4)

# if the optional channel argument is present
#  plot for this channel/time
if (range_plot != -1):
        if quantity_plot == 'channel':
                data_reduce = add.reduce(array=data[:,range_plot,:],axis=1)/len(range_plot)
                # plot data of given data vs. time
                if data_plot == 'a':
                        current_value = hypot((data_reduce[:,pol]).real,(data_reduce[:,pol]).imag)
                        if pol2:  current_value2 = hypot((data_reduce[:,pol2]).real,(data_reduce[:,pol2]).imag)
                elif data_plot == 'p':
                        current_value = arctan2((data_reduce[:,pol]).imag,(data_reduce[:,pol]).real)%2*pi-pi
                        if pol2:  current_value2 = arctan2((data_reduce[:,pol2]).imag,(data_reduce[:,pol2]).real)%2*pi-pi
                plot( time, current_value, "," )
                if pol2:  plot( time, current_value2, "," )
                if axis_range:  axis(axis_range)
                title("Time vs. Amplitude, Baseline " + \
                      sys.argv[2] + '-' + sys.argv[3] + ", Sub-band(" + subband +
                      ") " + " Chan0(" + str(range_plot[0]) + "), Nchan(" + str(len(range_plot)) + ")\n" + sys.argv[1] )
                xlabel("Time (MJD)")
                if event_time:  plot( event_time*ones(2), array([min(current_value),max(current_value)]))

        elif quantity_plot == 'time':
                data_reduce = add.reduce(array=data[range_plot,:,:],axis=0)/len(range_plot)
                # plot intensity of given data vs. channel
                if data_plot == 'a':
                        current_value = hypot((data_reduce[:,pol]).real,(data_reduce[:,pol]).imag)
                        if pol2:  current_value2 = hypot((data_reduce[:,pol2]).real,(data_reduce[:,pol2]).imag)
                elif data_plot == 'p':
                        current_value = arctan2((data_reduce[:,pol]).imag,(data_reduce[:,pol]).real)%2*pi-pi
                        if pol2: current_value2 = arctan2((data_reduce[:,pol2]).imag,(data_reduce[:,pol2]).real)%2*pi-pi
                plot( current_value, "," )
                if pol2: plot( current_value2, "," )
                if axis_range:  axis(axis_range)
                title("Channel vs. Amplitude, Baseline " + sys.argv[2] + '-' + sys.argv[3] + ", Sub-band(" + subband + "), " + \
                      " Time(" + str(time[range_plot[0]]) + " MJD) \n dtime(" + str(len(range_plot)) + " ints), " + sys.argv[1] )
                xlabel("Channel")

# otherwise, plot all channels/times
else:
        if quantity_plot == 'channel':
                current_value_min = 0.
                current_value_max = 0.
                # plot intensity of each channel vs. time
                for channel in range( nchannels ):
                        if data_plot == 'a':
                                current_value = hypot((data[:,channel,pol]).real,(data[:,channel,pol]).imag)
                                if pol2: current_value2 = hypot((data[:,channel,pol2]).real,(data[:,channel,pol2]).imag)
                                current_value_min = min(current_value_min,min(current_value))
                                current_value_max = max(current_value_max,max(current_value))
                        elif data_plot == 'p':
                                current_value = arctan2((data[:,channel,pol]).imag),array((data[:,channel,pol]).real)%2*pi-pi
                                if pol2: current_value2 = arctan2((data[:,channel,pol2]).imag,(data[:,channel,pol2]).real)%2*pi-pi
                        plot( time, current_value, "," )
                        if pol2: plot ( time, current_value2, ",")
                        if axis_range:  axis(axis_range)
                        title("Time vs. Amplitude, Baseline " + \
                              sys.argv[2] + '-' + sys.argv[3] + ", Sub-band(" + subband +
                              ") " + str(data.shape[1]) + " channels" + '\n' + sys.argv[1] )
                xlabel("Time (MJD)")
                if event_time:  plot( event_time*ones(2), array([current_value_min,current_value_max]))

        if quantity_plot == 'time':
                # plot intensity at each time vs. channel
                for t in range( len(time) ):
                        if data_plot == 'a':
                                current_value = hypot((data[t,:,pol]).real,(data[t,:,pol]).imag)
                                if pol2: current_value2 = hypot((data[t,:,pol2]).real,(data[t,:,pol2]).imag)
                        elif data_plot == 'p':
                                current_value = arctan2((data[t,:,pol]).imag,(data[t,:,pol]).real)%2*pi-pi
                                if pol2: current_value2 = arctan2((data[t,:,pol2]).imag,(data[t,:,pol2]).real)%2*pi-pi
                        plot( current_value, "," )
                        if pol2: plot( current_value2, ",")
                        if axis_range:  axis(axis_range)
                        title("Time vs. Amplitude, Baseline " + \
                              sys.argv[2] + '-' + sys.argv[3] + ", Sub-band(" + subband +
                              ") " + str(len(time)) + " times" + '\n' + sys.argv[1] )
                xlabel("Channel")

if data_plot == 'a':
        ylabel("Intensity")
elif data_plot == 'p':
        ylabel("Phase (rad)")
        
show()
