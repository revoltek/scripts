#!/usr/bin/env python

"""
uvplot.py

This script is intended to provide a way to quickly plot visibilities
(any of amplitude, phase, real part, or imaginary part can be plotted)
vs time/chan for LOFAR datasets. It uses pyrap and ppgplot (numpy version).
If you plot channel number or frequency on the x-axis, the data will
be averaged in time. If you plot time on the x-axis, the data will be
averaged in frequency.

If you use postscript device output, you may want to use 'embiggen.csh'
(located at /home/heald/bin/embiggen.csh) which makes the pointsizes bigger.

Written by George Heald, modified by Oscar Martinez
v1.0 completed 1/4/2010

version 1.1    2/4/2010: fix specification of axis limits when all data equal
                         add debug mode option to inputs
version 1.2    3/4/2010: add option to set axis limits
                         add ability to use flag information (if requested)
version 1.3    6/4/2010: fix subplot behavior
                         fix behavior when all data points flagged on a baseline
                         upgrade from numarray to numpy
                         make pyrap commands a bit smarter
                         add option to show autocorrelations
version 1.4    7/4/2010: add ability to select time,chan ranges, or by pol
                         add 'query' option to get quick info
version 1.5   11/4/2010: improve locations of XX,XY,YX,YY strings in plot title
                         improve printing of axis labels (*)
                         add option to unwrap phase
version 1.6   14/4/2010: add option to convert to stokes parameters
version 1.7   22/4/2010: add antenna selection option
version 1.8   12/5/2010: make plotting with antenna selection faster
version 1.9    2/7/2010: allow for use of old pyrap_tables module
version 1.10  14/7/2010: add gui frontend to interactively set options
                         add text line to plot indicating the input MS
version 1.11  16/7/2010: make gui persistent for re-use
                         fix some small gui loading issues
                         get rid of pyrap_tables which doesn't work
version 1.12  20/7/2010: add hh mm ss labeling for time axis
                         a bit of cleanup
                         correctly print number of points plotted per baseline
version 1.13  4/11/2010: fix slicing in time and frequency dimensions
                         improve axis limit behavior a bit
version 1.14  5/11/2010: fix ambiguity in new slicing behavior, allow negatives
                         add output of reference frequency information
version 1.15  8/11/2010: major update to use masked arrays for flagging
version 1.16 15/11/2010: avoid numpy.ma.average bug(?): imag part set to zero
                         allow plotting of real or imaginary parts
version 1.17 15/11/2010: add statistics monitoring: mean and standard deviation
                         of shown samples
                         add option of displaying XX-YY
version 1.18 22/11/2010: tune xxyydiff acquisition (proper treatment of flag masks)                        
version 1.19 17/01/2011: Possibility to select CORRECTED_DATA-MODEL_DATA as 
                         column to plot
                         Write the name of the inputMS in the title of plot
version 1.20 04/02/2011: Add in the title of the figures a data description
                         (C, D, M or C-M)
version 1.21 10/02/2011: Improve antenna selection method
version 1.22 15/02/2011: Allow user to stop plotting partway, nondestructively
version 1.23 08/03/2011: Plot XY.YX* (adding operation options)
version 1.24 15/03/2011: In operations the phase is always within +-PI
version 1.25 13/05/2011: Allow to select arbitrary column names as well as + and
                         - between two arbitrary columns
version 1.26 02/09/2011: Add hour angle for x-axis possibilities

To do:
- Implement time slices.
"""

import optparse
import signal
import sys
import ppgplot
# Note, for documentation on ppgplot, see
# http://www.astro.rug.nl/~breddels/python/ppgplot/ppgplot-doc.txt
import numpy.ma
from Tkinter import *
import tkFileDialog

try:
        import pyrap.tables as pt
except ImportError:
        print "Error: The pyrap tables module is not available."
        print "Perhaps you need to first type \'use LofIm\'?"
        exit()


version_string = 'v1.26, 2 September 2011\nWritten by George Heald, modified by Oscar Martinez'
print 'uvplot.py',version_string
print ''

keepPlotting = True

def main(options):

        global keepPlotting
        keepPlotting = True
        debug = options.debug
        inputMS = options.inms
        if inputMS == '':
                print 'Error: You must specify a MS name.'
                print '       Use "uvplot.py -h" to get help.'
                return
            
        if inputMS.endswith('/'):
            inputMS = inputMS[:-1]
        inputMSbasename = inputMS.split('/')[-1]
        if inputMSbasename == '':
            # The user has not specified the full path of the MS
            inputMSbasename = inputMS
        
        device = options.device
        if device=='?':
                ppgplot.pgldev()
                return
        xaxis = options.xaxis
        if xaxis == 'ha':
            print 'Adding derived columns to allow plotting hour angle...'
            try:
                pt.addDerivedMSCal(inputMS)
            except:
                print 'Failed, trying to remove and add columns...'
                try:
                    pt.removeDerivedMSCal(inputMS)
                    pt.addDerivedMSCal(inputMS)
                except:
                    print 'That failed too... plotting HA seems to not be possible.'
                    return
        yaxis = options.yaxis
        column = options.column
        nx, ny = options.nxy.split(',')
        axlimits = options.axlimits.split(',')
        if len(axlimits) == 4:
                xmin,xmax,ymin,ymax = axlimits
        else:
                print 'Error: You must specify four axis limits'
                return
        showFlags = options.flag
        flagCol = options.colflag
        showAutocorr = options.autocorr
        showStats = options.statistics
        timeslots = options.timeslots.split(',')
        if len(timeslots) != 2:
                print 'Error: Timeslots format is start,end'
                return
        for i in range(len(timeslots)): timeslots[i] = int(timeslots[i])
        antToPlotSpl = options.antennas.split(',')
        antToPlot = []
        for i in range(len(antToPlotSpl)):
                tmpspl = antToPlotSpl[i].split('..')
                if len(tmpspl) == 1:
                        antToPlot.append(int(antToPlotSpl[i]))
                elif len(tmpspl) == 2:
                        for j in range(int(tmpspl[0]),int(tmpspl[1])+1):
                                antToPlot.append(j)
                else:
                        print 'Error: Could not understand antenna list.'
                        return
        polarizations = options.polar.split(',')
        for i in range(len(polarizations)):
                polarizations[i] = int(polarizations[i])
        
        convertStokes = options.stokes        
        
        operation = options.operation
        if operation != '':
            operation = int(operation)
            if convertStokes:
                print 'Error: Stokes conversion is not compatible with special operations'
                return
        
        channels = options.channels.split(',')
        if len(channels) != 2:
                print 'Error: Channels format is start,end'
                return
        for i in range(len(channels)): channels[i] = int(channels[i])
        if channels[1] == -1:
                channels[1] = None # last element even if there is only one
        else:
                channels[1] += 1
        queryMode = options.query
        doUnwrap = options.wrap


        if not queryMode:
                # open the graphics device, use the right number of panels
                ppgplot.pgbeg(device, int(nx), int(ny))
                # set the font size
                ppgplot.pgsch(1.5)
                ppgplot.pgvstd()

        # open the main table and print some info about the MS
        t = pt.table(inputMS, readonly=True, ack=False)
        firstTime = t.getcell("TIME", 0)
        lastTime = t.getcell("TIME", t.nrows()-1)
        intTime = t.getcell("INTERVAL", 0)
        print 'Integration time:\t%f sec' % (intTime)
        nTimeslots = (lastTime - firstTime) / intTime
        if timeslots[1] == -1:
                timeslots[1] = nTimeslots
        else:
                timeslots[1] += 1
        print 'Number of timeslots:\t%d' % (nTimeslots)
        # open the antenna and spectral window subtables
        tant = pt.table(t.getkeyword('ANTENNA'), readonly=True, ack=False)
        tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
        numChannels = len(tsp.getcell('CHAN_FREQ',0))
        print 'Number of channels:\t%d' % (numChannels)
        print 'Reference frequency:\t%5.2f MHz' % (tsp.getcell('REF_FREQUENCY',0)/1.e6)

        # Station names
        antList = tant.getcol('NAME')
        if len(antToPlot)==1 and antToPlot[0]==-1:
                antToPlot = range(len(antList))
        print 'Station list (only starred stations will be plotted):'
        for i in range(len(antList)):
                star = ' '
                if i in antToPlot: star = '*'
                print '%s %2d\t%s' % (star, i, antList[i])

        # Bail if we're in query mode
        if queryMode:
                return

        # select by time from the beginning, and only use specified antennas
        tsel = t.query('TIME >= %f AND TIME <= %f AND ANTENNA1 IN %s AND ANTENNA2 IN %s' % (firstTime+timeslots[0]*intTime,firstTime+timeslots[1]*intTime,str(antToPlot),str(antToPlot)))

        # values to use for each polarization
        plotColors = [1,2,3,4]
        labXPositions = [0.35,0.45,0.55,0.65]
        labYPositions = [1.0,1.0,1.0,1.0]
        if convertStokes:
                polLabels = ['I','Q','U','V']
        else:
                polLabels = ['XX','XY','YX','YY']

        # define nicely written axis labels
        axisLabels = {'time': 'Time',
                      'ha': 'Hour angle',
                      'chan': 'Channel',
                      'freq': 'Frequency [MHz]',
                      'amp': 'Visibility amplitude',
                      'real': 'Real part of visibility',
                      'imag': 'Imaginary part of visibility',
                      'phase': 'Visibility phase [radians]'}

        # Now we loop through the baselines
        ppgplot.pgpage()
        for tpart in tsel.iter(["ANTENNA1","ANTENNA2"]):
                if not keepPlotting: return
                ant1 = tpart.getcell("ANTENNA1", 0)
                ant2 = tpart.getcell("ANTENNA2", 0)
                if ant1 not in antToPlot or ant2 not in antToPlot: continue
                if ant1 == ant2:
                        if not showAutocorr:
                                continue
                # Get the values to plot, strategy depends on axis type
                if xaxis == 'time' or xaxis == 'ha':
                        xaxisvals = getXAxisVals(tpart, xaxis, channels)
                        yaxisvals = getYAxisVals(tpart, yaxis, column, operation, showFlags, flagCol, channels, doUnwrap, convertStokes)
                else:
                        xaxisvals = getXAxisVals(tsp, xaxis, channels)
                        yaxisvals = getYAxisVals(tpart, yaxis, column, operation, showFlags, flagCol, channels, doUnwrap, convertStokes, xaxistype=1)
                if xaxisvals == None: # This baseline must be empty, go to next one
                        print 'No good data on baseline %s - %s' % (antList[ant1],antList[ant2])
                        continue
                    
                if debug:
                        print xaxisvals.shape
                        print yaxisvals.shape
                        for r in range(len(xaxisvals)):
                                print '%s'%yaxisvals[r]
                if len(xaxisvals) != len(yaxisvals): # something is wrong
                        print 'Error: X and Y axis types incompatible'
                        return

                # Plot the data, each polarization in a different color
                ppgplot.pgsci(1)
                if xmin == '':
                        minx = xaxisvals.min()
                else:
                        minx = float(xmin)
                if xmax == '':
                        maxx = xaxisvals.max()
                else:
                        maxx = float(xmax)
                if ymin == '':
                        miny = yaxisvals.min()
                        if numpy.ma.getmaskarray(yaxisvals.min()):
                                print 'All data flagged on baseline %s - %s' % (antList[ant1],antList[ant2])
                                continue
                else:
                        miny = float(ymin)
                if ymax == '':
                        maxy = yaxisvals.max()
                else:
                        maxy = float(ymax)
                if minx == maxx:
                        minx -= 1.0
                        maxx += 1.0
                else:
                        diffx = maxx - minx
                        minx -= 0.02*diffx
                        maxx += 0.02*diffx
                if miny == maxy:
                        miny -= 1.0
                        maxy += 1.0
                else:
                        diffy = maxy - miny
                        miny -= 0.02*diffy
                        maxy += 0.02*diffy
                #ppgplot.pgpage()
                ppgplot.pgswin(minx,maxx,miny,maxy)
                if xaxis == 'time' or xaxis == 'ha':
                        ppgplot.pgtbox('ZHOBCNST',0.0,0,'BCNST',0.0,0)
                else:
                        ppgplot.pgbox('BCNST',0.0,0,'BCNST',0.0,0)
                
                #ppgplot.pglab(axisLabels[xaxis], axisLabels[yaxis], '%s - %s'%(antList[ant1],antList[ant2]))
                #ppgplot.pgmtxt('T', 3.0, 0.5, 0.5, inputMSbasename)
                
                ppgplot.pglab(axisLabels[xaxis], axisLabels[yaxis], inputMSbasename + '(' + getDataDescription(column) + '): %s - %s'%(antList[ant1],antList[ant2]))
               
               
                if operation != 0:
                    # some operations is defined
                    if operation == 1:
                        label = 'XX-YY'
                    elif operation == 2:
                        label = 'XY.YX*'
                    else:
                        print 'Special operation not defined'
                        return
            
                    ppgplot.pgsci(plotColors[0])
                    tmpvals = yaxisvals
                    print 'Baseline',antList[ant1],'-',antList[ant2],': Plotting',len(tmpvals[~tmpvals.mask]),'points of ' + label
                    ppgplot.pgpt(xaxisvals[~tmpvals.mask], tmpvals[~tmpvals.mask], 1)
                            
                    addInfo(showStats, tmpvals[~tmpvals.mask], label, labXPositions[1], labYPositions[1])
                else:
                    for j in polarizations:
                        ppgplot.pgsci(plotColors[j])
                        tmpvals = yaxisvals[:,j]
                        if j == polarizations[0]:
                                print 'Baseline',antList[ant1],'-',antList[ant2],': Plotting',len(tmpvals[~tmpvals.mask]),'points per polarization'
                        ppgplot.pgpt(xaxisvals[~tmpvals.mask], tmpvals[~tmpvals.mask], 1)
                        
                        addInfo(showStats, tmpvals[~tmpvals.mask], polLabels[j], labXPositions[j], labYPositions[j])
                ppgplot.pgpage()

        # Close the PGPLOT device
        ppgplot.pgclos()

        if xaxis=='ha':
            print 'Removing derived columns...'
            pt.removeDerivedMSCal(inputMS)
                
# Add information, i.e. the label inte plot and the statistics if required
def addInfo(showStats, values, label, xPosition, yPosition):
        
        if showStats:
                # Compute mean and standard deviation  of yaxis values
                mean = values.mean()
                std = values.std()

                print "   " + label + ":" + " mean = %.5f" % mean + "\t std = %.5f" % std


        ppgplot.pgmtxt('T', yPosition, xPosition, 0.5, label)

def getData(table,column):
    
    splitColumn = column.split(',')
    
    if len(splitColumn) == 1:
         return table.getcol(column)
    elif len(splitColumn) == 3:
        columnAData = table.getcol(splitColumn[0])
        operation = splitColumn[1]
        columnBData = table.getcol(splitColumn[2])
        
        if operation == '-':
            return (columnAData-columnBData)
        elif operation == '+':
            return (columnAData+columnBData)
    
    # If we reach this point it means that the column is not correct
    print 'Column to plot: ' + column + ' is not correct!'
    return None
    
def getDataDescription(column):
    
    splitColumn = column.split(',')
    
    if len(splitColumn) == 1:
         return column[0]
    elif len(splitColumn) == 3:
        columnACapital = splitColumn[0][0]
        operation = splitColumn[1]
        columnBCapital = splitColumn[2][0]
        return (columnACapital + operation + columnBCapital)
    return None

# Get the x axis data out of the table (time, chan or freq)
def getXAxisVals(table, axisname, channels):
        
    if axisname == 'time':
        tmp = table.getcol('TIME')
        if tmp != None:
            return numpy.array(tmp-tmp.min(),dtype=numpy.float)
        return None
    if axisname == 'ha':
        tmp = table.getcol('HA')
        if tmp!= None:
            return numpy.array(tmp*180./numpy.pi*3600./15.,dtype=numpy.float)
        return None
    elif axisname == 'chan':
        tmp = table.getcol('CHAN_FREQ')[0]
        if tmp != None:
            return numpy.array(range(len(tmp))[channels[0]:channels[1]],dtype=numpy.float)+1
        return None
    elif axisname == 'freq':
        tmp = table.getcol('CHAN_FREQ')[0]
        if tmp != None:
            return numpy.array(tmp[channels[0]:channels[1]],dtype=numpy.float)/1.e6
        return None
    else:
        print 'Error: Requested axis not implemented'
        return None
    
def rangephase(phase):
    for i in range(len(phase)):
        element = phase[i]
        inrange = False
        while not inrange:
            if element > numpy.pi:
                element -= 2*numpy.pi
            elif element <= -numpy.pi:
                element += 2*numpy.pi
            else:
                inrange = True  
        phase[i] = element
    return phase


# Get the axis data out of the table(s)
def getYAxisVals(table, axisname, column, operation, showFlags, flagCol, channels, unwrap, stokes, xaxistype=0):
        
    if operation == 0:
        if axisname == 'amp':
            tmp = getData(table,column)
            flg = table.getcol(flagCol)
            if tmp != None:
                if showFlags:
                    tmp2 = numpy.ma.array(tmp, dtype=None, mask=False)
                else:
                    tmp2 = numpy.ma.array(tmp, dtype=None, mask=flg)
                tmp2[numpy.isnan(tmp2)]=numpy.ma.masked
                if stokes:
                    tmp2 = numpy.ma.transpose(numpy.ma.array([tmp2[:,:,0]+tmp2[:,:,3],tmp2[:,:,0]-tmp2[:,:,3],tmp2[:,:,1]+tmp2[:,:,2],numpy.complex(0.,-1.)*(tmp2[:,:,2]-tmp2[:,:,1])],dtype=None,mask=tmp2.mask),(1,2,0))
                if xaxistype == 0:
                    tmp3 = numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=1)
                else:
                    tmp3 = numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=0)
                return numpy.ma.sqrt((tmp3.real)**2+(tmp3.imag)**2)
            return None
        elif axisname == 'real':
            tmp = getData(table,column)
            flg = table.getcol(flagCol)
            if tmp != None:
                if showFlags:
                    tmp2 = numpy.ma.array(tmp, dtype=None, mask=False)
                else:
                    tmp2 = numpy.ma.array(tmp, dtype=None, mask=flg)
                tmp2[numpy.isnan(tmp2)]=numpy.ma.masked
                if stokes:
                    tmp2 = numpy.ma.transpose(numpy.ma.array([tmp2[:,:,0]+tmp2[:,:,3],tmp2[:,:,0]-tmp2[:,:,3],tmp2[:,:,1]+tmp2[:,:,2],numpy.complex(0.,-1.)*(tmp2[:,:,2]-tmp2[:,:,1])],dtype=None,mask=tmp2.mask),(1,2,0))
                if xaxistype == 0:
                    return numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=1).real
                else:
                    return numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=0).real
            return None
        elif axisname == 'imag':
            tmp = getData(table,column)
            flg = table.getcol(flagCol)
            if tmp != None:
                if showFlags:
                    tmp2 = numpy.ma.array(tmp, dtype=None, mask=False)
                else:
                    tmp2 = numpy.ma.array(tmp, dtype=None, mask=flg)
                tmp2[numpy.isnan(tmp2)]=numpy.ma.masked
                if stokes:
                    tmp2 = numpy.ma.transpose(numpy.ma.array([tmp2[:,:,0]+tmp2[:,:,3],tmp2[:,:,0]-tmp2[:,:,3],tmp2[:,:,1]+tmp2[:,:,2],numpy.complex(0.,-1.)*(tmp2[:,:,2]-tmp2[:,:,1])],dtype=None,mask=tmp2.mask),(1,2,0))
                if xaxistype == 0:
                    return numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=1).imag
                else:
                    return numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=0).imag
            return None
        elif axisname == 'phase':
            tmp = getData(table,column)
            flg = table.getcol(flagCol)
            if tmp != None:
                if showFlags:
                    tmp2 = numpy.ma.array(tmp, dtype=None, mask=False)
                else:
                    tmp2 = numpy.ma.array(tmp, dtype=None, mask=flg)
                tmp2[numpy.isnan(tmp2)]=numpy.ma.masked
                if stokes:
                    tmp2 = numpy.ma.transpose(numpy.ma.array([tmp2[:,:,0]+tmp2[:,:,3],tmp2[:,:,0]-tmp2[:,:,3],tmp2[:,:,1]+tmp2[:,:,2],numpy.complex(0.,-1.)*(tmp2[:,:,2]-tmp2[:,:,1])],dtype=None,mask=tmp2.mask),(1,2,0))
                if xaxistype == 0:
                    avgvals = numpy.ma.array(numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=1), mask=numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=1).mask)
                    tmp3 = numpy.ma.arctan2(avgvals.imag,avgvals.real)
                else:
                    avgvals = numpy.ma.array(numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=0), mask=numpy.ma.mean(tmp2[:,channels[0]:channels[1],:],axis=0).mask)
                    tmp3 = numpy.ma.arctan2(avgvals.imag,avgvals.real)
                if unwrap:
                    return numpy.ma.array(numpy.unwrap(tmp3), mask=tmp3.mask)
                else:
                    return tmp3
            return None                    
        else:
            print 'Error: Requested axis not implemented'
            return None
    elif operation == 1:

        # Special operation : XX-YY
        xxindex = 0
        yyindex = 3
            
        if axisname == 'real':
            
            real = getYAxisVals(table, 'real', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)    
            xxreal=numpy.ma.array(real[:, xxindex], mask=(real[:, xxindex]).mask)
            yyreal=numpy.ma.array(real[:, yyindex], mask=(real[:, yyindex]).mask)
            
            realop = xxreal-yyreal
            return numpy.ma.transpose(realop)
        elif axisname == 'imag':
            imag = getYAxisVals(table, 'imag', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)
            xximag=numpy.ma.array(imag[:, xxindex], mask=(imag[:, xxindex]).mask)
            yyimag=numpy.ma.array(imag[:, yyindex], mask=(imag[:, yyindex]).mask)
            
            imagop = xximag-yyimag
            return numpy.ma.transpose(imagop)
        elif axisname == 'amp':
            real = getYAxisVals(table, 'real', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)    
            xxreal=numpy.ma.array(real[:, xxindex], mask=(real[:, xxindex]).mask)
            yyreal=numpy.ma.array(real[:, yyindex], mask=(real[:, yyindex]).mask)
            
            imag = getYAxisVals(table, 'imag', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)
            xximag=numpy.ma.array(imag[:, xxindex], mask=(imag[:, xxindex]).mask)
            yyimag=numpy.ma.array(imag[:, yyindex], mask=(imag[:, yyindex]).mask)
            
            realop = xxreal-yyreal
            imagop = xximag-yyimag
            
            return numpy.ma.transpose(numpy.ma.sqrt((realop)**2+(imagop)**2))
        elif axisname == 'phase':
        
            real = getYAxisVals(table, 'real', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)    
            xxreal=numpy.ma.array(real[:, xxindex], mask=(real[:, xxindex]).mask)
            yyreal=numpy.ma.array(real[:, yyindex], mask=(real[:, yyindex]).mask)
            
            imag = getYAxisVals(table, 'imag', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)
            xximag=numpy.ma.array(imag[:, xxindex], mask=(imag[:, xxindex]).mask)
            yyimag=numpy.ma.array(imag[:, yyindex], mask=(imag[:, yyindex]).mask)
            
            realop = xxreal-yyreal
            imagop = xximag-yyimag
            
            phaseop = rangephase(numpy.ma.arctan2(imagop,realop))
            if unwrap:
                
                unwrappedarray = numpy.ma.transpose(phaseop)
                return numpy.ma.array(numpy.unwrap(unwrappedarray), mask=unwrappedarray.mask)
            else:
                return numpy.ma.transpose(phaseop)
        else:
            print 'Error: Requested axis not implemented'
            return None    
    elif operation == 2:

        # Special operation : XY . YX*
        
        xyindex = 1
        yxindex = 2
            
        if axisname == 'real':
            
            amp = getYAxisVals(table, 'amp', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)    
            phase = getYAxisVals(table, 'phase', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)
            xyamp=numpy.ma.array(amp[:, xyindex], mask=(amp[:, xyindex]).mask)
            yxamp=numpy.ma.array(amp[:, yxindex], mask=(amp[:, yxindex]).mask)
            xyphase=numpy.ma.array(phase[:, xyindex], mask=(phase[:, xyindex]).mask)
            yxphase=numpy.ma.array(phase[:, yxindex], mask=(phase[:, yxindex]).mask)
            
            ampop = xyamp*yxamp
            phaseop = xyphase - yxphase
            
            return numpy.ma.transpose(ampop * numpy.ma.cos(phaseop))
        elif axisname == 'imag':
            amp = getYAxisVals(table, 'amp', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)    
            phase = getYAxisVals(table, 'phase', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)
            xyamp=numpy.ma.array(amp[:, xyindex], mask=(amp[:, xyindex]).mask)
            yxamp=numpy.ma.array(amp[:, yxindex], mask=(amp[:, yxindex]).mask)
            xyphase=numpy.ma.array(phase[:, xyindex], mask=(phase[:, xyindex]).mask)
            yxphase=numpy.ma.array(phase[:, yxindex], mask=(phase[:, yxindex]).mask)
            
            ampop = xyamp*yxamp
            phaseop = xyphase - yxphase
            
            return numpy.ma.transpose(ampop * numpy.ma.sin(phaseop))
        elif axisname == 'amp':
            
            amp = getYAxisVals(table, 'amp', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)    
            xyamp=numpy.ma.array(amp[:, xyindex], mask=(amp[:, xyindex]).mask)
            yxamp=numpy.ma.array(amp[:, yxindex], mask=(amp[:, yxindex]).mask)
            
            ampop = xyamp*yxamp
            return numpy.ma.transpose(ampop)
        elif axisname == 'phase':
        
            phase = getYAxisVals(table, 'phase', column, 0, showFlags, flagCol, channels, unwrap, False, xaxistype)
            xyphase=numpy.ma.array(phase[:, xyindex], mask=(phase[:, xyindex]).mask)
            yxphase=numpy.ma.array(phase[:, yxindex], mask=(phase[:, yxindex]).mask)
            
            phaseop = rangephase(xyphase - yxphase)

            if unwrap:
                
                unwrappedarray = numpy.ma.transpose(phaseop)
                return numpy.ma.array(numpy.unwrap(unwrappedarray), mask=unwrappedarray.mask)
            else:
                return numpy.ma.transpose(phaseop)
        else:
            print 'Error: Requested axis not implemented'
            return None     
    else:
        print 'Error: Requested axis not implemented'
        return None      
def signal_handler(signal, frame):
        global keepPlotting
        try:
                ppgplot.pgend()
                keepPlotting = False
        except:
                sys.exit(0)

opt = optparse.OptionParser()
opt.add_option('-i','--inms',help='Input MS to plot [no default]',default='')
opt.add_option('-d','--device',help='PGPLOT device to use for the plotting [default /xwin], for options use the string "?"',default='/xwin')
opt.add_option('-x','--xaxis',help='X axis type (choose from time|ha|chan|freq) [default time]',default='time',type='choice',choices=['time','ha','chan','freq'])
opt.add_option('-y','--yaxis',help='Y axis type (choose from amp|phase|real|imag) [default amp]',default='amp',type='choice',choices=['amp','phase','real','imag'])
opt.add_option('-a','--axlimits',help='Axis limits (comma separated in order: xmin,xmax,ymin,ymax), leave any of them blank to use data min/max [default ",,,"]',default=',,,')
opt.add_option('-n','--nxy',help='Number of subplots in x,y [default 3,2]',default='3,2')
opt.add_option('-c','--column',help='Column to plot. The options are DATA|CORRECTED_DATA|MODEL_DATA but also other columns that the user may have created. It is also possible to specify a combination of two columns. For example the user may want the CORRECTED_DATA-MODEL_DATA. In such example the user should write: CORRECTED_DATA,-,MODEL_DATA. Currently + and - are supported. [default DATA]',default='DATA')
opt.add_option('-t','--timeslots',help='Timeslots to use (comma separated and zero-based: start,end[inclusive]) [default 0,-1 = all] Negative values work like python slicing, but please note that the second index here is inclusive. If plotting channel or frequency on x-axis, this parameter sets the time averaging interval.',default='0,-1')
opt.add_option('-s','--channels',help='Channels to use (comma separated and zero-based: start,end[inclusive]) [default 0,-1 = all] Negative values work like python slicing, but please note that the second index here is inclusive. If plotting time on x-axis, this parameter sets the channel averaging interval.',default='0,-1')
opt.add_option('-e','--antennas',help='Antennas to use (comma separated list, zero-based) [default -1=all] Use -q to see a list of available antennas. Only antennas in this list are plotted. To specify an inclusive range of antennas use .. format, e.g. -e 0..9 requests the first 10 antennas.',default='-1',type='string')
opt.add_option('-p','--polar',help='Polarizations to plot (it does not convert, so use integers as in the MS). [default is 0,1,2,3]',default='0,1,2,3')
opt.add_option('-w','--wrap',help='Unwrap phase? [default False]',default=False,action='store_true')
opt.add_option('-f','--flag',help='Show flagged data? [default False]',default=False,action='store_true')
opt.add_option('-g','--colflag',help='Column that contains flags [default FLAG]',default='FLAG')
opt.add_option('-k','--stokes',help='Convert to Stokes IQUV? [default False]',default=False,action='store_true')
opt.add_option('-u','--autocorr',help='Show autocorrelations? [default False]',default=False,action='store_true')
opt.add_option('-b','--debug',help='Run in debug mode? [default False]',default=False,action='store_true')
opt.add_option('-q','--query',help='Query mode (quits after reading dimensions, use for unfamiliar MSs) [default False]',default=False,action='store_true')
opt.add_option('--gui',help='Use GUI-based frontend to specify plotting options [default False]',default=False,action='store_true')
opt.add_option('-m','--statistics',help='Show statistics (mean and standard deviation)[default False]',default=False,action='store_true')
opt.add_option('-o','--operation',help='Plot an special operation over the polarizations. (choose from 0|1|2). 0 is for none operation (normal polarizations are plotted), 1 is XX-YY and 2 is XY.YX*. If some operation is specified the options polar and stokes are ignored. [default is 0, i.e. none operation]',default='0',type='choice',choices=['0','1','2'])
options, arguments = opt.parse_args()

signal.signal(signal.SIGINT, signal_handler)
class GuiFrontend:
        """
        This class defines the gui frontend to the plotting routine.
        Its main purpose is to provide a slightly easier way to specify 
        plotting options.
        """

        def __init__(self, master, options):
                self.myoptions = options
                self.frame = Frame(master, width=768, height=576)
                ###self.frame.grid_propagate(0)
                self.frame.pack()
                self.quitButton = Button(self.frame, text="Quit", fg="red", command=self.frame.quit)
                self.plotButton = Button(self.frame, text="Plot", fg="green", command=self.returnoptions)
                self.fileButton = Button(self.frame, text="Browse...", command=self.getmslocation)
                self.inputLabel = Label(self.frame, text="Input MS")
                self.inputEntry = Entry(self.frame)
                self.inputEntry.insert(0, self.myoptions.inms)
                self.deviceLabel = Label(self.frame, text="Plot device")
                self.deviceEntry = Entry(self.frame)
                self.deviceButton = Button(self.frame, text="Options?", command=self.listdevices)
                self.xaxisLabel = Label(self.frame, text="X-axis type")
                self.xvar = StringVar(self.frame)
                self.xvar.set(self.myoptions.xaxis)
                self.xaxisDropbox = OptionMenu(self.frame, self.xvar, "time", "ha", "chan", "freq")
                self.xaxisDropbox.pack()
                self.yaxisLabel = Label(self.frame, text="Y-axis type")
                self.yvar = StringVar(self.frame)
                self.yvar.set(self.myoptions.yaxis)
                self.yaxisDropbox = OptionMenu(self.frame, self.yvar, "amp", "phase","real","imag", "xxyydiff")
                self.yaxisDropbox.pack()
                self.xLimLabel = Label(self.frame, text="X-axis min,max")
                self.xLimEntry1 = Entry(self.frame)
                self.xLimEntry1.insert(0,self.myoptions.axlimits.split(',')[0])
                self.xLimEntry2 = Entry(self.frame)
                self.xLimEntry2.insert(0,self.myoptions.axlimits.split(',')[1])
                self.yLimLabel = Label(self.frame, text="Y-axis min,max")
                self.yLimEntry1 = Entry(self.frame)
                self.yLimEntry1.insert(0,self.myoptions.axlimits.split(',')[2])
                self.yLimEntry2 = Entry(self.frame)
                self.yLimEntry2.insert(0,self.myoptions.axlimits.split(',')[3])
                self.nxyLabel = Label(self.frame, text="Num. subplots x,y")
                self.nxyEntry = Entry(self.frame)
                self.nxyEntry.insert(0, self.myoptions.nxy)
                self.columnLabel = Label(self.frame, text="Column to plot")
                self.columnEntry = Entry(self.frame)
                self.columnEntry.insert(0, self.myoptions.column)
                self.timeLabel = Label(self.frame, text="Timeslots to plot/avg")
                self.timeEntry = Entry(self.frame)
                self.timeEntry.insert(0, self.myoptions.timeslots)
                self.chanLabel = Label(self.frame, text="Channels to plot/avg")
                self.chanEntry = Entry(self.frame)
                self.chanEntry.insert(0, self.myoptions.channels)
                self.antLabel = Label(self.frame, text="Antennas to plot")
                self.antListFrame = Frame(self.frame)
                self.antScrollbar = Scrollbar(self.antListFrame, orient=VERTICAL)
                self.antList = Listbox(self.antListFrame,selectmode=EXTENDED,yscrollcommand=self.antScrollbar.set)
                self.antScrollbar.config(command=self.antList.yview)
                self.antButton = Button(self.frame, text="Get list", command=self.getantlist)
                self.polLabel = Label(self.frame, text="Polarizations")
                self.polEntry = Entry(self.frame)
                self.polEntry.insert(0, self.myoptions.polar)
                self.opLabel = Label(self.frame, text="Operation")
                self.opEntry = Entry(self.frame)
                self.opEntry.insert(0, self.myoptions.operation)
                self.polvar = IntVar()
                self.polvar.set(self.myoptions.stokes)
                self.stokesCheck = Checkbutton(self.frame, text="Convert to IQUV?", variable=self.polvar)
                self.stokesCheck.var = self.polvar
                self.wrapvar = IntVar()
                self.wrapvar.set(self.myoptions.wrap)
                self.wrapCheck = Checkbutton(self.frame, text="Unwrap phases?", variable=self.wrapvar)
                self.wrapCheck.var = self.wrapvar
                self.flagLabel = Label(self.frame, text="Flag column")
                self.flagEntry = Entry(self.frame)
                self.flagEntry.insert(0,self.myoptions.colflag)
                self.flagvar = IntVar()
                self.flagvar.set(self.myoptions.flag)
                self.flagCheck = Checkbutton(self.frame, text="Show flagged data?", variable=self.flagvar)
                self.flagCheck.var = self.flagvar
                self.acorvar = IntVar()
                self.acorvar.set(self.myoptions.autocorr)
                self.acorCheck = Checkbutton(self.frame, text="Show autocorrelations?", variable=self.acorvar)
                self.acorCheck.var = self.acorvar
                self.statsvar = IntVar()
                self.statsvar.set(self.myoptions.statistics)
                self.statsCheck = Checkbutton(self.frame, text="Show statistics?", variable=self.statsvar)
                self.statsCheck.var = self.statsvar

                guirow = 0
                self.inputLabel.grid(row=guirow,column=0,sticky=W)
                self.inputEntry.grid(row=guirow,column=1,columnspan=2,sticky=W+E)
                self.fileButton.grid(row=guirow,column=3)
                guirow += 1
                self.deviceLabel.grid(row=guirow,column=0,sticky=W)
                self.deviceEntry.grid(row=guirow,column=1,columnspan=2,sticky=W+E)
                self.deviceEntry.insert(0, self.myoptions.device)
                self.deviceButton.grid(row=guirow,column=3)
                guirow += 1
                self.xaxisLabel.grid(row=guirow,column=0,sticky=W)
                self.xaxisDropbox.grid(row=guirow,column=1,sticky=W)
                guirow += 1
                self.yaxisLabel.grid(row=guirow,column=0,sticky=W)
                self.yaxisDropbox.grid(row=guirow,column=1,sticky=W)
                self.wrapCheck.grid(row=guirow,column=2,sticky=W)
                guirow += 1
                self.xLimLabel.grid(row=guirow,column=0,sticky=W)
                self.xLimEntry1.grid(row=guirow,column=1)
                self.xLimEntry2.grid(row=guirow,column=2)
                guirow += 1
                self.yLimLabel.grid(row=guirow,column=0,sticky=W)
                self.yLimEntry1.grid(row=guirow,column=1)
                self.yLimEntry2.grid(row=guirow,column=2)
                guirow += 1
                self.nxyLabel.grid(row=guirow,column=0,sticky=W)
                self.nxyEntry.grid(row=guirow,column=1)
                guirow += 1
                self.columnLabel.grid(row=guirow,column=0,sticky=W)
                self.columnEntry.grid(row=guirow,column=1)
                self.acorCheck.grid(row=guirow,column=2,sticky=W)
                guirow += 1
                self.timeLabel.grid(row=guirow,column=0,sticky=W)
                self.timeEntry.grid(row=guirow,column=1)
                self.statsCheck.grid(row=guirow,column=2,sticky=W)
                guirow += 1
                self.chanLabel.grid(row=guirow,column=0,sticky=W)
                self.chanEntry.grid(row=guirow,column=1)
                guirow += 1
                self.antLabel.grid(row=guirow,column=0,sticky=W)
                self.antList.pack(side=LEFT,fill=BOTH,expand=1)
                self.antScrollbar.pack(side=LEFT,fill=BOTH)
                self.antListFrame.grid(row=guirow,column=1,columnspan=2,sticky=W+E+N+S)
                self.antButton.grid(row=guirow,column=3)
                guirow += 1
                self.polLabel.grid(row=guirow,column=0,sticky=W)
                self.polEntry.grid(row=guirow,column=1)
                self.stokesCheck.grid(row=guirow,column=2,sticky=W)
                guirow += 1
                self.opLabel.grid(row=guirow,column=0,sticky=W)
                self.opEntry.grid(row=guirow,column=1)
                guirow += 1
                self.flagLabel.grid(row=guirow,column=0,sticky=W)
                self.flagEntry.grid(row=guirow,column=1)
                self.flagCheck.grid(row=guirow,column=2,sticky=W)
                guirow += 1
                self.plotButton.grid(row=guirow,column=0, pady=10)
                self.quitButton.grid(row=guirow,column=3, pady=10)

                if self.myoptions.inms != '':
                        self.getantlist()

        def returnoptions(self):
                self.myoptions.inms = self.inputEntry.get()
                self.myoptions.device = self.deviceEntry.get()
                self.myoptions.xaxis = self.xvar.get()
                self.myoptions.yaxis = self.yvar.get()
                self.myoptions.wrap = bool(self.wrapvar.get())
                self.myoptions.axlimits = self.xLimEntry1.get()+','+self.xLimEntry2.get()+','+self.yLimEntry1.get()+','+self.yLimEntry2.get()
                self.myoptions.nxy = self.nxyEntry.get()
                self.myoptions.column = self.columnEntry.get()
                self.myoptions.timeslots = self.timeEntry.get()
                self.myoptions.channels = self.chanEntry.get()
                tmplist = map(int,self.antList.curselection())
                if tmplist == []:
                        tmpstring = '-1,'
                else:
                        tmpstring = ''
                        for a in tmplist:
                                if a != '':
                                        tmpstring += '%d,'%a
                self.myoptions.antennas = tmpstring[:-1]
                self.myoptions.polar = self.polEntry.get()
                self.myoptions.operation = self.opEntry.get()
                self.myoptions.stokes = bool(self.polvar.get())
                self.myoptions.colflag = self.flagEntry.get()
                self.myoptions.flag = bool(self.flagvar.get())
                self.myoptions.autocorr = bool(self.acorvar.get())
                self.myoptions.statistics = bool(self.statsvar.get())
                print 'Plotting with GUI-specified options ...\n'
                if self.myoptions.debug:
                        print self.myoptions
                main(self.myoptions)
                print '\nPlotting routine completed.'
                print 'You can change options and plot again, or quit.\n'

        def getmslocation(self):
                dirname = tkFileDialog.askdirectory(initialdir='.',title='Please select the MS to plot')
                self.inputEntry.delete(0, END)
                self.inputEntry.insert(0, dirname)

        def listdevices(self):
                print 'Listing of available PGPLOT devices:\n'
                ppgplot.pgldev()
                print '\n'

        def getantlist(self):
                print 'Listing antennas in MS '+self.inputEntry.get()+'\n'
                ttmp = pt.table(self.inputEntry.get(),readonly=True,ack=False)
                tant = pt.table(ttmp.getkeyword('ANTENNA'),readonly=True,ack=False)
                antlist = tant.getcol('NAME')
                self.antList.delete(0,END)
                ninserted = 0
                for ant in antlist:
                        self.antList.insert(END, ant)
                        self.antList.selection_set(ninserted)
                        ninserted += 1

def loader(options):
        if (options.gui or options.colflag == 'ui') and not options.query:
                if options.colflag == 'ui':
                        print 'Warning: You used \'-gui\' instead of \'--gui\'.'
                        print 'Resetting the colflag parameter to FLAG and using the gui interface.'
                        print 'If you really have a flag column called ui, this is a problem!\n'
                        options.colflag = 'FLAG'
                root = Tk()
                root.title('uvplot.py GUI frontend')
                # initialize the gui with any user-selected options
                mygui = GuiFrontend(root, options)
                # start the gui
                print '--------------------------------------------------------'
                print 'WARNING!!!! The gui frontend is still under development.'
                print 'Use caution - double check your inputs!!'
                print '--------------------------------------------------------\n\n'
                root.mainloop()
                # update the options with what was specified in the gui
                options = mygui.myoptions
        else:
                main(options)

loader(options)


