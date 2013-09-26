#!/usr/bin/env python

"""
uvcov.py

This script is intended to provide a way to quickly plot uv coverage
for LOFAR datasets. It uses pyrap and ppgplot (numpy version).

If you use postscript device output, you may want to use 'embiggen.csh'
(located at /home/heald/bin/embiggen.csh) which makes the pointsizes bigger.

If you want square plots, you have to modify environment variables:
PGPLOT_PS_WIDTH
PGPLOT_PS_HEIGHT
Values of 7800 (equivalent to a physical size of 7.8 inches) work fine for me.

Written by George Heald
v1.0 completed 3/6/2010

10 June 2010  v1.1  Add choice to plot in kilolambda
 4 Aug  2010  v1.2  Add choice to assume same u,v in meters
10 Feb  2011  v1.3  Add option to specify title of plot
                    Add ability to plot broadband uv coverage from one MS
                    Allow antenna ranges instead of only lists (in -e)
 7 July 2011  v1.4  Add multicolor plot for different LOFAR baseline types

To do:
- Fix antenna selection for plotting >1 MS

"""

import optparse
import glob
import signal
import sys
import ppgplot
# Note, for documentation on ppgplot, see
# http://www.astro.rug.nl/~breddels/python/ppgplot/ppgplot-doc.txt
import numpy
import pyrap.tables as pt

version_string = 'v1.4, 7 July 2011'
print 'uvcov.py',version_string
print ''

def main(options):

	debug = options.debug
        MSlist = []
	device = options.device
	if device=='?':
		ppgplot.pgldev()
		return
        for inmspart in options.inms.split(','):
                for msname in glob.iglob(inmspart):
	                MSlist.append(msname)
	if len(MSlist) == 0:
		print 'Error: You must specify at least one MS name.'
		print '       Use "uvplot.py -h" to get help.'
		return
        if len(MSlist) > 1:
                print 'WARNING: Antenna selection (other than all) may not work well'
                print '         when plotting more than one MS. Carefully inspect the'
                print '         listings of antenna numbers/names!'
        if options.title == '':
                plottitle = options.inms
        else:
                plottitle = options.title
	axlimits = options.axlimits.split(',')
	if len(axlimits) == 4:
		xmin,xmax,ymin,ymax = axlimits
	else:
		print 'Error: You must specify four axis limits'
		return
	timeslots = options.timeslots.split(',')
	if len(timeslots) != 3:
		print 'Error: Timeslots format is start,skip,end'
		return
	for i in range(len(timeslots)):
		timeslots[i] = int(timeslots[i])
		if timeslots[i] < 0:
			print 'Error: timeslots values must not be negative'
			return
        doPlotColors = options.colors
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
	queryMode = options.query
        plotLambda = options.kilolambda

        badval = 0.0
        xaxisvals0 = numpy.array([])
        yaxisvals0 = numpy.array([])
        xaxisvals1 = numpy.array([])
        yaxisvals1 = numpy.array([])
        xaxisvals2 = numpy.array([])
        yaxisvals2 = numpy.array([])
        xaxisvals3 = numpy.array([])
        yaxisvals3 = numpy.array([])
        xaxisvals4 = numpy.array([])
        yaxisvals4 = numpy.array([])
        xaxisvals5 = numpy.array([])
        yaxisvals5 = numpy.array([])
        savex0 = numpy.array([])
        savey0 = numpy.array([])
        savex1 = numpy.array([])
        savey1 = numpy.array([])
        savex2 = numpy.array([])
        savey2 = numpy.array([])
        savex3 = numpy.array([])
        savey3 = numpy.array([])
        savex4 = numpy.array([])
        savey4 = numpy.array([])
        savex5 = numpy.array([])
        savey5 = numpy.array([])
        numPlotted = 0
        ptcolor = 0
        for inputMS in MSlist:
	        # open the main table and print some info about the MS
                print 'Getting info for', inputMS
	        t = pt.table(inputMS, readonly=True, ack=False)
                tfreq = pt.table(t.getkeyword('SPECTRAL_WINDOW'),readonly=True,ack=False)
                ref_freq = tfreq.getcol('REF_FREQUENCY',nrow=1)[0]
                ch_freq = tfreq.getcol('CHAN_FREQ',nrow=1)[0]
                print 'Reference frequency:\t%f MHz' % (ref_freq/1.e6)
                if options.wideband:
                        ref_wavelength = 2.99792458e8/ch_freq
                else:
                        ref_wavelength = [2.99792458e8/ref_freq]
                print 'Reference wavelength:\t%f m' % (ref_wavelength[0])
                if options.sameuv and numPlotted > 0:
                        print 'Assuming same uvw as first MS!'
                        if plotLambda:
                                for w in ref_wavelength:
                                        xaxisvals0 = numpy.append(xaxisvals0,[savex0/w/1000.,-savex0/w/1000.])
                                        yaxisvals0 = numpy.append(yaxisvals0,[savey0/w/1000.,-savey0/w/1000.])
                                        xaxisvals1 = numpy.append(xaxisvals1,[savex1/w/1000.,-savex1/w/1000.])
                                        yaxisvals1 = numpy.append(yaxisvals1,[savey1/w/1000.,-savey1/w/1000.])
                                        xaxisvals2 = numpy.append(xaxisvals2,[savex2/w/1000.,-savex2/w/1000.])
                                        yaxisvals2 = numpy.append(yaxisvals2,[savey2/w/1000.,-savey2/w/1000.])
                                        xaxisvals3 = numpy.append(xaxisvals3,[savex3/w/1000.,-savex3/w/1000.])
                                        yaxisvals3 = numpy.append(yaxisvals3,[savey3/w/1000.,-savey3/w/1000.])
                                        xaxisvals4 = numpy.append(xaxisvals4,[savex4/w/1000.,-savex4/w/1000.])
                                        yaxisvals4 = numpy.append(yaxisvals4,[savey4/w/1000.,-savey4/w/1000.])
                                        xaxisvals5 = numpy.append(xaxisvals5,[savex5/w/1000.,-savex5/w/1000.])
                                        yaxisvals5 = numpy.append(yaxisvals5,[savey5/w/1000.,-savey5/w/1000.])
                        else:
                                print 'Plotting more than one MS with same uv, all in meters... do you want -k?'
                                xaxisvals0 = numpy.append(xaxisvals0,[savex0,-savex0])
                                yaxisvals0 = numpy.append(yaxisvals0,[savey0,-savey0])
                                xaxisvals1 = numpy.append(xaxisvals1,[savex1,-savex1])
                                yaxisvals1 = numpy.append(yaxisvals1,[savey1,-savey1])
                                xaxisvals2 = numpy.append(xaxisvals2,[savex2,-savex2])
                                yaxisvals2 = numpy.append(yaxisvals2,[savey2,-savey2])
                                xaxisvals3 = numpy.append(xaxisvals3,[savex3,-savex3])
                                yaxisvals3 = numpy.append(yaxisvals3,[savey3,-savey3])
                                xaxisvals4 = numpy.append(xaxisvals4,[savex4,-savex4])
                                yaxisvals4 = numpy.append(yaxisvals4,[savey4,-savey4])
                                xaxisvals5 = numpy.append(xaxisvals5,[savex5,-savex5])
                                yaxisvals5 = numpy.append(yaxisvals5,[savey5,-savey5])
                        continue
                        
	        firstTime = t.getcell("TIME", 0)
	        lastTime = t.getcell("TIME", t.nrows()-1)
	        intTime = t.getcell("INTERVAL", 0)
	        print 'Integration time:\t%f sec' % (intTime)
	        nTimeslots = (lastTime - firstTime) / intTime
	        print 'Number of timeslots:\t%d' % (nTimeslots)
                if timeslots[1] == 0:
                        if nTimeslots >= 100:
                                timeskip = int(nTimeslots/100)
                        else:
                                timeskip = 1
                else:
                        timeskip = int(timeslots[1])
                print 'For each baseline, plotting one point every %d samples' % (timeskip)
       	        if timeslots[2] == 0:
        		timeslots[2] = nTimeslots
        	# open the antenna subtable
        	tant = pt.table(t.getkeyword('ANTENNA'), readonly=True, ack=False)
        
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
        	tsel = t.query('TIME >= %f AND TIME <= %f AND ANTENNA1 IN %s AND ANTENNA2 IN %s' % (firstTime+timeslots[0]*intTime,firstTime+timeslots[2]*intTime,str(antToPlot),str(antToPlot)), columns='ANTENNA1,ANTENNA2,UVW')

        	# Now we loop through the baselines
                i = 0
                nb = (len(antToPlot)*(len(antToPlot)-1))/2
                sys.stdout.write('Reading uvw for %d baselines: %04d/%04d'%(nb,i,nb))
                sys.stdout.flush()
	        for tpart in tsel.iter(["ANTENNA1","ANTENNA2"]):
        		ant1 = tpart.getcell("ANTENNA1", 0)
        		ant2 = tpart.getcell("ANTENNA2", 0)
                        if ant1 not in antToPlot or ant2 not in antToPlot: continue
        		if ant1 == ant2: continue
                        i += 1
                        sys.stdout.write('\b\b\b\b\b\b\b\b\b%04d/%04d'%(i,nb))
                        sys.stdout.flush()
                        if doPlotColors:
                                stNameStr = antList[ant1][0]+antList[ant2][0]
                                if   stNameStr == 'CC': ptcolor = 0
                                elif stNameStr == 'RR': ptcolor = 1
                                elif 'C' in stNameStr and 'R' in stNameStr: ptcolor = 2
                                elif 'C' in stNameStr: ptcolor = 3
                                elif 'R' in stNameStr: ptcolor = 4
                                else: ptcolor = 5
        		# Get the values to plot
                        uvw = tpart.getcol('UVW', rowincr=timeskip)
                        if numPlotted == 0:
                                savex0 = numpy.append(savex0,[uvw[:,0],-uvw[:,0]])
                                savey0 = numpy.append(savey0,[uvw[:,1],-uvw[:,1]])
                                savex1 = numpy.append(savex1,[uvw[:,0],-uvw[:,0]])
                                savey1 = numpy.append(savey1,[uvw[:,1],-uvw[:,1]])
                                savex2 = numpy.append(savex2,[uvw[:,0],-uvw[:,0]])
                                savey2 = numpy.append(savey2,[uvw[:,1],-uvw[:,1]])
                                savex3 = numpy.append(savex3,[uvw[:,0],-uvw[:,0]])
                                savey3 = numpy.append(savey3,[uvw[:,1],-uvw[:,1]])
                                savex4 = numpy.append(savex4,[uvw[:,0],-uvw[:,0]])
                                savey4 = numpy.append(savey4,[uvw[:,1],-uvw[:,1]])
                                savex5 = numpy.append(savex5,[uvw[:,0],-uvw[:,0]])
                                savey5 = numpy.append(savey5,[uvw[:,1],-uvw[:,1]])
                        if plotLambda:
                                for w in ref_wavelength:
                                        if ptcolor == 0:
                                                xaxisvals0 = numpy.append(xaxisvals0,[uvw[:,0]/w/1000.,-uvw[:,0]/w/1000.])
                                                yaxisvals0 = numpy.append(yaxisvals0,[uvw[:,1]/w/1000.,-uvw[:,1]/w/1000.])
                                        elif ptcolor == 1:
                                                xaxisvals1 = numpy.append(xaxisvals1,[uvw[:,0]/w/1000.,-uvw[:,0]/w/1000.])
                                                yaxisvals1 = numpy.append(yaxisvals1,[uvw[:,1]/w/1000.,-uvw[:,1]/w/1000.])
                                        elif ptcolor == 2:
                                                xaxisvals2 = numpy.append(xaxisvals2,[uvw[:,0]/w/1000.,-uvw[:,0]/w/1000.])
                                                yaxisvals2 = numpy.append(yaxisvals2,[uvw[:,1]/w/1000.,-uvw[:,1]/w/1000.])
                                        elif ptcolor == 3:
                                                xaxisvals3 = numpy.append(xaxisvals3,[uvw[:,0]/w/1000.,-uvw[:,0]/w/1000.])
                                                yaxisvals3 = numpy.append(yaxisvals3,[uvw[:,1]/w/1000.,-uvw[:,1]/w/1000.])
                                        elif ptcolor == 4:
                                                xaxisvals4 = numpy.append(xaxisvals4,[uvw[:,0]/w/1000.,-uvw[:,0]/w/1000.])
                                                yaxisvals4 = numpy.append(yaxisvals4,[uvw[:,1]/w/1000.,-uvw[:,1]/w/1000.])
                                        elif ptcolor == 5:
                                                xaxisvals5 = numpy.append(xaxisvals5,[uvw[:,0]/w/1000.,-uvw[:,0]/w/1000.])
                                                yaxisvals5 = numpy.append(yaxisvals5,[uvw[:,1]/w/1000.,-uvw[:,1]/w/1000.])
                        else:
                                if ptcolor == 0:
                                        xaxisvals0 = numpy.append(xaxisvals0,[uvw[:,0],-uvw[:,0]])
                                        yaxisvals0 = numpy.append(yaxisvals0,[uvw[:,1],-uvw[:,1]])
                                elif ptcolor == 1:
                                        xaxisvals1 = numpy.append(xaxisvals1,[uvw[:,0],-uvw[:,0]])
                                        yaxisvals1 = numpy.append(yaxisvals1,[uvw[:,1],-uvw[:,1]])
                                elif ptcolor == 2:
                                        xaxisvals2 = numpy.append(xaxisvals2,[uvw[:,0],-uvw[:,0]])
                                        yaxisvals2 = numpy.append(yaxisvals2,[uvw[:,1],-uvw[:,1]])
                                elif ptcolor == 3:
                                        xaxisvals3 = numpy.append(xaxisvals3,[uvw[:,0],-uvw[:,0]])
                                        yaxisvals3 = numpy.append(yaxisvals3,[uvw[:,1],-uvw[:,1]])
                                elif ptcolor == 4:
                                        xaxisvals4 = numpy.append(xaxisvals4,[uvw[:,0],-uvw[:,0]])
                                        yaxisvals4 = numpy.append(yaxisvals4,[uvw[:,1],-uvw[:,1]])
                                elif ptcolor == 5:
                                        xaxisvals5 = numpy.append(xaxisvals5,[uvw[:,0],-uvw[:,0]])
                                        yaxisvals5 = numpy.append(yaxisvals5,[uvw[:,1],-uvw[:,1]])
        		#if debug:
                        #        print uvw.shape
        		#	print xaxisvals.shape
        		#	print yaxisvals.shape
                        #else:
                        #        sys.stdout.write('.')
                        #        sys.stdout.flush()
                sys.stdout.write(' Done!\n')
                numPlotted += 1

        print 'Plotting uv points ...'
	# open the graphics device, using only one panel
	ppgplot.pgbeg(device, 1, 1)
	# set the font size
	ppgplot.pgsch(1)
	ppgplot.pgvstd()

        xaxisvals = numpy.append(xaxisvals0,numpy.append(xaxisvals1,numpy.append(xaxisvals2,numpy.append(xaxisvals3,numpy.append(xaxisvals4,xaxisvals5)))))
        yaxisvals = numpy.append(yaxisvals0,numpy.append(yaxisvals1,numpy.append(yaxisvals2,numpy.append(yaxisvals3,numpy.append(yaxisvals4,yaxisvals5)))))
        tmpvals0 = numpy.sqrt(xaxisvals0**2+yaxisvals0**2)
        tmpvals1 = numpy.sqrt(xaxisvals1**2+yaxisvals1**2)
        tmpvals2 = numpy.sqrt(xaxisvals2**2+yaxisvals2**2)
        tmpvals3 = numpy.sqrt(xaxisvals3**2+yaxisvals3**2)
        tmpvals4 = numpy.sqrt(xaxisvals4**2+yaxisvals4**2)
        tmpvals5 = numpy.sqrt(xaxisvals5**2+yaxisvals5**2)
	# Plot the data
        if debug:
                print xaxisvals0[tmpvals0!=badval]
                print yaxisvals0[tmpvals0!=badval]
	ppgplot.pgsci(1)
        uvmax = max(xaxisvals.max(),yaxisvals.max())
        uvmin = min(xaxisvals.min(),yaxisvals.min())
        uvuplim = 0.02*(uvmax-uvmin)+uvmax
        uvlolim = uvmin-0.02*(uvmax-uvmin)
	if xmin == '':
		minx = uvlolim
	else:
		minx = float(xmin)
	if xmax == '':
		maxx = uvuplim
	else:
		maxx = float(xmax)
	if ymin == '':
		miny = uvlolim
	else:
		miny = float(ymin)
	if ymax == '':
		maxy = uvuplim
	else:
		maxy = float(ymax)
	if minx == maxx:
		minx = -1.0
		maxx = 1.0
	if miny == maxy:
		miny = -1.0
		maxy = 1.0
        ppgplot.pgpage()
	ppgplot.pgswin(minx,maxx,miny,maxy)
        ppgplot.pgbox('BCNST',0.0,0,'BCNST',0.0,0)
        if plotLambda:
	        ppgplot.pglab('u [k\gl]', 'v [k\gl]', '%s'%(plottitle))
        else:
	        ppgplot.pglab('u [m]', 'v [m]', '%s'%(plottitle))
	ppgplot.pgpt(xaxisvals0[tmpvals0!=badval], yaxisvals0[tmpvals0!=badval], 1)
        #if doPlotColors: ppgplot.pgmtxt('T', 1, 0.35, 0.5, 'C-C')
        ppgplot.pgsci(2)
	ppgplot.pgpt(xaxisvals1[tmpvals1!=badval], yaxisvals1[tmpvals1!=badval], 1)
        #if doPlotColors: ppgplot.pgmtxt('T', 1, 0.50, 0.5, 'R-R')
        ppgplot.pgsci(4)
	ppgplot.pgpt(xaxisvals2[tmpvals2!=badval], yaxisvals2[tmpvals2!=badval], 1)
        #if doPlotColors: ppgplot.pgmtxt('T', 1, 0.65, 0.5, 'C-R')
        ppgplot.pgsci(3)
	ppgplot.pgpt(xaxisvals3[tmpvals3!=badval], yaxisvals3[tmpvals3!=badval], 1)
        #if doPlotColors: ppgplot.pgmtxt('T', 1, 0.55, 0.5, 'C-I')
        ppgplot.pgsci(5)
	ppgplot.pgpt(xaxisvals4[tmpvals4!=badval], yaxisvals4[tmpvals4!=badval], 1)
        #if doPlotColors: ppgplot.pgmtxt('T', 1, 0.65, 0.5, 'R-I')
        ppgplot.pgsci(6)
	ppgplot.pgpt(xaxisvals5[tmpvals5!=badval], yaxisvals5[tmpvals5!=badval], 1)
        #if doPlotColors: ppgplot.pgmtxt('T', 1, 0.75, 0.5, 'I-I')

	# Close the PGPLOT device
	ppgplot.pgclos()
		

def signal_handler(signal, frame):
        sys.exit(0)

opt = optparse.OptionParser()
opt.add_option('-i','--inms',help='Input MS(s) to plot [no default]. Multiple MS names can be given together, separated by commas. Wildcards are also accepted in order to make it easier to plot more than one MS at a time.',default='')
opt.add_option('-d','--device',help='PGPLOT device to use for the plotting [default /xs], for options use the string "?"',default='/xs')
opt.add_option('-a','--axlimits',help='Axis limits (comma separated in order: xmin,xmax,ymin,ymax), leave any of them blank to use data min/max [default ",,,"]',default=',,,')
opt.add_option('-t','--timeslots',help='Timeslots to use (comma separated and zero-based: start,skip,end) [default 0,0,0 = full time range, skipping such that 100 points are plotted per baseline]',default='0,0,0')
opt.add_option('-e','--antennas',help='Antennas to use (comma separated list, zero-based) [default -1=all] Use -q to see a list of available antennas. Only antennas in this list are plotted. When plotting more than one MS this may not work well. To specify an inclusive range of antennas use .. format, e.g. -e 0..9 requests the first 10 antennas.',default='-1',type='string')
opt.add_option('-k','--kilolambda',help='Plot in kilolambda rather than meters? [default False]',default=False,action='store_true')
opt.add_option('-w','--wideband',help='Plot each channel separately? Only useful with -k. [default False]',default=False,action='store_true')
opt.add_option('-s','--sameuv',help='Assume same uv coordinates (in meters) for multiple MSs? This is useful if all input MSs are SBs of a single observation. It is NOT useful when combining MSs from different timeranges. [default False]',default=False,action='store_true')
opt.add_option('--title',help="Plot title [default ''=MS name] If you need no title at all, use --title=' '",default='')
opt.add_option('-b','--debug',help='Run in debug mode? [default False]',default=False,action='store_true')
opt.add_option('-q','--query',help='Query mode (quits after reading dimensions, use for unfamiliar MSs) [default False]',default=False,action='store_true')
opt.add_option('-c','--colors',help='Plot different LOFAR baseline classes with different colors? [default False]',default=False,action='store_true')
options, arguments = opt.parse_args()

signal.signal(signal.SIGINT, signal_handler)

main(options)

