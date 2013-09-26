#!/usr/bin/env python

# Transfer solutions from calibrator scans to the other scans

# * Fill solutions that have amp=1 and phase=0 with the nearest good solution using interpolation when possible
# * Handle any combination of correlations
# * Handle any number of calibrators/targets and find the correct match
# * Handle multichannel solutions
# * Handle phasors solutions

# send bug report to: fdg@mpa-garching.mpg.de

# TODO:
# Flag amplitude outliers
# Carefully check time of cal solutions (we use MS time that is not necessarily the same)

import os,sys,copy
import numpy as np
import lofar.parmdb as pdb
import matplotlib.pyplot as plt
import pyrap.tables as pt
import scipy.interpolate as interp
from caltar import *

#################################################################
# Look for the closet calibrators around a target and return them
# if only one calibrator found return it in both gmin_f and gmax_f
def search(target, calibs, corr, ant):
	# global values
	gmin = 0
	gmin_f = ''
	gmax = 1e30
	gmax_f = ''
	for calib in calibs:
		if calib.GetTime('min') < target.GetTime('min') and calib.GetTime('min') > gmin and \
				calib.GetFreq('min') == target.GetFreq('min') and calib.GetFreq('max') == target.GetFreq('max'):
			gmin = calib.GetTime('min')
			gmin_f = calib
			prop_in_freq = False
		if calib.GetTime('max') > target.GetTime('max') and calib.GetTime('max') < gmax and \
				calib.GetFreq('min') == target.GetFreq('min') and calib.GetFreq('max') == target.GetFreq('max'):
			gmax = calib.GetTime('max')
			gmax_f = calib
			prop_in_freq = False
		if calib.GetFreq('min') < target.GetFreq('min') and calib.GetFreq('min') > gmin and \
                       		calib.GetTime('min') == target.GetTime('min') and calib.GetTime('max') == target.GetTime('max'):
                        gmin = calib.GetFreq('min')
                        gmin_f = calib
			prop_in_freq = True
                if calib.GetFreq('max') > target.GetFreq('max') and calib.GetFreq('max') < gmax and \
                		calib.GetTime('min') == target.GetTime('min') and calib.GetTime('max') == target.GetTime('max'):
                        gmax = calib.GetFreq('max')
                        gmax_f = calib
			prop_in_freq = True

	# if no calibrators found print a message
	if gmin_f == '' and gmax_f == '': print "No calibrator found for ", target.name
	# if only one calibrator founf return it as min AND max
	elif gmin_f == '': gmin_f = gmax_f
	elif gmax_f == '': gmax_f = gmin_f
	return gmin_f, gmax_f, prop_in_freq

########################################
# Plot solution of a target and the 1 or 2 calibrator(s) involved
def plot_sol(target,cals,ant,prop_in_freq):
	fig = plt.figure(figsize=(8, 8))
	ax1 = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)
	if prop_in_freq:
	  for c in cals: # plot only for the first time? is the [0]
	    num = len(c.GetAmp(ant,'0:0')[0])
	    ax1.plot(np.linspace(c.GetFreq('min'),c.GetFreq('max'),num), c.GetAmp(ant,'0:0')[0], 'ro')
	    ax1.plot(np.linspace(c.GetFreq('min'),c.GetFreq('max'),num), c.GetAmp(ant,'1:1')[0], 'bo')
	    ax2.plot(np.linspace(c.GetFreq('min'),c.GetFreq('max'),num), c.GetPh(ant,'0:0')[0], 'ro')
	    ax2.plot(np.linspace(c.GetFreq('min'),c.GetFreq('max'),num), c.GetPh(ant,'1:1')[0], 'bo')
	  num = len(target.GetAmp(ant,'0:0')[0])
	  ax1.plot(np.linspace(target.GetFreq('min'),target.GetFreq('max'),num), target.GetAmp(ant,'0:0')[0], 'r*')
	  ax1.plot(np.linspace(target.GetFreq('min'),target.GetFreq('max'),num), target.GetAmp(ant,'1:1')[0], 'b*')
	  ax2.plot(np.linspace(target.GetFreq('min'),target.GetFreq('max'),num), target.GetPh(ant,'0:0')[0], 'r*')
	  ax2.plot(np.linspace(target.GetFreq('min'),target.GetFreq('max'),num), target.GetPh(ant,'1:1')[0], 'b*')
	else:
	  for c in cals:
	    num = len(c.GetAmp(ant,'0:0').transpose()[0])
	    ax1.plot(np.linspace(c.GetTime('min'),c.GetTime('max'),num), c.GetAmp(ant,'0:0').transpose()[0], 'ro')
	    ax1.plot(np.linspace(c.GetTime('min'),c.GetTime('max'),num), c.GetAmp(ant,'1:1').transpose()[0], 'bo')
	    ax2.plot(np.linspace(c.GetTime('min'),c.GetTime('max'),num), c.GetPh(ant,'0:0').transpose()[0], 'ro')
	    ax2.plot(np.linspace(c.GetTime('min'),c.GetTime('max'),num), c.GetPh(ant,'1:1').transpose()[0], 'bo')
	  num = len(target.GetAmp(ant,'0:0').transpose()[0])
	  ax1.plot(np.linspace(target.GetTime('min'),target.GetTime('max'),num), target.GetAmp(ant,'0:0').transpose()[0], 'r*')
	  ax1.plot(np.linspace(target.GetTime('min'),target.GetTime('max'),num), target.GetAmp(ant,'1:1').transpose()[0], 'b*')
	  ax2.plot(np.linspace(target.GetTime('min'),target.GetTime('max'),num), target.GetPh(ant,'0:0').transpose()[0], 'r*')
	  ax2.plot(np.linspace(target.GetTime('min'),target.GetTime('max'),num), target.GetPh(ant,'1:1').transpose()[0], 'b*')
	  ax2.set_xlabel('Time')
	
	ax1.set_ylabel('Amplitude')
	ax2.set_ylabel('Phase')
	fig.show()
	var = raw_input("Press Enter to continue...")
	fig.savefig(target.GetName()+'_'+ant+'.png')


########################################
# Return the last time slot af a 2-D array (first axis time (external array), second axis freq (internal array)) in a 2-D array
def getLastTime(arr):
	return arr[-1,:]
def getFirstTime(arr):
	return arr[0,:]
def getLastFreq(arr):
	return arr.transpose()[-1,:].transpose()
def getFirstFreq(arr):
	return arr.transpose()[0,:].transpose()

#########################################################################
# Propagate solution to "target" from "cal_low" and "cal_hi" in frequency
def propagate_f(target, cal_low, cal_hi, amp_mod, ph_mod, corr, ant):
#	print cal_low.name, '+', cal_hi.name, '->', target.name
	# Amplitude
	if amp_mod == '0': new_amp = [[1]]
	else:
		cal_amp_hi = getFirstFreq(cal_hi.GetAmp(ant,corr))
		cal_amp_low = getLastFreq(cal_low.GetAmp(ant,corr))
		freq_hi = cal_hi.GetFreq('min')
		freq_low = cal_low.GetFreq('max')
		# handle single calibrator (cal_low==cal_hi)
		if cal_low == cal_hi:
			# single calibrator AFTER the target (cal_amp_low and cal_amp_hi set at the First cal value)
			if cal_low.GetFreq('centre') > target.GetFreq('centre'):
				cal_amp_low = cal_amp_hi
				freq_low = target.GetFreq('min')-1
			# single calibrator BEFORE the target (cal_amp_low and cal_amp_hi set at the Last cal value)
			else:
				cal_amp_hi = cal_amp_low
				freq_hi = target.GetFreq('max')+1

		if amp_mod == 'near': interpkind='nearest'
		else: interpkind = 'linear'
		cal_amp_low = cal_amp_low.transpose()
		cal_amp_hi = cal_amp_hi.transpose()
		interpolation = interp.interp1d([freq_low, freq_hi], [cal_amp_low, cal_amp_hi], axis=0, kind=interpkind)
		new_amp = interpolation(np.arange(target.GetFreq('min'),target.GetFreq('max'),target.GetFreq('step'))).transpose()

	# Phase
	if ph_mod == '0': new_ph = [[0]]
	else:
		cal_ph_hi = getFirstFreq(cal_hi.GetPh(ant,corr))
		cal_ph_low = getLastFreq(cal_low.GetPh(ant,corr))
		freq_hi = cal_hi.GetFreq('min')
		freq_low = cal_low.GetFreq('max')
		if cal_low == cal_hi:
			if cal_low.GetFreq('centre') > target.GetFreq('centre'):
				cal_ph_low = cal_ph_hi
				freq_low = target.GetFreq('min')-1
			else:
				cal_ph_hi = cal_ph_low
				freq_hi = target.GetFreq('max')+1

	        if ph_mod == 'near': interpkind='nearest'
		else: interpkind = 'linear'
		cal_ph_low = cal_ph_low.transpose()
		cal_ph_hi = cal_ph_hi.transpose()
		interpolation = interp.interp1d([freq_low, freq_hi], [cal_ph_low, cal_ph_hi], axis=0, kind=interpkind)
		new_ph = interpolation(np.arange(target.GetFreq('min'),target.GetFreq('max'),target.GetFreq('step'))).transpose()

	target.SetCom(new_amp, new_ph, corr, ant, [target.GetFreq('min'),target.GetFreq('max'),cal_low.GetTime('min'),cal_low.GetTime('max')])

####################################################################
# Propagate solution to "target" from "cal_low" and "cal_hi" in time
def propagate_t(target, cal_low, cal_hi, amp_mod, ph_mod, corr, ant):
#	print cal_low.name, '+', cal_hi.name, '->', target.name
	# Amplitude
	if amp_mod == '0': new_amp = [[1]]
	else:
		# use the last value of cal_low and the first of cal_hi
		cal_amp_hi = getFirstTime(cal_hi.GetAmp(ant,corr))
		cal_amp_low = getLastTime(cal_low.GetAmp(ant,corr))
		time_hi = cal_hi.GetTime('min')
		time_low = cal_low.GetTime('max')
		# handle single calibrator (cal_low==cal_hi)
		if cal_low == cal_hi:
			# single calibrator AFTER the target (cal_amp_low and cal_amp_hi set at the First cal value)
			if cal_low.GetTime('centre') > target.GetTime('centre'):
				cal_amp_low = cal_amp_hi
				time_low = target.GetTime('min')-1
			# single calibrator BEFORE the target (cal_amp_low and cal_amp_hi set at the Last cal value)
			else:
				cal_amp_hi = cal_amp_low
				time_hi = target.GetTime('max')+1

		if amp_mod == 'near': interpkind='nearest'
		else: interpkind = 'linear'
		interpolation = interp.interp1d([time_low, time_hi], [cal_amp_low, cal_amp_hi], axis=0, kind=interpkind)
		new_amp = interpolation(np.arange(target.GetTime('min'),target.GetTime('max'),target.GetTime('step')))

	# Phase
	if ph_mod == '0': new_ph = [[0]]
	else:
		cal_ph_hi = getFirstTime(cal_hi.GetPh(ant,corr))
		cal_ph_low = getLastTime(cal_low.GetPh(ant,corr))
		time_hi = cal_hi.GetTime('min')
		time_low = cal_low.GetTime('max')
		if cal_low == cal_hi:
			if cal_low.GetTime('centre') > target.GetTime('centre'):
				cal_ph_low = cal_ph_hi
				time_low = target.GetTime('min')-1
			else:
				cal_ph_hi = cal_ph_low
				time_hi = target.GetTime('max')+1

	        if ph_mod == 'near': interpkind = 'nearest'
		else: interpkind = 'linear'
		interpolation = interp.interp1d([time_low, time_hi], [cal_ph_low, cal_ph_hi], axis=0, kind=interpkind)
		new_ph = interpolation(np.arange(target.GetTime('min'),target.GetTime('max'),target.GetTime('step')))

	target.SetCom(new_amp, new_ph, corr, ant, [cal_low.GetTime('min'),cal_low.GetTime('max'),target.GetFreq('min'),target.GetFreq('max')])

######################################
# Options
import optparse
opt = optparse.OptionParser(usage="%prog -c calSBs -t targetSBs", version="%prog 0.2")
opt.add_option('-c', '--calibs', help='calibrator SBs (comma separated)', type="string")
opt.add_option('-t', '--targets', help='target SBs (comma separated)', type="string")
opt.add_option('-p', '--phase', help='lin=lin-interp (default), near=nearest value, 0=do not correc', type="string", default='lin')
opt.add_option('-a', '--amp', help='lin=lin-interp (default), near=nearest value, 0=do not correct', type="string", default='lin')
opt.add_option('-r', '--refant', help='[NOT YET IMPLEMENTED] Reference antenna (example=RS406LBA; default='')', type="string", default='')
opt.add_option('-m', '--mean', help='Must be [tf][ap], average in time(t)/freq(f) all the values of a calibrator amp(a)/phases(p) (example=aft; default='')', default='')
opt.add_option('-k', '--checksol', help='If True, check for bad solutions and remove them (default=False)', action="store_true", default=False)
opt.add_option('-u', '--update', help='If True, update the average value of already existant target solution (default=False)', action="store_true", default=False)
opt.add_option('-d', '--plots', help='List of antenna names (comma separated) for which plots have to be generated (example=RS406LBA; default='')', type="string", default='')
(options, args) = opt.parse_args()

# Options check
if options.calibs == None or options.targets == None:
	sys.exit("ERROR: need at lest 1 target and 1 calibrator")
if ( options.phase != 'lin' and options.phase != 'near' and options.phase != '0') or ( options.amp != 'lin' and options.amp != 'near' and options.amp != '0' ):
	sys.exit("ERROR: phase/amp options must be lin|near|0")
calibs_f = options.calibs.split(',')
targets_f = options.targets.split(',')

amp_mod = options.amp
if amp_mod != '0' and amp_mod != 'near': amp_mod == 'lin'
print "Type of solution propagation for amplitudes: ", amp_mod
ph_mod = options.phase
if ph_mod != 'lin' and ph_mod != 'near': ph_mod == '0'
print "Type of solution propagation for phases: ", ph_mod
refant = options.refant
if refant != '': print "Reference antenna: ", refant
mean_mode = options.mean
do_checksol = options.checksol
do_update = options.update
ant_plots = options.plots.split(',')

#####################
# compatibility tests
# check if all the calibrators have the same kind of solutions:
calibs = [Calib(i) for i in calibs_f]
names=calibs[0].db.getNames()
for calib in calibs:
  if names != calib.db.getNames():
    sys.exit("ERROR: calibrators do not have the same parameters names")
# check also for target if in the UPDATE mode
if do_update == True:
  targets = [TargetU(i) for i in targets_f]
  for target in targets:
    if names != target.db.getNames():
      sys.exit("ERROR: targets do not have the same parameters names of calibrators")
else: targets = [Target(i) for i in targets_f]

# set calibrators parameters
for calib in calibs:
	calib.MeanMode(mean_mode)
	if do_checksol == True: calib.SetIsCheckSol()
	if refant != '': calib.SetRef(refant)

print "Antennas: "+str(calibs[0].GetAntlist())

# find target solutions
for target in targets:
  print "Working Target:", target.GetName()
  for corr in calibs[0].GetCorrlist():
    i = 0
    for ant in calibs[0].GetAntlist():
	i += 1
	sys.stdout.write('\r')
    	print "Corr:", corr,
	sys.stdout.write(" [%-20s] %d%%" % ('='*(i*20/len(calibs[0].GetAntlist())), i*100/len(calibs[0].GetAntlist())))
	sys.stdout.flush()
	cal_low, cal_hi, prop_in_freq = search(target, calibs, corr, ant)
	if prop_in_freq:
		if cal_low == '' and cal_hi == '': pass
		else: propagate_f(target, cal_low, cal_hi, amp_mod, ph_mod, corr, ant)
	else:
		if cal_low == '' and cal_hi == '': pass
		else: propagate_t(target, cal_low, cal_hi, amp_mod, ph_mod, corr, ant)
    print ""

# write parmdb file and exec parmdb
for target in targets:
	target.RunParmdb()

# plots
targetsR = [SBr(i) for i in targets_f]
for target in targetsR:
	for ant in calibs[0].GetAntlist():
		if ant in ant_plots:
			cal_low, cal_hi, prop_in_freq = search(target, calibs, corr, ant)
			plot_sol(target,[cal_low,cal_hi],ant,prop_in_freq)

# DEBUG
#targets = [Calib(i) for i in targets_f]
#print calibs[0].GetAmp('CS001LBA','0:0')
#print calibs[1].GetAmp('CS001LBA','0:0')
#print targets[0].GetAmp('CS001LBA','0:0')
#print calibs[0].GetPh('CS001LBA','0:0')
#print calibs[1].GetPh('CS001LBA','0:0')
#print targets[0].GetPh('CS001LBA','0:0')
#fig = plt.figure(figsize=(8, 8))
#ax = fig.add_subplot(111)
#ax.set_xlabel(r'Time')
#ax.set_ylabel(r'Amp')
#ax.plot(np.arange(calibs[0].GetTime('min'),calibs[0].GetTime('max')-1,calibs[0].GetTime('step')), calibs[0].GetAmp('CS001LBA','0:0'), linestyle='-', color='blue')
#ax.plot(np.arange(calibs[1].GetTime('min'),calibs[1].GetTime('max')-1,calibs[1].GetTime('step')), calibs[1].GetAmp('CS001LBA','0:0'), linestyle='-', color='blue')
#ax.plot(np.arange(targets[0].GetTime('min'),targets[0].GetTime('max'),targets[0].GetTime('step')), targets[0].GetAmp('CS001LBA','0:0'), linestyle='-', color='red')
#fig.savefig('debug.png')

print "Done."
