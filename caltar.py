
import os,sys,copy
from subprocess import *
import numpy as np
import lofar.parmdb as pdb
import pyrap.tables as pt

#################################
# general class to handle MS
class MS:
	def __init__(self, f):
		self.time = {'min':0,'max':0,'step':0,'centre':0}
		self.freq = {'min':0,'max':0,'step':0,'centre':0}
		self.name = f
		if not os.path.isdir(f): sys.exit('ERROR: wrong MS name '+f)
		# find freq
		t = pt.table(f+'/SPECTRAL_WINDOW',ack=False,readonly=True)
		freq_central = t.getcell('REF_FREQUENCY', 0)
		bw = t.getcell('TOTAL_BANDWIDTH', 0)
		self.freq['max'] = freq_central+bw/2 #TODO: check if this is correct
		self.freq['min'] = freq_central-bw/2
		self.freq['step'] = t.getcell('CHAN_WIDTH', 0)[0]
		self.freq['centre'] = freq_central
		self.freq['n'] = t.getcol('CHAN_FREQ')[0].shape[0]
		t.close()
		# find time
		t = pt.table(f,ack=False,readonly=True)
		self.time['step'] = t.getcell('EXPOSURE',0)
		self.time['min'] = min(t.getcol('TIME'))-self.time['step']/2 #TODO: check if this is correct
		self.time['max'] = max(t.getcol('TIME'))+self.time['step']/2
		self.time['centre'] = self.time['min'] + (self.time['min'] - self.time['max'])/2
		self.time['n'] = t.getcol('TIME').shape[0]
		t.close()
		self.amp = {}
		self.ph = {}

	def GetTime(self, p=''):
		if p != '': return self.time[p]
		else: return self.time
	def GetFreq(self, p=''):
		if p != '': return self.freq[p]
		else: return self.freq
	def GetName(self): return self.name

####################################
# class that defines a parmdb (reading)
class ParmdbR():
	def __init__(self, f):
		self.db = pdb.parmdb(f)
		gtypescorrs,parms,ants = np.array([ i.rsplit(':',2) for i in self.db.getNames() ]).transpose()
		gtypes, corrs = np.array([ i.split(':',1) for i in gtypescorrs ]).transpose()
		self.gtype = list(set(gtypes))
		if len(self.gtype) > 1: sys.exit('ERROR: cannot handle more than one solution type')
		else: self.gtype = self.gtype[0]
		self.corrlist = list(set(corrs))
		self.parmlist = list(set(parms))
		self.antlist = list(set(ants))
	
	#  Set wheter perform a check on solutions
	def SetIsCheckSol(self):
		self.ischecksol = True
	def IsCheckSol(self):
		if hasattr(self,'ischecksol'): return True
		else: return False

	# Set whether do a freq/time average of the solutions
	def MeanMode(self, mean_mode):
		self.mean_mode = mean_mode
	def GetMeanMode(self):
		if hasattr(self,'mean_mode'): return self.mean_mode
		else: return ''

	def GetAntlist(self):
		return self.antlist
	def GetCorrlist(self):
		return self.corrlist
	def GetParmlist(self):
		return self.parmlist
	def GetGType(self):
		return self.gtype

	# Return amp and ph of a given antenna/corr combination, return None if it does not exist
	def GetAmp(self, ant, corr):
		if not (ant in self.GetAntlist() and corr in self.GetCorrlist()): sys.exit('ERROR: bad Amp request, Ant:', ant, 'Corr:', corr)
		else: amp, ph = self._GetSol(ant, corr)
		return amp
	def GetPh(self, ant, corr):
		if not (ant in self.GetAntlist() and corr in self.GetCorrlist()): sys.exit('ERROR: bad Amp request, Ant:', ant, 'Corr:', corr)
		else: amp, ph = self._GetSol(ant, corr)
		return ph

	# Extract the solutions and put them in two arrays amp and ph
	def _GetSol(self, ant, corr):
		if 'Real' in self.GetParmlist() and 'Imag' in self.GetParmlist():
			re = self.db.getValuesGrid(self.GetGType()+':'+corr+':Real:'+ant)[self.GetGType()+':'+corr+':Real:'+ant]['values']
			im = self.db.getValuesGrid(self.GetGType()+':'+corr+':Imag:'+ant)[self.GetGType()+':'+corr+':Imag:'+ant]['values']
		elif 'Ampl' in self.GetParmlist() and 'Phase' in self.GetParmlist():
			amp = self.db.getValuesGrid(self.GetGType()+':'+corr+':Ampl:'+ant)[self.GetGType()+':'+corr+':Ampl:'+ant]['values']
			ph = self.db.getValuesGrid(self.GetGType()+':'+corr+':Phase:'+ant)[self.GetGType()+':'++corr+':Phase:'+ant]['values']
			re = amp*np.cos(ph)
			im = amp*np.sin(ph)
		else: sys.exit('ERROR: invalid parameter names:', self.GetParmlist())

		# Check first in time (fix freq) and then in freq (fix time)
		if self.IsCheckSol() == True:
			re = re.transpose()
			im = im.transpose()
			for f in xrange(len(re)):
				re[f], im[f] = self.__CheckSol(re[f],im[f])
			re = re.transpose()
			im = im.transpose()
			for t in xrange(len(re)):
				re[t], im[t] = self.__CheckSol(re[t],im[t])
		
		amp = np.sqrt((re**2)+(im**2))
		ph = np.arctan2(im, re)
		# Perform flagging
		#if IsFlagSol() == True: amp = self.__FlagSol(amp) # TODO
		# Eventually, apply a mean
		mean_mode = self.GetMeanMode()
		if 'a' in mean_mode and 't' in mean_mode: amp = np.array([np.mean(amp,0)])
		if 'a' in mean_mode and 'f' in mean_mode: amp = np.array([[i] for i in np.mean(amp,1)])
		if 'p' in mean_mode and 't' in mean_mode: ph = np.array([np.mean(ph,0)])
		if 'p' in mean_mode and 'f' in mean_mode: ph = np.array([[i] for i in np.mean(ph,1)])
		return amp, ph

	# Check 2 1-D arrays (re,im) and fill bad values with the closest value possible (near on the edges or interp)
	def __CheckSol(self, re, im):
		# Only bad solutions, cannot do anything
		if (np.all(re == 1) and np.all(im == 0)) or (np.all(re == 0) and np.all(im == 0)):
			return re, im
		# if a solution failed: re=1, im=0 (XX/YY) or re=0, im=0 (XY,YX)
		# the closest good solution in applyed (interpolating if possible)
		# first positions
		if (re[0] == 1. and im[0] == 0.) or (re[0] == 0 and im[0] == 0):
			for i in xrange(len(re)):
				if (re[i] != 1. or im[i] != 0) and (re[i] != 0 or im[i] != 0):
					re[0:i]=re[i]
					im[0:i]=im[i]
					break
		# last positions
		if (re[-1] == 1. and im[-1] == 0) or (re[-1] == 0 and im[-1] == 0):
			for i in reversed(xrange(len(re))):
				if (re[i] != 1. or im[i] != 0) and (re[i] != 0 or im[i] != 0):
					re[i:]=re[i]
					im[i:]=im[i]
					break
		# intermediate positions -> interplate
		for i in xrange(len(re)):
			if (re[i] == 1. and im[i] == 0) or (re[i] == 0 and im[i] == 0):
				for j in xrange(i+1,len(re)):
					if (re[j] != 1. or im[j] != 0) and (re[j] != 0 or im[j] != 0):
						re_hi = re[j]
						im_hi = im[j]
						j_hi = j
						break
				for j in reversed(xrange(0,i)):
					if (re[j] != 1. or im[j] != 0) and (re[j] != 0 or im[j] != 0):
						re_low = re[j]
						im_low = im[j]
						j_low = j
						break
				for j in xrange(j_low+1,j_hi):
					re[j] = re_low+((re_hi-re_low)/(j_hi-j_low))*(j-j_low)
					im[j] = im_low+((im_hi-im_low)/(j_hi-j_low))*(j-j_low)
		return re, im

	# Recompute the phase setting "refant" as the reference antenna [TO BE IMPLEMENTED]
	def SetRef(self, refant):
		pass

####################################
# class that defines a calibrator SB
class Calib(MS, ParmdbR):
	def __init__(self, f):
		MS.__init__(self, f)
		ParmdbR.__init__(self, f+'/instrument')
		print "Activating calibrator:", self.GetName()

###########################################################
# class that defines a target SB where to transfer the data
class Target(MS):
	def __init__(self, f):
		MS.__init__(self, f)
		self.com = 'create tablename=\''+f+'/instrument\'\n'
		self.com += 'adddef Gain:0:0:Ampl values=1.0\n'
		self.com += 'adddef Gain:1:1:Ampl values=1.0\n'
		self.com += 'adddef Gain:0:0:Real values=1.0\n'
		self.com += 'adddef Gain:1:1:Real values=1.0\n'
		self.com += 'adddef DirectionalGain:0:0:Ampl values=1.0\n'
		self.com += 'adddef DirectionalGain:1:1:Ampl values=1.0\n'
		self.com += 'adddef DirectionalGain:0:0:Real values=1.0\n'
		self.com += 'adddef DirectionalGain:1:1:Real values=1.0\n'
		self.com += 'set stepx='+str(self.freq['step'])+'\n'
		self.com += 'set stepy='+str(self.time['step'])+'\n'
		self.amp = {}
		self.ph = {}

	# String formatter for parmdbm
	def FormatValues(self, arr):
	        string = '['
	        for i in xrange(len(arr)):
	          for j in xrange(len(arr[i])):
	            string += str(arr[i,j])+','
	        return string[0:-1]+']'

	# Add a new command to the command list for parmdbm
	def SetCom(self, amp, ph, corr, ant, domain):
	        real = amp*(np.cos(ph))
       		imag = amp*(np.sin(ph))
        	shapef = real.shape[0]
        	shapet = real.shape[1]
        	self.com += 'add Gain:'+corr+':Real:'+ant+' domain='+str(domain)+', shape=['+str(shapet)+', '+str(shapef)+'], values='+self.FormatValues(real)+'\n'
        	self.com += 'add Gain:'+corr+':Imag:'+ant+' domain='+str(domain)+', shape=['+str(shapet)+', '+str(shapef)+'], values='+self.FormatValues(imag)+'\n'

	# Run parmdbm using self.com as the input file
	def RunParmdb(self):
		print "Writing parmdb of "+self.GetName()
		self.com += 'quit'
		fw = file(self.GetName()+'-caltable', 'w')
		fw.write(self.com)
		fw.close()
		check_call(['parmdbm < '+self.GetName()+'-caltable >/dev/null 2>&1'], shell=True)
		check_call(['rm '+self.GetName()+'-caltable'], shell=True)

###########################################################
# class that defines a target SB where to update the data
class TargetU(Target, ParmdbR):
	def __init__(self, f):
		Target.__init__(self, f)
		ParmdbR.__init__(self, f+'/instrument')

        # Override Target.SetCom method to just perform amplitude rescaling
        def SetCom(self, amp, ph, corr, ant, domain):
		myamp, null = self.GetAmp(ant,corr)
		amp = myamp*(np.mean(amp)/np.mean(myamp))
		Target.SetCom(amp, ph, corr, ant, domain)

###########################################################
# class that defines a target SB where to read data
class SBr(MS, ParmdbR):
	def __init__(self, f):
		MS.__init__(self, f)
		ParmdbR.__init__(self, f+'/instrument')
