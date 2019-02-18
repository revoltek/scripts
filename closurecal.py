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

import os, sys, logging, itertools
import pyrap.tables as pt
import numpy as np
import scipy.optimize

logging.basicConfig(level=logging.DEBUG)

ms = sys.argv[1] 
antRef = 0
plotph = False
plotamp = False
plotavg = False
plotall = True
plotTEC = False
mode = 'double'
timeavg = 1
freqavg = 4
solvetec = False

def getPh(phase, antIdx, ant):
    """
    Get the phases relative to an antenna towards all the other antennas "Phi 1->2"
    The antenna is assumed to be at the first position of the BL ordinament, sign is corrected accordingly
    """
    phase[antIdx[1] == ant] *= -1 # account for "worng" BL direction
    p = phase[(antIdx[0] == ant) | (antIdx[1] == ant)]
    phase[antIdx[1] == ant] *= -1 # correct back the values
    return p


def getAmp(amp, antIdx, ant, ant2 = None):
    """
    Get the amps relative to an antenna towards all the other antennas "Lambda_12"
    if ant2 != None: then return only that BL
    """
    if ant2 == None:
        return amp[(antIdx[0] == ant) | (antIdx[1] == ant)]
    else:
        return amp[((antIdx[0] == ant) & (antIdx[1] == ant2)) | ((antIdx[0] == ant2) & (antIdx[1] == ant))]


def getWe(weight, antIdx, ant, ant2 = None):
    """
    Get the weight relative to an antenna
    if ant2 != None: then return only that BL
    """
    if ant2 == None:
        return weight[(antIdx[0] == ant) | (antIdx[1] == ant)]
    else:
        return weight[((antIdx[0] == ant) & (antIdx[1] == ant2)) | ((antIdx[0] == ant2) & (antIdx[1] == ant))]


def norm(phase):
    """
    Normalize phases in [-pi, +pi]
    """
    out = np.fmod(phase, 2. * np.pi)

    # Convert to range [-pi, pi]
    out[out < -np.pi] += 2. * np.pi
    out[out > np.pi] -= 2. * np.pi
    return out


def angMean(angs, weights):
    """
    Find the weighted mean of a series of angles
    """
    #assert len(angs) == len(weight)
    # normalization is unnecessary as we deal with just the angle
    return np.angle( np.sum( weights * np.exp(1j*np.array(angs)) ))# / ( len(angs) * sum(weight) ) )


def angRMS(angs, weights):
    """
    Find the weighted rms of a series of angles
    """
    diff = angs - angMean(angs, weights)
    diff[diff < -np.pi] += 2*np.pi
    diff[diff > np.pi] -= 2*np.pi
    return np.sqrt( angMean(diff**2, weights) ) # weighted std dev


def findtec(phases, weights, freq, time, ant):
    """
    Find tec
    time, ant are just for plotting purposes
    """
    # TODO: add weights
    par1complex = lambda p, freq, phases, weights: ( abs(np.cos(8.44797245e9*p[0]/freq) - np.cos(phases)) +\
                                            abs(np.sin(8.44797245e9*p[0]/freq) - np.sin(phases)) ) * weights
#    par2complex = lambda p, freq, phases, weights: ( abs(np.cos(2.*8.44797245e9*p[0]/freq + p[1]) - np.cos(phases)) +\
#                                                     abs(np.sin(2.*8.44797245e9*p[0]/freq + p[1]) - np.sin(phases)) ) * weights
    fitresult, success = scipy.optimize.leastsq(par1complex, [0], args=(freq, phases, weights), maxfev=10000)
#    fitresult, success = scipy.optimize.basinhopping(par1complex, [0], T=1., minimizer_kwargs={'args':(freq, phases, weights)})
#    fitresult, success = scipy.optimize.leastsq(par1complex, [0,0], args=(freq, phases, weights), maxfev=10000)
    logging.debug("leastsqr"+str(fitresult))
    if plotTEC:
        fig.clf()
        ax = fig.add_subplot(110)
        fitfuncfastplot = lambda p, freq: np.mod(8.44797245e9*p[0]/freq + 1.*np.pi, 2.*np.pi) - np.pi
        ax.plot(freq, np.mod(phases + np.pi, 2.*np.pi) - np.pi, 'or' )
        TEC = np.mod((-8.44797245e9*(fitresult[0])/freq)+np.pi, 2*np.pi) - np.pi
        residual = np.mod(phases-TEC+np.pi,2.*np.pi)-np.pi
        ax.plot(freq, residual, '.', color='yellow')
        ax.plot(freq, fitfuncfastplot(fitresult, freq), "r-")
        plt.savefig(ant+'_T'+str(time)+'.png')
    return fitresult[0]


if plotph or plotamp or plotavg or plotall:
    import matplotlib as mpl
    mpl.rc('font',size =8 )
    mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22)
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    cmap = cm.get_cmap('Spectral')
    fig = plt.figure()
    fig.subplots_adjust(wspace=0)


# get antenna names
logging.info('Get antenna names')
tant = pt.table(ms+'/ANTENNA', readonly=True, ack=False)
antNames = tant.getcol('NAME')
Nant = len(antNames)
tant.close()

# get freq
logging.info('Get frequency information')
tspw = pt.table(ms+'/SPECTRAL_WINDOW', ack=False)
chans = tspw.getcol('CHAN_FREQ')[0]
Nfreq = len(chans)
assert Nfreq%freqavg == 0
tspw.close()
if solvetec: freqavg = Nfreq

logging.info('Open table and fetch data')
tms = pt.table(ms, readonly=True, ack=False)

# get time
Ntime = len(set(tms.getcol('TIME')))
assert Ntime%timeavg == 0

# array with solutions
solall = {'amp':np.zeros( (Ntime/timeavg,Nfreq/freqavg,Nant), dtype=np.float64), 'phase':np.zeros( (Ntime/timeavg,Nfreq/freqavg,Nant), dtype=np.float64)}
# in these array I store the solution at each time and freq, every timeavg times then I combine them. I need to store all the frequencies.
solsblock = {'amp':np.zeros( (timeavg,Nfreq,Nant), dtype=np.float64), 'phase':np.zeros( (timeavg,Nfreq,Nant), dtype=np.float64)}
solsblock_w = {'amp':np.zeros( (timeavg,Nfreq,Nant), dtype=np.float64), 'phase':np.zeros( (timeavg,Nfreq,Nant), dtype=np.float64)}

for t, ts in enumerate(tms.iter('TIME')):
    logging.info('Working on time: '+str(t))
    time = ts.getcell('TIME',0)
    # shape: ant, chan, pol
    weight = ts.getcol('WEIGHT_SPECTRUM')
    flags = ts.getcol('FLAG')
    weight[flags == True] = 0 # weight flagged data 0
    data = ts.getcol('SMOOTHED_DATA')
    data[ weight == 0 ] = 1. # remove nans
    data_m = ts.getcol('MODEL_DATA')
    ants1 = ts.getcol('ANTENNA1')
    ants2 = ts.getcol('ANTENNA2')
    ants = np.array(list(set(ants1)))
    antIdx = np.array([ants1,ants2])

    for f, freq in enumerate(chans):
        logging.info('Working on freq: '+str(f))

        # scalar
        #data_amp = np.absolute(data[:,f,0])+np.absolute(data[:,f,3])
        #data_ph = norm( np.angle(data[:,f,0])+np.angle(data[:,f,3]) )
        #data_ph_m = norm( np.angle(data_m[:,f,0])+np.angle(data_m[:,f,3]) )
        #weight = ( weight[:,f,0] + weight[:,f,3] )/2. # note that flags are not propagated in pol
    
        # single pol
        data_ph = norm( np.angle(data_m[:,f,0]) - np.angle(data[:,f,0]) )
        data_amp = np.abs( data_m[:,f,0] ) / np.abs ( data[:,f,0]  )
        data_we = weight[:,f,0]
        
        # TODO: if ref ant is flagged?
    
        # cycle on antenna to solve for
        for s, antSol in enumerate(ants):
            #logging.info('Working on antenna: '+str(antSol))
    
            if antSol != antRef: # leave 0 in the solutions
    
                # PHASES
                if mode == 'double':
                    # double closure
                    ph_ref = getPh(data_ph, antIdx, antRef)
                    ph_sol = getPh(data_ph, antIdx, antSol)
                    # (ph_ref - ph_1) - (ph_sol - ph_1)
                    sols = norm( ph_ref - ph_sol )
        
                    # calculate weights
                    we_ref = getWe(data_we, antIdx, antRef)
                    we_sol = getWe(data_we, antIdx, antSol)
                    sols_w = (we_ref + we_sol ) /2.
        
                    # if antSol = ant1: p_rs + p_ss = p_rs (single, remove)
                    sols[antSol] = 0
                    sols_w[antSol] = 0
                    # if antRef = ant1: p_rr + p_rs = p_rs (single, keep)
                    sols_w[antRef] = we_ref[antSol] # autocorr gives 0 weight, no /2
        
                elif mode == 'triple':
                    ph_ref = getPh(data_ph, antIdx, antRef)
                    ph_sol = getPh(data_ph, antIdx, antSol)
                    we_ref = getWe(data_we, antIdx, antRef)
                    we_sol = getWe(data_we, antIdx, antSol)
                    # triple closure
                    sols = []
                    sols_w = []
                    for at, ant2 in enumerate(ants):
                        if ant2 == antRef: continue # p_r1 + p_1r + p_rs = p_rs (single)
                        if ant2 == antSol: continue # p_r1 + p_1s + p_ss = p_r1 + p_1s (double with 1)
                        # if ant1 == ant2: fall back in double -> p_r1 + p_11 + p_1s = p_r1 + p_1s (double with 1==2, keep)
                        
                        # (ph_ref - ph_2) + (ph_2 - ph_1) - (ph_sol - ph_1)
                        ph_tri = getPh(data_ph, antIdx, ant2)
                        sols.append( norm( ph_ref[at] + ph_tri - ph_sol ) )
                        we_tri = getWe(data_we, antIdx, antRef)
                        sols_w.append( (we_ref[at] + we_sol + we_tri ) /3. )
        
                if not (np.array(sols_w).flatten() == 0).all():
                    avg = angMean( np.array(sols).flatten(), weights=np.array(sols_w).flatten() ) # weighted angular mean
                    solsblock['phase'][t%timeavg,f,s] = avg 
                    solsblock_w['phase'][t%timeavg,f,s] = 1./angRMS( np.array(sols).flatten(), np.array(sols_w).flatten() ) # weighted std dev

                # Debug plots
                if plotph and ( antNames[antSol] == 'CS002LBA' or antNames[antSol] == 'RS310LBA' or antNames[antSol] == 'RS106LBA' ):
                    fig.clf()
                    ax = fig.add_subplot(111)
                    ax.plot(range(len(sols)), sols, 'ro')
                    ax.set_title( "Antenna "+antNames[antSol]+" rms: "+str(1./solsblock_w['phase'][t%timeavg,f,s]) )
                    ax.plot([0,36],[solsblock['phase'][t%timeavg,f,s],solsblock['phase'][t%timeavg,f,s]], 'k-')
                    ax.set_ylim(ymin=-np.pi, ymax=np.pi)
                    ax.set_xlim(xmin=-1, xmax=36)
                    logging.debug('Plotting ph_T%d_F%d_%s.png' % (time, freq, antNames[antSol]))
                    plt.savefig('ph_T%d_F%d_%s.png' % (time, freq, antNames[antSol]), bbox_inches='tight')

            if solvetec : continue # skip amp if TEC solve
        
            # AMPLITUDES
            # a1S*aS3/a13 = e1 eS eS e2 / e1 e2 = e2**2
            # TODO: convert to log space
            amp_sol = getAmp(data_amp, antIdx, antSol)
            we_sol = getWe(data_we, antIdx, antSol)
            sols = []
            sols_w = []
            for ant1 in ants:
                if ant1 == antSol: continue # skip if 1==S
                amp_1 = getAmp(data_amp, antIdx, ant1)
                we_1 = getWe(data_we, antIdx, ant1)
                amp_1S = getAmp(data_amp, antIdx, ant1, ant2=antSol)
                we_1S = getWe(data_we, antIdx, ant1, ant2=antSol)
                sols.append(1./np.sqrt(amp_1S * amp_sol / amp_1))
                sols_w.append( (we_1S + we_sol + we_1) /3.)
    
                # if any antenna of the closure relation is flagged or an autocorrelation, set the weight to 0
                sols_w[-1][ we_1 == 0 ] = 0
                sols_w[-1][ we_sol == 0 ] = 0
                if we_1S == 0: sols_w[-1] = np.zeros_like(sols_w[-1])
    
            if not (np.array(sols_w).flatten() == 0).all():
                sols = np.log10(sols) # for amplitude work in log space
                avg = np.average( np.array(sols).flatten(), weights=np.array(sols_w).flatten() ) # weighted avg
                solsblock['amp'][t%timeavg,f,s] = avg 
                solsblock_w['amp'][t%timeavg,f,s] = 1./np.sqrt( np.average( ( np.array(sols).flatten() - avg )**2,\
                    weights=np.array(sols_w).flatten()) ) # weighted std dev

            # Debug plots
            if plotamp and ( antNames[antSol] == 'CS002LBA' or antNames[antSol] == 'RS310LBA' or antNames[antSol] == 'RS106LBA' ):
                fig.clf()
                ax = fig.add_subplot(111)
                for a in range(len(sols)):
                    ax.plot(range(len(sols[a][(sols_w[a] != 0)])), sols[a][(sols_w[a] != 0)], 'bo')
                ax.set_title( "Antenna "+antNames[antSol]+" rms: "+str(1./solsblock_w['amp'][t%timeavg,f,s]) )
                ax.plot([0,36],[solsblock['amp'][t%timeavg,f,s],solsblock['amp'][t%timeavg,f,s]], 'k-')
                logging.debug('Plotting amp_T%d_F%d_%s.png' % (time, freq, antNames[antSol]))
                plt.savefig('amp_T%d_F%d_%s.png' % (time, freq, antNames[antSol]), bbox_inches='tight')
    
        # end freq cycle

    # save actual solutions by re-averaging inside the freq/time steps
    if (t+1) % timeavg == 0:
        for s in range(Nant):

            if plotavg: 
                fig.clf()

            for f in range(Nfreq/freqavg):
                if solvetec:
                    solall['phase'][t/timeavg,f,s] = findtec( solsblock['phase'][:,f*freqavg:(f+1)*freqavg,s].flatten(),\
                        weights=solsblock_w['phase'][:,f*freqavg:(f+1)*freqavg,s].flatten(), freq=chans, time = t/timeavg, ant = antNames[s])
                else:
                    solall['phase'][t/timeavg,f,s] = angMean( solsblock['phase'][:,f*freqavg:(f+1)*freqavg,s].flatten(),\
                        weights=solsblock_w['phase'][:,f*freqavg:(f+1)*freqavg,s].flatten() )
                    # convert back from log space
                    solall['amp'][t/timeavg,f,s] = 10**np.average( solsblock['amp'][:,f*freqavg:(f+1)*freqavg,s].flatten(),\
                        weights=solsblock_w['amp'][:,f*freqavg:(f+1)*freqavg,s].flatten() )

                # Debug plots
                # color: freq, xaxis: time, table: ant
                if plotph or plotamp: 
                    times = list(range(solsblock['amp'].shape[0]))
                    ax = fig.add_subplot(121)
                    ax.set_title("PHASE - Antenna "+antNames[s])
                    ax.set_xlim(xmin=-0.5, xmax=len(times)-0.5)
                    for i in range(f*freqavg,(f+1)*freqavg):
                        ax.errorbar(times, solsblock['phase'][:,i,s], yerr=1./solsblock_w['phase'][:,i,s], c=cmap(float(i)/freqavg), fmt='o')
                    ax.plot([times[0],times[-1]], [solall['phase'][t/timeavg,f,s], solall['phase'][t/timeavg,f,s]], 'k-')
            
                    ax = fig.add_subplot(122)
                    ax.set_title("AMP - Antenna "+antNames[s])
                    for i in range(f*freqavg,(f+1)*freqavg):
                        ax.errorbar(times, solsblock['amp'][:,i,s], yerr=1./solsblock_w['amp'][:,i,s], c=cmap(float(i)/freqavg), fmt='o')
                    ax.plot([times[0],times[-1]], np.log10([solall['amp'][t/timeavg,f,s], solall['amp'][t/timeavg,f,s]]), 'k-')

                    logging.debug('Plotting Fin_T%d_F%d_%s.png' % (t/timeavg, f, antNames[s]))
                    plt.savefig('Fin_T%d_F%d_%s.png' % (t/timeavg, f, antNames[s]), bbox_inches='tight')

    #if t == 20: break
    # end time cycle

if plotall:
    for a, ant in enumerate(antNames):
        fig.clf()
        ax = fig.add_subplot(211)
        if solvetec:
            ax.set_title("TEC - Antenna "+ant)
            ax.plot( solall['phase'][:,:,a], 'o', markersize=3 )
        else:
            ax.set_title("PHASE - Antenna "+ant)
            ax.plot( solall['phase'][:,:,a], 'o', markersize=3 )
            ax = fig.add_subplot(212)
            ax.set_title("AMP - Antenna "+antNames[s])
            ax.plot( solall['amp'][:,:,a], '-', markersize=3 )
        logging.debug('Plotting '+ant+'.png')
        plt.savefig(ant+'.png', bbox_inches='tight')

tms.close()
