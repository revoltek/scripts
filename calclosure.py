#!/usr/bin/python

import os, sys, logging, itertools
import pyrap.tables as pt
import numpy as np

logging.basicConfig(level=logging.DEBUG)

ms = '3C295_SB193.MS' 
antRef = 0
plot = False
plotall = True
mode = 'triple'

def getPh(phase, antIdx, ant):
    """
    Get the phases relative to an antenna towards all the other antennas "Phi 1->2"
    The antenna is assumed to be at the first position of the BL ordinament, sign is corrected accordingly
    """
    phase[antIdx[1] == ant] *= -1 # account for "worng" BL direction
    p = phase[(antIdx[0] == ant) | (antIdx[1] == ant)]
    phase[antIdx[1] == ant] *= -1 # correct back the values
    return p


def getWe(weight, antIdx, ant):
    """
    Get the weight relative to an antenna
    """
    return weight[(antIdx[0] == ant) | (antIdx[1] == ant)]


def norm(phase):
    out = np.fmod(phase, 2. * np.pi)

    # Convert to range [-pi, pi]
    out[out < -np.pi] += 2. * np.pi
    out[out > np.pi] -= 2. * np.pi
    return out


def angMean(angs, weight):
    """
    Find the weighted mean of a series of angles
    """
    #assert len(angs) == len(weight)
    # normalization is unnecessary as we deal with just the angle
    return np.angle( np.sum( weight * np.exp(1j*np.array(angs)) ))# / ( len(angs) * sum(weight) ) )


if plot or plotall:
    import matplotlib as mpl
    mpl.rc('font',size =8 )
    mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure()
    fig.subplots_adjust(wspace=0)

# get antenna names
tant = pt.table(ms+'/ANTENNA', readonly=True, ack=False)
antNames = tant.getcol('NAME')
Nant = len(antNames)
tant.close()

# open ms and grab a timeslice
logging.info('Open table and fetch data')
tms = pt.table(ms, readonly=True, ack=False)
Ntime = len(set(tms.getcol('TIME')))

# array with solutions
solall = np.zeros( (Ntime,Nant), dtype=np.float64)
solall_t = np.zeros( (Ntime,Nant), dtype=np.float64)

for t, ts in enumerate(tms.iter('TIME')):
    logging.info('Working on time: '+str(t))
    time = ts.getcell('TIME',0) # shape: ant, chan, pol
    data = ts.getcol('SMOOTHED_DATA')[:,0,:] # shape: ant, chan, pol
    data_m = ts.getcol('MODEL_DATA')[:,0,:]
    ants1 = ts.getcol('ANTENNA1')
    ants2 = ts.getcol('ANTENNA2')
    ants = np.array(list(set(ants1)))
    antIdx = np.array([ants1,ants2])
    flags = ts.getcol('FLAG')[:,0,:]
    weight = ts.getcol('WEIGHT_SPECTRUM')[:,0,:]

    # scalar
    #data_amp = np.absolute(data[:,0])+np.absolute(data[:,3])
    #data_ph = norm( np.angle(data[:,0])+np.angle(data[:,3]) )
    #data_ph_m = norm( np.angle(data_m[:,0])+np.angle(data_m[:,3]) )
    #flags = np.logical_or(flags[:,0],flags[:,3])
    #weight = ( weight[:,0] + weight[:,3] )/2.

    # single pol
    data_ph = norm( np.angle(data_m[:,0]) - np.angle(data[:,0]) ) # -ph_sol is the real phase error
    data_we = weight[:,0]
    data_we[flags[:,0] == True] = 0 # weight flagged data 0
    
    # TODO: if ref ant is flagged?

    # cycle on antenna to solve for
    for s, antSol in enumerate(ants):

        if antSol == antRef: continue # always gives 0
        
        #logging.info('Working on antenna: '+str(antSol))

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
            solall[t,s] = angMean(sols, sols_w)

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

        #print sols, sols_w
        #sys.exit(1)

        solall[t,s] = angMean( np.array(sols).flatten(), np.array(sols_w).flatten() )

        #logging.debug("Mean: "+str(solall[t,s]))

        if plot: 
            fig.clf()
            figgrid, ax = plt.subplots(1, 1, figsize=(13,10), sharex=True, sharey=True)
            if mode == 'double':
                ax.plot(ants[(sols_w != 0)], sols[(sols_w != 0)], 'bo')
                ax.plot(ants[(sols_w == 0)], sols[(sols_w == 0)], 'ro')
            elif mode == 'triple':
                ax.plot(xrange(len(sols)), sols, 'yo')
            ax.set_title("Antenna "+antNames[antSol])
            ax.plot([0,36],[solall[t,s],solall[t,s]], 'k-')
            ax.set_ylim(ymin=-np.pi, ymax=np.pi)
            ax.set_xlim(xmin=-1, xmax=36)
            logging.debug('Plotting %d_%02i.png' % (time, antSol))
            plt.savefig('%d_%02i.png' % (time, antSol), bbox_inches='tight')
         
    if t == 500: break

if plotall:
    time = sorted(set(tms.getcol('TIME'))) # shape: ant, chan, pol
    for a, ant in enumerate(antNames):
        fig.clf()
        figgrid, ax = plt.subplots(1, 1, figsize=(13,10), sharex=True, sharey=True)
        #ax = fig.add_subplot(111)
        ax.set_title("Antenna "+ant)
        ax.plot( solall[0:501,a], 'ro', markersize=3)
        #ax.plot( time, solall[:,a], 'ro')
        ax.set_ylim(ymin=-np.pi, ymax=np.pi)
        ax.set_xlim(xmin=0, xmax=550)
        logging.debug('Plotting '+ant+'.png')
        plt.savefig(ant+'.png', bbox_inches='tight')

tms.close()







