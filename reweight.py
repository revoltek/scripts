#!/usr/bin/env python

# Plot LOFAR weights
# Update LOFAR weights using residual visibilities
#
# Author: Francesco de Gasperin
# Credits: Frits Sweijen, Etienne Bonnassieux

import os, sys, logging, time
import numpy as np
from casacore.tables import taql, table
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

from lib_timer import Timer

class MShandler():
    def __init__(self, ms_files, wcolname, dcolname):
        """
        ms_files: a list of MeasurementSet files
        wcolname: name of the weight column
        dcolname: name of the residual column
        """
        logging.info('Reading: %s', ','.join(ms_files))
        logging.info('Weight column: %s', wcolname)
        logging.info('Data column: %s', dcolname)
        self.ms_files = ms_files
        self.wcolname = wcolname
        self.dcolname = dcolname

        # if more than one MS, virtualconcat
        #if len(ms_files) > 1:
        #    self.ms = taql('select TIME, FLAG, ANTENNA1, ANTENNA2, %s, %s from %s' % ( wcolname, dcolname, ','.join(ms_files) ))
        #else:
        #    self.ms = taql('select TIME, FLAG, ANTENNA1, ANTENNA2, %s, %s from %s' % ( wcolname, dcolname, ms_files[0] ))
        self.ms = table(ms_files[0], readonly=False, ack=False)

    def get_antennas(self):
        return taql('select NAME from %s/ANTENNA' % (self.ms_files[0]) )

    def get_freqs(self):
        # TODO: check if there is a smarter way to do it with concat.MS
        freqs = []
        for ms_file in self.ms_files:
            freqs += list( taql('SELECT CHAN_FREQ FROM %s/SPECTRAL_WINDOW' % ms_file)[0]['CHAN_FREQ'] * 1e-6 ) # in MHz
        return freqs

    def get_time(self):
        # TODO: fix if multiple time MSs passed
        ms_avgbl = taql('SELECT TIME FROM %s GROUPBY TIME' %(self.ms_files[0]))
        return ms_avgbl.getcol('TIME')

    def get_elev(self):
        # TODO: fix if multiple time MSs passed
        ms_avgbl = taql('SELECT TIME, MEANS(GAGGR(MSCAL.AZEL1()[1]), 0) AS ELEV FROM %s GROUPBY TIME' %(self.ms_files[0]))
        return ms_avgbl.getcol('ELEV')

    def iter_antenna(self, antennas=None):
        """
        Iterator to get all visibilities of each antenna
        It should return an array of Ntimes x Nfreqs
        """
        ms = self.ms # to use the "$" in taql

        for ant_id, ant_name in enumerate(self.get_antennas().getcol('NAME')):
            if antennas is not None and ant_name not in antennas: continue
            logging.info('Workign on antenna: %s', ant_name)
            yield ant_id, ant_name, taql('select TIME, GAGGR(FLAG) as GFLAG, ANTENNA1, ANTENNA2, GAGGR(%s) as GDATA, GAGGR(%s) as GWEIGHT \
                                          from $ms \
                                          where (ANTENNA1=%i or ANTENNA2=%i) and (ANTENNA1 != ANTENNA2) \
                                          groupby TIME' \
                    % (self.dcolname, self.wcolname, ant_id, ant_id) )


def reweight(MSh, mode):

    # get data per antenna
    var_antenna = {}
    med_antenna = {}
    for ant_id, ant_name, ms_ant in MSh.iter_antenna():

        with Timer('Get data'):
            data = ms_ant.getcol('GDATA') # axes: time, ant, freq, pol
            flags = ms_ant.getcol('GFLAG') # axes: time, ant, freq, pol

            # put flagged data to NaNs
            data[flags] = np.nan

            # if completely flagged set variance to 1 and continue
            if np.all(flags):
                var_antenna[ant_id] = None
                med_antenna[ant_id] = None
                continue

        with Timer('Prepare data'):

            # data column is updated subtracting adjacent channels
            if mode == 'subchan':
                data_shifted_l = np.roll(data, -1, axis=2)
                data_shifted_r = np.roll(data, +1, axis=2)
                # if only 2 freq it's aleady ok, subtracting one from the other
                if data.shape[2] > 2:
                    data_shifted_l[:,:,-1,:] = data_shifted_l[:,:,-3,:] # last chan uses the one but last
                    data_shifted_r[:,:,0,:] = data_shifted_r[:,:,2,:] # first chan uses third
                # get the "best" shift, either on the right or left. This is to avoid propagating bad channels (e.g. with RFI)
                ratio_l = np.nanvar(data_shifted_l, axis=(0,1,3))/np.nanmean(data_shifted_l, axis=(0,1,3))
                ratio_l[ np.isnan(ratio_l) ] = np.inf
                ratio_r = np.nanvar(data_shifted_l, axis=(0,1,3))/np.nanmean(data_shifted_l, axis=(0,1,3))
                ratio_r[ np.isnan(ratio_r) ] = np.inf
                data = np.where( ( ratio_l < ratio_r )[np.newaxis,np.newaxis,:,np.newaxis], data - data_shifted_l, data - data_shifted_r)

            # data column is updated subtracting adjacent times
            if mode == 'subtime':
                data_shifted_l = np.roll(data, -1, axis=0)
                data_shifted_r = np.roll(data, +1, axis=0)
                # if only 2 freq it's aleady ok, subtracting one from the other
                if data.shape[2] > 2:
                    data_shifted_l[-1,:,:,:] = data_shifted_l[-3,:,:,:] # last timeslot uses the one but last
                    data_shifted_r[0,:,:,:] = data_shifted_r[2,:,:,:] # first timeslot uses third
                # get the "best" shift, either on the right or left. This is to avoid propagating bad channels (e.g. with RFI)
                ratio_l = np.nanvar(data_shifted_l, axis=(1,2,3))/np.nanmean(data_shifted_l, axis=(1,2,3))
                ratio_l[ np.isnan(ratio_l) ] = np.inf
                ratio_r = np.nanvar(data_shifted_l, axis=(1,2,3))/np.nanmean(data_shifted_l, axis=(1,2,3))
                ratio_r[ np.isnan(ratio_r) ] = np.inf
                data = np.where( ( ratio_l < ratio_r )[:,np.newaxis,np.newaxis,np.newaxis], data - data_shifted_l, data - data_shifted_r)

            # use residual data, nothing to do here
            elif mode == 'residual':
                pass

        with Timer('Calc variances'):
            # find variance per time/freq for each antenna
            var_freqs = np.nanvar( data, axis=(1,2) ) # time x pol
            var_times = np.nanvar( data, axis=(0,1) ) # freq x pol
            var_antenna[ant_id] = var_freqs[:, np.newaxis]+var_times # sum of the time/freq variances - axes: time,freq,pol
            med_freqs = np.abs( np.nanmean( data, axis=(1,2) )**2 ) # time x pol
            med_times = np.abs( np.nanmean( data, axis=(0,1) )**2 ) # freq x pol
            med_antenna[ant_id] = med_freqs[:, np.newaxis]+med_times # sum of the time/freq mean - axes: time,freq,pol

    # reconstruct BL weights from antenna variance
    for ms_bl in MSh.ms.iter(["ANTENNA1","ANTENNA2"]):
        ant_id1 = ms_bl.getcol('ANTENNA1')[0]
        ant_id2 = ms_bl.getcol('ANTENNA2')[0]

        if var_antenna[ant_id1] is None or var_antenna[ant_id2] is None: continue

        w = 1.e-11/( var_antenna[ant_id1]*med_antenna[ant_id2] + var_antenna[ant_id2]*med_antenna[ant_id1] \
               + var_antenna[ant_id1]*var_antenna[ant_id2] )
        f = ms_bl.getcol('FLAG')
        # find how many unflagged weights are nans
        ntoflag = np.count_nonzero(np.isnan(w[~f]))
        logging.debug( 'BL: %i - %i: created %i new flags (%f%%)' % ( ant_id1, ant_id2, ntoflag, (100.*ntoflag)/np.size(w) ) )
        ms_bl.putcol(MSh.wcolname, w)
        ms_bl.flush()
        # flag weights that are nans
        taql('update $ms_bl set WEIGHT_SPECTRUM[isnan(WEIGHT_SPECTRUM)]=0')
        ms_bl.flush()

def plot(MSh, antennas):

    logging.info('Getting time/freq aggregated values...')
    freqs = MSh.get_freqs()
    elev = MSh.get_elev()
    time = MSh.get_time()
    time -= time[0]
    time /= 3600. # in h from the beginning of the obs

    fig = plt.figure(figsize=(15,15))

    for ant_id, ant_name, ms_ant in MSh.iter_antenna(antennas):
        
        fig.suptitle(ant_name, fontweight='bold')
        fig.subplots_adjust(wspace=0)
        #figgrid, axa = plt.subplots(Nr, Nc, sharex=True, sharey=True, figsize=figSize)
        axt = plt.subplot2grid((4, 2), (0, 0), colspan=2)
        axf = plt.subplot2grid((4, 2), (1, 0), colspan=2)
        ax1 = plt.subplot2grid((4, 2), (2, 0))
        ax2 = plt.subplot2grid((4, 2), (2, 1))
        ax3 = plt.subplot2grid((4, 2), (3, 0))
        ax4 = plt.subplot2grid((4, 2), (3, 1))

        ms_ant_avgbl = taql('SELECT MEANS(GAGGR(GWEIGHT[GFLAG]),1) AS WEIGHT, ALLS(GAGGR(GFLAG),1) as FLAG from $ms_ant') # return (1,time,freq,pol)

        w = ms_ant_avgbl.getcol('WEIGHT')[0] # axis: time, freq, pol
        # put flagged data to NaNs and skip if completely flagged
        if np.all(ms_ant_avgbl.getcol('FLAG')[0]):
            continue

        logging.info('Plotting %s...' % ant_name)

        w_f = np.nanmean(w, axis=0) # average in time
        w_t = np.nanmean(w, axis=1) # average in freq

        # 3D plot
        bbox = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        aspect = (time[-1]-time[0])*bbox.height/((freqs[-1]-freqs[0])*bbox.width)

        im = ax1.imshow(w[...,0].T, origin='lower', interpolation="none", cmap=plt.cm.jet, \
                        extent=[time[0],time[-1],freqs[0],freqs[-1]], aspect=str(aspect))#, vmin=0, vmax=1e-17)
        im = ax2.imshow(w[...,1].T, origin='lower', interpolation="none", cmap=plt.cm.jet, \
                        extent=[time[0],time[-1],freqs[0],freqs[-1]], aspect=str(aspect))#, vmin=0, vmax=1e-17)
        im = ax3.imshow(w[...,2].T, origin='lower', interpolation="none", cmap=plt.cm.jet, \
                        extent=[time[0],time[-1],freqs[0],freqs[-1]], aspect=str(aspect))#, vmin=0, vmax=1e-17)
        im = ax4.imshow(w[...,3].T, origin='lower', interpolation="none", cmap=plt.cm.jet, \
                        extent=[time[0],time[-1],freqs[0],freqs[-1]], aspect=str(aspect))#, vmin=0, vmax=1e-17)

        ax2.tick_params(labelleft='off')
        ax4.tick_params(labelleft='off')
        ax3.set_xlabel('Time [h]')
        ax4.set_xlabel('Time [h]')
        ax1.set_ylabel('Frequency [MHz]')
        ax3.set_ylabel('Frequency [MHz]')

        # Elevation
        ax_elev = axt.twinx()
        ax_elev.plot(time, elev * 180/np.pi, 'k', linestyle=':', linewidth=1, label='Elevation')
        ax_elev.set_ylabel('Elevation [deg]')

        axt.scatter(time, w_t[:,0], marker='.', alpha=0.25, color='red', label='XX Weights')
        axt.scatter(time, w_t[:,1], marker='.', alpha=0.25, color='green', label='XY Weights')
        axt.scatter(time, w_t[:,2], marker='.', alpha=0.25, color='orange', label='YX Weights')
        axt.scatter(time, w_t[:,3], marker='.', alpha=0.25, color='blue', label='YY Weights')
        axt.set_xlim(np.min(time), np.max(time))
        axt.set_ylim(np.nanmin(w_t), np.nanmax(w_t))
        axt.set_xlabel('Time [h]')
        #axt.set_ylabel('Weight')

        axf.scatter(freqs, w_f[:,0], marker='.', alpha=0.25, color='red', label='XX Weights')
        axf.scatter(freqs, w_f[:,1], marker='.', alpha=0.25, color='green', label='XY Weights')
        axf.scatter(freqs, w_f[:,2], marker='.', alpha=0.25, color='orange', label='YX Weights')
        axf.scatter(freqs, w_f[:,3], marker='.', alpha=0.25, color='blue', label='YY Weights')
        axf.set_xlim(np.min(freqs), np.max(freqs))
        axf.set_ylim(np.nanmin(w_f), np.nanmax(w_f))
        axf.set_xlabel('Frequency [MHz]')
        #axf.set_ylabel('Weight')

        handles, labels = axt.get_legend_handles_labels()
        handles2, labels2 = ax_elev.get_legend_handles_labels()
        leg = axt.legend(handles+handles2, labels+labels2, loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=5, borderaxespad=0.0)

        imagename = ant_name+'.png'
        logging.info('Save file: %s' % (imagename))
        fig.savefig(imagename, bbox_inches='tight', additional_artists=leg, dpi=250)
        fig.clf()

def readArguments():
    import argparse
    parser=argparse.ArgumentParser("Plot/Update weights for LOFAR")
    parser.add_argument("-v", "--verbose", help="Be verbose. Default is False", required=False, action="store_true")
    parser.add_argument("-p", "--plot", help="Plot the weights. Default is False", required=False, action="store_true")
    parser.add_argument("-m", "--mode", type=str, help="Mode can be: 'residual' if dcolname contain residual data; 'subchan'/'subtime' if adjacent-channel or adjacent-time subtraction has to be performed to remove the signal in dcolname. If not given do not update weights. Default: do not update", required=False, default=None)
    parser.add_argument("-d", "--dcolname", type=str, help="Name of the data column. Default: DATA.", required=False, default='DATA')
    parser.add_argument("-w", "--wcolname", type=str, help="Name of the weights column. Default: WEIGHT_SPECTRUM.", required=False, default="WEIGHT_SPECTRUM")
    parser.add_argument("-a", "--antennas", type=str, help="List of antennas to plot (comma separated). Default: all antennas", required=False, default=None)
    parser.add_argument("ms_files", type=str, help="MeasurementSet name(s).", nargs="+")
    args=parser.parse_args()
    return vars(args)

if __name__=="__main__":
    start_time = time.time()

    args         = readArguments()
    verbose      = args["verbose"]
    do_plot      = args["plot"]
    mode         = args["mode"]
    ms_files     = args["ms_files"]
    wcolname     = args["wcolname"]
    dcolname     = args["dcolname"]

    if len(ms_files) > 1:
        logging.error('More than 1 MS not implemented.')
        sys.exit()

    if args["antennas"] is not None:
        antennas = args["antennas"].replace(' ','').split(',')
    else: antennas = None

    if mode != 'residual' and mode != 'subchan' and mode != 'subtime' and mode is not None:
        logging.error('Unknown mode: %s' % mode)
        sys.exit()

    if verbose: logging.basicConfig(level=logging.DEBUG)
    else: logging.basicConfig(level=logging.INFO)

    logging.info('Reading MSs...')
    MSh = MShandler(ms_files, wcolname, dcolname)

    if mode is not None:
        logging.info('Computing weights...')
        reweight(MSh, mode)

    if do_plot:
        logging.info('Plotting...')
        plot(MSh, antennas)

    logging.debug('Running time %.0f s' % (time.time()-start_time))
