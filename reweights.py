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

        self.antennas = taql('select NAME from %s/ANTENNA' % self.ms_files[0])

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

    def iter_antenna(self):
        """
        Iterator to get all visibilities of each antenna
        It should return an array of Ntimes x Nfreqs
        """
        ms = self.ms # to use the "$" in taql

        for ant_id, ant_name in enumerate(self.antennas.getcol('NAME')):
            logging.info('Workign on antenna: %s', ant_name)
            yield ant_id, ant_name, taql('select TIME, GAGGR(FLAG) as GFLAG, ANTENNA1, ANTENNA2, GAGGR(%s) as GDATA, GAGGR(%s) as GWEIGHT \
                                          from $ms \
                                          where (ANTENNA1=%i or ANTENNA2=%i) and (ANTENNA1 != ANTENNA2) \
                                          groupby TIME' \
                    % (self.dcolname, self.wcolname, ant_id, ant_id) )


def reweight(MSh, mode):

    # get data per antenna
    var_antenna = {}
    for ant_id, ant_name, ms_ant in MSh.iter_antenna():

        data = ms_ant.getcol('GDATA') # axes: time, ant, freq, pol
        print data.shape

        # put flagged data to NaNs
        data[ms_ant.getcol('GFLAG')] = np.nan

        # if completely flagged set variance to 1 and continue
        if np.isnan(data).all():
            var_antenna[ant_id] = None
            continue

        # data column is updated subtracting adjacent channels
        # in case of mode=residuals there's nothing to do
        # TODO: get the best var from the sub on the right and on the left - this shoudl avoid contamination from a bed channel
        if mode == 'subchan':
            data_shifted = np.roll(data, -1, axis=2)
            # if only 2 freq it's aleady ok, subtracting one from the other
            if data.shape[2] > 2:
                data_shifted[:,:,-1,:] = data_shifted[:,:,-3,:] # last chan uses the one but last (they switch)
            data -= data_shifted

        # find variance per time/freq for each antenna
        var_freqs = np.nanvar( data, axis=(1,2) ) # time x pol
        var_times = np.nanvar( data, axis=(0,1) ) # freq x pol
        var_antenna[ant_id] = var_freqs[:, np.newaxis]+var_times # sum of the time/freq variances - axes: time,freq,pol
        #print var_antenna[ant_id].shape

    # reconstruct BL weights from antenna variance
    for ms_bl in MSh.ms.iter(["ANTENNA1","ANTENNA2"]):
        ant_id1 = ms_bl.getcol('ANTENNA1')[0]
        ant_id2 = ms_bl.getcol('ANTENNA2')[0]

        if var_antenna[ant_id1] is None or var_antenna[ant_id2] is None: continue

        w = 1./(var_antenna[ant_id1] + var_antenna[ant_id2]) # TODO: what is the right equation here?
        logging.debug( 'nan count (BL: %i - %i): %i' % ( ant_id1, ant_id2, np.sum(np.isnan(w))) )
        ms_bl.putcol(MSh.wcolname, w)
        ms_bl.flush()

def plot(MSh):

    logging.info('Getting time/freq aggregated values...')
    freqs = MSh.get_freqs()
    elev = MSh.get_elev()
    time = MSh.get_time()
    time -= time[0]
    time /= 3600. # in h from the beginning of the obs

    fig = plt.figure(figsize=(15,10))

    for ant_id, ant_name, ms_ant in MSh.iter_antenna():
        
        fig.suptitle(ant_name, fontweight='bold')
        fig.subplots_adjust(wspace=0)
        axt = fig.add_subplot(211)
        axf = fig.add_subplot(212)

        ms_ant_avgbl = taql('SELECT MEANS(GAGGR(GWEIGHT),1) AS WEIGHT, ALLS(GAGGR(GFLAG),1) as FLAG from $ms_ant') # return (1,time,freq,pol)

        w = ms_ant_avgbl.getcol('WEIGHT')[0] # axis: time, freq, pol
        # put flagged data to NaNs and skip if completely flagged
        if (ms_ant_avgbl.getcol('FLAG')[0] == True).all():
            continue

        logging.info('Plotting %s...' % ant_name)

        w_f = np.nanmean(w, axis=0) # average in time
        w_t = np.nanmean(w, axis=1) # average in freq

        # Elevation
        ax_elev = axt.twinx()
        ax_elev.plot(time, elev * 180/np.pi, 'k', linestyle=':', linewidth=1, label='Elevation')
        ax_elev.set_ylabel('Elevation [deg]')

        axt.scatter(time, w_t[:,0], marker='.', alpha=0.25, color='red', label='XX Weights')
        axt.scatter(time, w_t[:,1], marker='.', alpha=0.25, color='green', label='XY Weights')
        axt.scatter(time, w_t[:,2], marker='.', alpha=0.25, color='orange', label='YX Weights')
        axt.scatter(time, w_t[:,3], marker='.', alpha=0.25, color='blue', label='YY Weights')
        axt.set_xlim(np.nanmin(time), np.nanmax(time))
        axt.set_ylim(np.nanmin(w_t), np.nanmax(w_t))
        axt.set_xlabel('Time [h]')
        #axt.set_ylabel('Weight')

        axf.scatter(freqs, w_f[:,0], marker='.', alpha=0.25, color='red', label='XX Weights')
        axf.scatter(freqs, w_f[:,1], marker='.', alpha=0.25, color='green', label='XY Weights')
        axf.scatter(freqs, w_f[:,2], marker='.', alpha=0.25, color='orange', label='YX Weights')
        axf.scatter(freqs, w_f[:,3], marker='.', alpha=0.25, color='blue', label='YY Weights')
        axf.set_xlim(np.nanmin(freqs), np.nanmax(freqs))
        axf.set_ylim(np.nanmin(w_f), np.nanmax(w_f))
        axf.set_xlabel('Frequency [MHz]')
        #axf.set_ylabel('Weight')

        handles, labels = axt.get_legend_handles_labels()
        handles2, labels2 = ax_elev.get_legend_handles_labels()
        leg = axt.legend(handles+handles2, labels+labels2, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=5, borderaxespad=0.0)

        imagename = ant_name+'.png'
        logging.info('Save file: %s' % (imagename))
        fig.savefig(imagename, bbox_inches='tight', additional_artists=leg, dpi=250)
        fig.clf()

def readArguments():
    import argparse
    parser=argparse.ArgumentParser("Plot/Update weights for LOFAR")
    parser.add_argument("-v", "--verbose", help="Be verbose. Default is False", required=False, action="store_true")
    parser.add_argument("-p", "--plot", help="Plot the weights. Default is False", required=False, action="store_true")
    parser.add_argument("-m", "--mode", type=str, help="Mode can be: 'residual' if dcolname contain residual data; 'subchan' if adjacent-channel subtraction has to be performed to remove the signal in dcolname. If not given do not update weights. Default: do not update", required=False, default=None)
    parser.add_argument("-d", "--dcolname", type=str, help="Name of the data column. Default: DATA.", required=False, default='DATA')
    parser.add_argument("-w", "--wcolname", type=str, help="Name of the weights column. Default: WEIGHT_SPECTRUM.", required=False, default="WEIGHT_SPECTRUM")
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

    if mode != 'residual' and mode != 'subchan' and mode is not None:
        logging.error('Unknown mode: %s' % mode)
        sys.exit()

    if verbose: logging.basicConfig(level=logging.DEBUG)

    logging.info('Reading MSs...')
    MSh = MShandler(ms_files, wcolname, dcolname)

    if mode is not None:
        logging.info('Computing weights...')
        reweight(MSh, mode)

    if do_plot:
        logging.info('Plotting...')
        plot(MSh)

    logging.debug('Running time %.0f s' % (time.time()-start_time))
