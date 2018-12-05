#!/usr/bin/env python

# Plot LOFAR weights
# Update LOFAR weights using residual visibilities
#
# Author: Francesco de Gasperin

import os, sys, logging, time
import numpy as np
from casacore.tables import taql
logging.basicConfig(level=logging.DEBUG)

class MShandler():
    def __init__(self, ms_files, wcolname, rcolname):
        """
        ms_files: a list of MeasurementSet files
        wcolname: name of the weight column
        rcolname: name of the residual column
        """
        logging.debug('Reading: %s', ','.join(ms_files))
        logging.debug('Weight column: %s', wcolname)
        logging.debug('Residual column: %s', rcolname)
        self.ms_files = ms_files

        # if more than one MS, virtualconcat
        if len(ms_files) > 1:
            self.ms = taql('select FLAG, ANTENNA1, ANTENNA2, %s, %s from %s where ANTENNA1 != ANTENNA2' % ( wcolname, rcolname, ','.join(ms_files) ))
        else:
            self.ms = taql('select FLAG, ANTENNA1, ANTENNA2, %s, %s from %s where ANTENNA1 != ANTENNA2' % ( wcolname, rcolname, ms_files[0] ))

    def iter_antenna(self):
        """
        Iterator to get all visibilities of each antenna
        It should return an array of Ntimes x Nfreqs
        """
        antennas = taql('select NAME from %s/ANTENNA' % self.ms[0])
        for antenna in antennas:
            yield taql('select FLAG, ANTENNA1, ANTENNA2, %s from $self.ms where ANTENNA1=%i or ANTENNA2=%i' % (rcolname, antenna, antenna) )

    def plot(self, time=None, freq=None):
        pass

def reweight(ms_files, wcolname='WEIGHT_SPECTRUM', rcolname='RESIDUAL_DATA'):

    # read residuals
    logging.info('Reading MSs...')
    MSh = MShandler(ms_files, wcolname, rcolname)

    # get residuals per antenna
    logging.info('Computing weights...')
    var_antenna = {}
    for ms_ant in MSh.iter_antenna():

        residuals = ms_ant.getcol(rcolname)

        # put flagged data to NaNs
        residuals[ms_ant.getcol('FLAG')] = np.nan

        # if completely flagged set variance to 1 and continue
        if np.isnan(residuals).all():
            var_antenna[ant] = 1.
            continue

        # find variance per time/freq for each antenna
        var_times = np.nanvar(residuals, axis=0)
        var_freqs = np.nanvar(residuals, axis=1)
        var_antenna[ant] = var_times[:, np.newaxis]+var_freqs # sum of the time/freq variances

    # reconstruct BL weights from antenna variance
    for ms_bl in MSh.ms.iter(["ANTENNA1","ANTENNA2"]):
        ant1 = ms_bl['ANTENNA1'][0]
        ant2 = ms_bl['ANTENNA2'][0]
        w = var_antenna[ant1] var_antenna[ant2]
        ms_bl[wcolname] = w

def readArguments():
    import argparse
    parser=argparse.ArgumentParser("Plot/Update weights for LOFAR")
    parser.add_argument("-v", "--verbose", help="Be verbose. Default is False", required=False, action="store_true")
    parser.add_argument("-w", "--wcolname", type=str, help="Name of the weights column. Default: WEIGHT_SPECTRUM.", required=False, default="WEIGHT_SPECTRUM")
    parser.add_argument("-r", "--rcolname", type=str, help="Name of the residuals column. Default: RESIDUAL_DATA.", required=False, default="RESIDUAL_DATA")
    parser.add_argument("ms_files", type=str, help="MeasurementSet name(s).", nargs="+")
    args=parser.parse_args()
    return vars(args)

if __name__=="__main__":
    args         = readArguments()
    ms_files     = args["ms_files"]
    wcolname     = args["wcolname"]
    rcolname     = args["rcolname"]

    start_time = time.time()
    reweight(ms_files, wcolname, rcolname)
    logging.debug('Running time %.0f s' % (time.time()-start_time))
