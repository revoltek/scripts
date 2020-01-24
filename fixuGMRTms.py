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

"""
Francesco de Gasperin
based on work by Martijn Oei

Adapt the MS format of uGMRT data to one usable by LOFAR software.
"""

import os, sys, logging, time
import numpy as np
from casacore import tables


class MS(object):

    def __init__(self, ms_file):
        logging.info("Starting work on MS at '" + ms_file + "'...")

        self.ms_file = ms_file
        self.t       = tables.table(ms_file,                      readonly = False, ack = False)
        self.tpol    = tables.table(ms_file + "/POLARIZATION",    readonly = False, ack = False)
        self.tspect  = tables.table(ms_file + "/SPECTRAL_WINDOW", readonly = False, ack = False)


    def close(self):
        '''
        Close tables opened in '__init__'.
        '''
        logging.info("Closing tables...")

        self.t.close()
        self.tpol.close()
        self.tspect.close()


    def columnExists(self, columnName):
        '''
        Check whether a column with name 'columnName' exists.
        '''
        columnNames = self.t.colnames()

        return (columnName in columnNames)


    def removeColumns(self):
        '''
        Remove columns that are never used by the LOFAR software, and are thus a waste of disk space (e.g. "SIGMA_SPECTRUM"),
        or that are generated (in the right shape) later in the pipeline (e.g. "MODEL_DATA").
        Note: removal of columns can give errors when executing LOFAR command 'msoverview'.
        '''
        logging.info("- Removal of unnecessary data columns -")

        for columnName in ["SIGMA_SPECTRUM", "MODEL_DATA"]: # This list could possibly be expanded.
            if (self.columnExists(columnName)):
                self.t.removecols(columnName)


    def updatePolarisation(self):
        '''
        Make sure that the MS contains 4 polarisations.
        "CORR_TYPE"    column description comment: 'The polarization type for each correlation product, as a Stokes enum.'
        "CORR_PRODUCT" column description comment: 'Indices describing receptors of feed going into correlation'
        '''
        logging.info("- Adaptation of polarisation metadata -")

        correlationTypesNew      = np.array([[5, 6, 7, 8]])
        correlationProductsNew   = np.array([[[0, 0], [0, 1], [1, 0], [1, 1]]])
        numberOfCorrelationsNew  = 4

        self.tpol.putcol("CORR_TYPE",    correlationTypesNew)
        self.tpol.putcol("CORR_PRODUCT", correlationProductsNew)
        self.tpol.putcol("NUM_CORR",     numberOfCorrelationsNew)


    def updateFreqMetadata(self):
        '''
        Flip the frequency channel order to match the LOFAR convention of ascending frequencies.
        We first check whether the flip is necessary at all, to make the code robust to running this program twice on the same MS.
        '''
        frequencies = self.tspect.getcol("CHAN_FREQ")

        if (frequencies[0, 0] > frequencies[0, -1]):
            logging.info("- Adaptation of frequency metadata -")
            self.tspect.putcol("CHAN_FREQ", np.fliplr(frequencies))
            return True
        else:
            logging.info("Frequency order already correct.")
            return False


    def updateFieldMetadata(self):
        '''
        MSs that were originally multi-field, have been split up in single-field MSs.
        Adapt the field metadata accordingly.
        '''
        logging.info("- Adaptation of field information -")

        pathMS      = self.ms_file
        pathMSField = self.ms_file + "/FIELD"

        # Remove metadata of other fields in the FIELD subtable.
        tables.taql("delete from $pathMSField where rownr() not in (select distinct FIELD_ID from $pathMS)")

        # Set 'SOURCE_ID' to 0 in the FIELD subtable.
        tables.taql("update $pathMSField set SOURCE_ID=0")

        # Set 'FIELD_ID' to 0 in the main table.
        tables.taql("update $pathMS set FIELD_ID=0")


    def updateIntervals(self):
        '''
        Update the INTERVAL and TIME columns,
        so that differences between time stamps are always integer multiples of a fixed constant.
        '''
        logging.info("- Adaptation of time intervals -")

        # Calculate the new interval.
        pathMS          = self.ms_file
        times           = (tables.taql("select distinct TIME from $pathMS")).getcol("TIME")
        intervals       = times[1 : ] - times[ : -1]
        intervals       = intervals[intervals < 1.5 * np.min(intervals)] # Select only intervals that do not correspond with big jumps in time.
        intervalPrecise = np.mean(intervals)

        # Update INTERVAL column.
        intervalsOld    = self.t.getcol("INTERVAL")
        intervalsNew    = np.ones_like(intervalsOld) * intervalPrecise
        self.t.putcol("INTERVAL", intervalsNew)

        # Update TIME column.
        timesOld        = self.t.getcol("TIME")
        timesNew        = timesOld[0] + np.round((timesOld - timesOld[0]) / intervalPrecise, 0) * intervalPrecise
        self.t.putcol("TIME", timesNew)


        logging.info("Interval set to %f" % intervalPrecise)
        logging.debug("Time intervals (old; possibly unequal):")
        logging.debug(np.unique(timesOld)[1 : ] - np.unique(timesOld)[ : -1])
        logging.debug("Time intervals (new; should be equal):")
        logging.debug(np.unique(timesNew)[1 : ] - np.unique(timesNew)[ : -1])


    def updateColumns(self, updateFreq):
        '''
        Update DATA, FLAG and WEIGHT_SPECTRUM columns.
        '''
        logging.info("- Change existing (or create alternative) columns for data, flags and weights -")

        logging.info("Loading visibilities...")
        visibilities      = self.t.getcol("DATA")

        if (updateFreq): # If the frequency channel order was flipped in 'updateFreqMetadata', also flip the data.
            visibilities = visibilities[ : , : : -1, : ]

        if (visibilities.shape[2] == 2):
            visibilitiesNew             = np.zeros((visibilities.shape[0], visibilities.shape[1], 4), dtype = np.complex128)
            visibilitiesNew[ : , : , 0] = visibilities[ : , : , 0]
            visibilitiesNew[ : , : , 3] = visibilities[ : , : , 1]
            visibilities                = visibilitiesNew

        keywordNames      = self.t.colkeywordnames("DATA")
        columnDescription = self.t.getcoldesc("DATA")
        dataManagerInfo   = self.t.getdminfo("DATA")

        dataManagerInfo["NAME"]                                  = "TiledDATAMartijn"
        dataManagerInfo["SPEC"]["DEFAULTTILESHAPE"]              = np.array([4, 40, 819],                                      dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["TileShape"] = np.array([4, 40, 819],                                      dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CubeShape"] = np.array([4, visibilities.shape[1], visibilities.shape[0]], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CellShape"] = np.array([4, visibilities.shape[1]],                        dtype = np.int32)

        #logging.debug("keywordNames:")
        #logging.debug(keywordNames)
        #logging.debug("columnDescription:")
        #logging.debug(columnDescription)
        #logging.debug("dataManagerInfo:")
        #logging.debug(dataManagerInfo)
        #logging.debug("dataManagerInfo (updated):")
        #logging.debug(dataManagerInfo)

        logging.info("Removing column 'DATA', if it exists...")
        if (self.columnExists('DATA')):
            self.t.removecols('DATA')

        logging.info("Adding column 'DATA'...")
        self.t.addcols(tables.makecoldesc('DATA', columnDescription), dataManagerInfo)

        logging.info("Filling column 'DATA'...")
        self.t.putcol('DATA', visibilities)


        logging.info("Loading flags...")
        flags             = self.t.getcol("FLAG")

        if (updateFreq): # If the frequency channel order was flipped in 'updateFreqMetadata', also flip the data.
            flags = flags[ : , : : -1, : ]

        if (flags.shape[2] == 2):
            flagsNew          = np.zeros((flags.shape[0], flags.shape[1], 4), dtype = np.bool_)
            flagsNew[ : , : , 0] = flags[ :, : , 0]
            flagsNew[ : , : , 1] = flags[ :, : , 0] # Take over flags from LL correlation
            flagsNew[ : , : , 2] = flags[ :, : , 0] # Take over flags from LL correlation
            flagsNew[ : , : , 3] = flags[ :, : , 1]
            flags = flagsNew

        keywordNames      = self.t.colkeywordnames("FLAG")
        columnDescription = self.t.getcoldesc("FLAG")
        dataManagerInfo   = self.t.getdminfo("FLAG")

        dataManagerInfo["NAME"]                                  = "TiledFlagMartijn"
        dataManagerInfo["SPEC"]["DEFAULTTILESHAPE"]              = np.array([4, 40, 819],                        dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["TileShape"] = np.array([4, 40, 819],                        dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CubeShape"] = np.array([4, flags.shape[1], flags.shape[0]], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CellShape"] = np.array([4, flags.shape[1]],                 dtype = np.int32)

        #logging.debug("keywordNames:")
        #logging.debug(keywordNames)
        #logging.debug("columnDescription:")
        #logging.debug(columnDescription)
        #logging.debug("dataManagerInfo:")
        #logging.debug(dataManagerInfo)
        #logging.debug("dataManagerInfo (updated):")
        #logging.debug(dataManagerInfo)

        logging.info("Removing column 'FLAG', if it exists...")
        if (self.columnExists('FLAG')):
            self.t.removecols('FLAG')

        logging.info("Adding column 'FLAG'...")
        self.t.addcols(tables.makecoldesc('FLAG', columnDescription), dataManagerInfo)

        logging.info("Filling column 'FLAG'...")
        self.t.putcol('FLAG', flags)


        logging.info("Loading weights...")
        weights           = self.t.getcol("WEIGHT_SPECTRUM")

        if (updateFreq): # If the frequency channel order was flipped in 'updateFreqMetadata', also flip the data.
            weights = weights[ : , : : -1, : ]

        if (weights.shape[2] == 2):
            weightsNew               = np.zeros((weights.shape[0], weights.shape[1], 4), dtype = np.float64)
            weightsNew[ : , : , 0]   = weights[ : , : , 0]
            weightsNew[ : , : , 3]   = weights[ : , : , 1]
            weights                  = weightsNew

        keywordNames      = self.t.colkeywordnames("WEIGHT_SPECTRUM")
        columnDescription = self.t.getcoldesc("WEIGHT_SPECTRUM")
        dataManagerInfo   = self.t.getdminfo("WEIGHT_SPECTRUM")

        dataManagerInfo["NAME"]                                  = "TiledWgtSpectrumMartijn"
        dataManagerInfo["SPEC"]["DEFAULTTILESHAPE"]              = np.array([4, 40, 819],                            dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["TileShape"] = np.array([4, 40, 819],                            dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CubeShape"] = np.array([4, weights.shape[1], weights.shape[0]], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CellShape"] = np.array([4, weights.shape[1]],                   dtype = np.int32)

        #logging.debug("keywordNames:")
        #logging.debug(keywordNames)
        #logging.debug("columnDescription:")
        #logging.debug(columnDescription)
        #logging.debug("dataManagerInfo:")
        #logging.debug(dataManagerInfo)
        #logging.debug("dataManagerInfo (updated):")
        #logging.debug(dataManagerInfo)

        logging.info("Removing column 'WEIGHT_SPECTRUM', if it exists...")
        if (self.columnExists('WEIGHT_SPECTRUM')):
            self.t.removecols('WEIGHT_SPECTRUM')

        logging.info("Adding column 'WEIGHT_SPECTRUM'...")
        self.t.addcols(tables.makecoldesc('WEIGHT_SPECTRUM', columnDescription), dataManagerInfo)

        logging.info("Filling column 'WEIGHT_SPECTRUM'...")
        self.t.putcol('WEIGHT_SPECTRUM', weights)



def readArguments():
    import argparse
    parser = argparse.ArgumentParser("Adapt uGMRT MS to LOFAR MS format.")
    parser.add_argument("-v", "--verbose", help = "Be verbose. Default is False", required = False, action = "store_true")
    parser.add_argument("ms_files",        help = "MeasurementSet name(s).", type = str, nargs = "+")
    args   = parser.parse_args()
    return vars(args)



if (__name__ == "__main__"):
    start_time   = time.time()

    args         = readArguments()
    verbose      = args["verbose"]
    ms_files     = args["ms_files"]

    if (verbose):
        logging.basicConfig(level = logging.DEBUG)
    else:
        logging.basicConfig(level = logging.INFO)

    logging.info('Reading MSs...')
    MSs          = []
    for ms_file in ms_files:
        MSs.append(MS(ms_file))

    for MS in MSs:
        MS.removeColumns()
        MS.updatePolarisation()
        updateFreq = MS.updateFreqMetadata()
        MS.updateFieldMetadata()
        MS.updateIntervals()
        MS.updateColumns(updateFreq)
        MS.close()

    logging.debug('Running time %.0f s' % (time.time() - start_time))
    logging.info('Done.')
