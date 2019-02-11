#!/usr/bin/env python

"""
Francesco de Gasperin
based on work by Martijn Oei

Adapt the MS format of uGMRT data to one usable by LOFAR software.
"""

import os, sys, logging, time
import numpy as np
from casacore import tables
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

class MS( object ):

    def __init__(self, ms_file):
        self.ms_file = ms_file
        self.t = tables.table(ms_file, readonly=False, ack=False)
        self.tpol = tables.table(ms_file+"/POLARIZATION", readonly=False, ack=False)
        self.tspect = tables.table(ms_file+"/SPECTRAL_WINDOW", readonly=False, ack=False)

        logging.info("Starting work on MS at '" + ms_file + "'...")

    def columnExists(self, columnName):
        '''
        Check whether a column with name 'columnName' exists.
        '''
        columnNames = self.t.colnames()
        return (columnName in columnNames)

    def removeColumns(self):
        # Removal of columns can give errors when executing LOFAR command 'msoverview'.
        logging.info("- Removal of unnecessary data columns -")

        # Old array: ["EXPOSURE", "SIGMA_SPECTRUM"]. The column "EXPOSURE" seems necessary for DPPP (and LOFAR's msoverview), though.
        for columnName in ["SIGMA_SPECTRUM"]: # This list could possibly be expanded.
            if (self.columnExists(columnName)):
                t.removecols(columnName)

    def updatePolarisation(self):
        # CORR_TYPE    column description comment: 'The polarization type for each correlation product, as a Stokes enum.'
        # CORR_PRODUCT column description comment: 'Indices describing receptors of feed going into correlation'
        logging.info("- Adaptation of polarisation metadata -")

        correlationTypesNew      = np.array([[5, 6, 7, 8]])
        correlationProductsNew   = np.array([[[0, 0], [0, 1], [1, 0], [1, 1]]])
        numberOfCorrelationsNew  = 4

        self.tpol.putcol("CORR_TYPE",    correlationTypesNew)
        self.tpol.putcol("CORR_PRODUCT", correlationProductsNew)
        self.tpol.putcol("NUM_CORR",     numberOfCorrelationsNew)

    def updateFreqMetadata(self):

        frequencies = self.tspect.getcol("CHAN_FREQ")

        # reorder if flipped
        if frequencies[0,0] > frequencies[0,-1]:
            logging.info("- Adaptation of frequency metadata -")
            self.tspect.putcol("CHAN_FREQ", np.fliplr(frequencies) )
            return True
        else:
            logging.info('Frequency already correct.')
            return False

    def updateFieldMetadata(self):
        logging.info("- Adaptation of field information -")

        pathMS = self.ms_file
        pathMSField = self.ms_file + "/FIELD"

        # Remove metadata of other fields in the FIELD subtable
        tables.taql("delete from $pathMSField where rownr() not in (select distinct FIELD_ID from $pathMS)")

        # Set 'SOURCE_ID' to 0 in the FIELD subtable
        tables.taql("update $pathMSField set SOURCE_ID=0")

        # Set 'FIELD_ID' to 0 in the main table
        tables.taql("update $pathMS set FIELD_ID=0")

    def updateIntervals(self):
        logging.info("- Adaptation of intervals -")

        pathMS = self.ms_file
        times = (tables.taql("select distinct TIME from $pathMS")).getcol("TIME")
        intervalPrecise = times[1] - times[0]

        # Open the MS as a table in a way that changes can be made.
        intervals       = self.t.getcol("INTERVAL")
        intervalsNew    = np.ones_like(intervals) * intervalPrecise
        self.t.putcol("INTERVAL", intervalsNew)

        logging.debug("Time intervals (should be equal):")
        logging.debug(times[1:] - times[:-1])

    def updateColumns(self, updateFreq):
        logging.info("- Change existing (or create alternative) columns for data, flags and weights -")

        # First, change the 'direction' of the frequency axis.
        # If the visibilities were sorted in descending frequency order (id est given for high frequencies first),
        # they are converted to ascending frequency order. Note: applying this operation twice returns the MS to its old state!
        logging.info("Loading visibilities...")
        visibilities = self.t.getcol("DATA")

        if updateFreq: visibilities = visibilities[:,::-1,:]
        if len(visibilities.shape) == 2:
            visibilitiesNew          = np.zeros((visibilities.shape[0], visibilities.shape[1], 4), dtype = np.complex128)
            visibilitiesNew[:, :, 0] = visibilities[:, :, 0]
            visibilitiesNew[:, :, 3] = visibilities[:, :, 1]
            visibilities = visibilitiesNew

        keywordNames             = self.t.colkeywordnames("DATA")
        columnDescription        = self.t.getcoldesc("DATA")
        dataManagerInfo          = self.t.getdminfo("DATA")

        logging.debug("keywordNames:")
        logging.debug(keywordNames)
        logging.debug("columnDescription:")
        logging.debug(columnDescription)
        logging.debug("dataManagerInfo:")
        logging.debug(dataManagerInfo)

        dataManagerInfo["NAME"]                                  = "TiledDATAFix"
        dataManagerInfo["SPEC"]["DEFAULTTILESHAPE"]              = np.array([4, 40, 819], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["TileShape"] = np.array([4, 40, 819], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CubeShape"] = np.array([4, visibilities.shape[1], visibilities.shape[0]], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CellShape"] = np.array([4, visibilities.shape[1]], dtype = np.int32)

        logging.debug("dataManagerInfo (updated):")
        logging.debug(dataManagerInfo)

        logging.info("Removing column 'DATA', if it exists...")
        if self.columnExists('DATA'):
            self.t.removecols('DATA')

        logging.info("Adding column 'DATA'...")
        self.t.addcols(tables.makecoldesc('DATA', columnDescription), dataManagerInfo)

        logging.info("Filling column 'DATA'...")
        self.t.putcol('DATA', visibilities)

        ###
        logging.info("Loading flags...")
        flags = self.t.getcol("FLAG")

        if updateFreq: flags = flags[:,::-1,:]
        if len(flags.shape) == 2:
            flagsNew                 = np.zeros((flags.shape[0], flags.shape[1], 4), dtype = np.bool_)
            flagsNew[:, :, 0]        = flags[:, :, 0]
            flagsNew[:, :, 1]        = flags[:, :, 0] # Take over flags from LL correlation
            flagsNew[:, :, 2]        = flags[:, :, 0] # Take over flags from LL correlation
            flagsNew[:, :, 3]        = flags[:, :, 1]
            flags = flagsNew

        keywordNames             = self.t.colkeywordnames("FLAG")
        columnDescription        = self.t.getcoldesc("FLAG")
        dataManagerInfo          = self.t.getdminfo("FLAG")

        logging.debug("keywordNames:")
        logging.debug(keywordNames)
        logging.debug("columnDescription:")
        logging.debug(columnDescription)
        logging.debug("dataManagerInfo:")
        logging.debug(dataManagerInfo)

        dataManagerInfo["NAME"]                                  = "TiledFlagMartijn"
        dataManagerInfo["SPEC"]["DEFAULTTILESHAPE"]              = np.array([4, 40, 819], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["TileShape"] = np.array([4, 40, 819], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CubeShape"] = np.array([4, flags.shape[1], flags.shape[0]], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CellShape"] = np.array([4, flags.shape[1]], dtype = np.int32)

        logging.debug("dataManagerInfo (updated):")
        logging.debug(dataManagerInfo)

        logging.info("Removing column 'FLAG', if it exists...")
        if self.columnExists('FLAG'):
            self.t.removecols('FLAG')

        logging.info("Adding column 'FLAG'...")
        self.t.addcols(tables.makecoldesc('FLAG', columnDescription), dataManagerInfo)

        logging.info("Filling column 'FLAG'...")
        self.t.putcol('FLAG', flags)

        ###
        logging.info("Loading weights...")
        weights = self.t.getcol("WEIGHT_SPECTRUM")

        if updateFreq: weights = weights[:,::-1,:]
        if len(flags.shape) == 2:
            weightsNew               = np.zeros((weights.shape[0], weights.shape[1], 4), dtype = np.float64)
            weightsNew[:, :, 0]      = weights[:, :, 0]
            weightsNew[:, :, 3]      = weights[:, :, 1]
            weights = weightsNew

        keywordNames             = self.t.colkeywordnames("WEIGHT_SPECTRUM")
        columnDescription        = self.t.getcoldesc("WEIGHT_SPECTRUM")
        dataManagerInfo          = self.t.getdminfo("WEIGHT_SPECTRUM")

        logging.debug("keywordNames:")
        logging.debug(keywordNames)
        logging.debug("columnDescription:")
        logging.debug(columnDescription)
        logging.debug("dataManagerInfo:")
        logging.debug(dataManagerInfo)

        dataManagerInfo["NAME"]                                  = "TiledWgtSpectrumMartijn"
        dataManagerInfo["SPEC"]["DEFAULTTILESHAPE"]              = np.array([4, 40, 819], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["TileShape"] = np.array([4, 40, 819], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CubeShape"] = np.array([4, weights.shape[1], weights.shape[0]], dtype = np.int32)
        dataManagerInfo["SPEC"]["HYPERCUBES"]["*1"]["CellShape"] = np.array([4, weights.shape[1]], dtype = np.int32)

        logging.debug("dataManagerInfo (updated):")
        logging.debug(dataManagerInfo)

        logging.info("Removing column 'WEIGHT_SPECTRUM', if it exists...")
        if self.columnExists('WEIGHT_SPECTRUM'):
            self.t.removecols('WEIGHT_SPECTRUM')

        logging.info("Adding column 'WEIGHT_SPECTRUM'...")
        self.t.addcols(tables.makecoldesc('WEIGHT_SPECTRUM', columnDescription), dataManagerInfo)

        logging.info("Filling column 'WEIGHT_SPECTRUM'...")
        self.t.putcol('WEIGHT_SPECTRUM', weights)



def readArguments():
    import argparse
    parser=argparse.ArgumentParser("Adapt GMRT MS to LOFAR format.")
    parser.add_argument("-v", "--verbose", help="Be verbose. Default is False", required=False, action="store_true")
    parser.add_argument("ms_files", type=str, help="MeasurementSet name(s).", nargs="+")
    args=parser.parse_args()
    return vars(args)

if __name__=="__main__":
    start_time = time.time()

    args         = readArguments()
    verbose      = args["verbose"]
    ms_files     = args["ms_files"]

    if verbose: logging.basicConfig(level=logging.DEBUG)
    else: logging.basicConfig(level=logging.INFO)

    logging.info('Reading MSs...')
    MSs = []
    for ms_file in ms_files:
        MSs.append( MS(ms_file) )

    for MS in MSs:
        MS.removeColumns()
        MS.updatePolarisation()
        updateFreq = MS.updateFreqMetadata()
        MS.updateFieldMetadata()
        MS.updateIntervals()
        MS.updateColumns(updateFreq)

    logging.debug('Running time %.0f s' % (time.time()-start_time))
    logging.info('Done.')
