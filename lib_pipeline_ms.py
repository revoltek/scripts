#!/usr/bin/python

import os, sys
import logging
import numpy as np
import lofar.parmdb as parmdb

def merge_parmdb(parmdb_p, parmdb_a, parmdb_out, clobber=False):
    """
    Merges facet selfcal parmdbs into a parmdb for a single band

    Parameters
    ----------
    parmdb_p : str
        Filename of CommonScalarPhase and TEC parmdb
    parmdb_a : str
        Filename of Gain parmdb. The nearset match in frequency to that of the
        input band will be used
    parmdb_out : str
        Filename of output file
    clobber : bool, optional
        If True, overwrite existing output file

    """
    if type(clobber) is str:
        if clobber.lower() == 'true':
            clobber = True
        else:
            clobber = False

    if os.path.exists(parmdb_out):
        if clobber:
            shutil.rmtree(parmdb_out)
        else:
            return
    pdb_out = parmdb.parmdb(parmdb_out, create=True)

    # Copy over the CommonScalar phases and TEC
    pdb_p = parmdb.parmdb(parmdb_p)
    for parmname in pdb_p.getNames():
        parms = pdb_p.getValuesGrid(parmname)
        pdb_out.addValues(parms)

    # Copy over the Gains
    pdb_a = parmdb.parmdb(parmdb_a)
    for parmname in pdb_a.getNames():
        parms = pdb_a.getValuesGrid(parmname)
        pdb_out.addValues(parms)

    # Write values
    pdb_out.flush()


def find_nchan(ms):
    """
    Find number of channel in this ms
    """
    import pyrap.tables as tb
    t = tb.table(ms+'/SPECTRAL_WINDOW', ack=False)
    nchan = t.getcol('NUM_CHAN')
    t.close()
    assert (nchan[0] == nchan).all() # all spw have same channels?
    logging.debug('Channel in '+ms+': '+str(nchan[0]))
    return nchan[0]


def find_timeint(ms):
    """
    Get time interval in seconds
    """
    import pyrap.tables as tb
    t = tb.table(ms, ack=False)
    Ntimes = len(set(t.getcol('TIME')))
    t.close()
    t = tb.table(ms+'/OBSERVATION', ack=False)
    deltat = (t.getcol('TIME_RANGE')[0][1]-t.getcol('TIME_RANGE')[0][0])/Ntimes
    t.close()
    logging.debug('Time interval for '+ms+': '+str(deltat))
    return deltat


def get_phase_centre(ms):
    """
    Get the phase centre of the first source (is it a problem?) of an MS
    """
    import pyrap.tables as pt
    field_table = pt.table(ms + '/FIELD', ack=False)
    field_no = 0
    ant_no = 0
    direction = field_table.getcol('PHASE_DIR')
    ra = direction[ ant_no, field_no, 0 ]
    dec = direction[ ant_no, field_no, 1 ]
    return (ra*180/np.pi, dec*180/np.pi)
