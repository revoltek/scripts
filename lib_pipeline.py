#!/usr/bin/python

from threading import Thread
import subprocess
from Queue import Queue
import lofar.parmdb as parmdb
from scipy.interpolate import interp1d
import numpy as np
import logging, time

def thread_cmd(action_list, max_threads = 12):
    def worker(queue):
        for cmd in iter(queue.get, None):
            subprocess.call(cmd, shell=True)

    q = Queue()
    threads = [Thread(target=worker, args=(q,)) for _ in range(max_threads)]

    for i, t in enumerate(threads): # start workers
        t.daemon = True
        t.start()

    for action in action_list:
        q.put_nowait(action)
        # CASA and other tasks may conflict if initialized to close together
        time.sleep(1)
    for _ in threads: q.put(None) # signal no more commands
    for t in threads: t.join()

def add_coloring_to_emit_ansi(fn):
    # add methods we need to the class
    def new(*args):
        levelno = args[1].levelno
        if(levelno>=50):
            color = '\x1b[31m' # red
        elif(levelno>=40):
            color = '\x1b[31m' # red
        elif(levelno>=30):
            color = '\x1b[33m' # yellow
        elif(levelno>=20):
            color = '\x1b[32m' # green 
        elif(levelno>=10):
            color = '\x1b[35m' # pink
        else:
            color = '\x1b[0m' # normal
        args[1].msg = color + args[1].msg +  '\x1b[0m'  # normal
        #print "after"
        return fn(*args)
    return new

def set_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # get rid of all other loggers imported by modules
    for l in logger.handlers: l.setLevel('ERROR')
    logging.StreamHandler.emit = add_coloring_to_emit_ansi(logging.StreamHandler.emit)
    # create file handler which logs even debug messages
    check_rm('pipeline.logging')
    fh = logging.FileHandler('pipeline.logging')
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def merge_parmdb(parmdb_gain, parmdb_csp, parmdb_empty, parmdb_out):
    """
    parmdb_gain: slow varying gain solutions
    parmdb_csp: fast varying common scalar phase + TEC solutions
    parmdb_out: empty parmdb with gain and csp fast varying
    """
    
    pdb_gain = parmdb.parmdb(parmdb_gain)
    pdb_csp = parmdb.parmdb(parmdb_csp)
    pdb_empty = parmdb.parmdb(parmdb_empty)

    parms_gain = pdb_gain.getValuesGrid('*')
    parms_csp = pdb_csp.getValuesGrid('*')
    parms_empty = pdb_empty.getValuesGrid('*')

    # copy over the common scalar phases and TEC
    for key in parms_csp.keys():
        parms_empty[key]['values'] = parms_csp[key]['values']

    for key in parms_gain.keys():
        t_gain = parms_gain[key]['times']
        t_empty = parms_empty[key]['times']
        v_gain = parms_gain[key]['values'][:,0]

        # interpolate gain values
        parms_empty[key]['values'][:,0] = interp1d(t_gain, v_gain, kind='nearest', bounds_error=False)(t_empty)

        # put nearest values on the boundaries
        parms_empty[key]['values'][:,0][ np.where(t_empty<t_gain[0]) ] = v_gain[0]
        parms_empty[key]['values'][:,0][ np.where(t_empty>t_gain[-1]) ] = v_gain[-1]

    pdbnew = parmdb.parmdb(parmdb_out, create=True)
    pdbnew.addValues(parms_empty)
    pdbnew.flush()

#def check_rm(filename):
#    """
#    Check if file exists and remove it
#    Handle reg exp of glob
#    """
#    import os, errno, glob
#    for f in glob.glob(filename):
#        try:
#            os.system('rm -r '+f)
#        except OSError as e: # this would be "except OSError, e:" before Python 2.6
#            if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
#                raise # re-raise exception if a different error occured

def check_rm(regexp):
    """
    Check if file exists and remove it
    Handle reg exp of glob and spaces
    """
    import os, glob
    filenames = regexp.split(' ')
    for filename in filenames:
        # glob is used to check if file exists
        for f in glob.glob(filename):
            os.system('rm -r '+f)

