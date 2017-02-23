#!/usr/bin/python

import os, sys, time, pickle, random, shutil
import subprocess
import logging
from threading import Thread
from Queue import Queue
from scipy.interpolate import interp1d
import numpy as np
import lofar.parmdb as parmdb

from lib_pipeline_ms import *
from lib_pipeline_img import *

def get_cluster():
    """
    Find in which computing cluster the pipeline is running       
    """
    import socket
    hostname = socket.gethostname()
    if hostname == 'lgc1' or hostname == 'lgc2': return 'Hamburg'
    elif 'leidenuniv' in hostname: return 'Leiden'
    elif hostname[0:3] == 'lof': return 'CEP3'
    else: 
        logging.error('Hostname %s unknown.' % hostname)
        return 'Unknown'

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
        return fn(*args)
    return new


def set_logger(filename='pipeline.logging'):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # get rid of all other loggers imported by modules
    for l in logger.handlers: l.setLevel('CRITICAL')
    logging.StreamHandler.emit = add_coloring_to_emit_ansi(logging.StreamHandler.emit)
    # create file handler which logs even debug messages
    check_rm(filename)
    fh = logging.FileHandler(filename)
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


def size_from_facet(img, c_coord, pixsize):
    """
    Given an image, a new centre find the smallest image size which cover the whole image.
    img = CASA-image name
    c_coord = [ra,dec] in degrees, the wanted image centre
    pixsize = in arcsec, the final image will have this pixels size, so a rescaling might be needed
    """
    import pyrap.images
    img = pyrap.images.image(img)
    c = img.coordinates()
    # assumes images in a standard casa shape
    assert c.get_axes()[2] == ['Declination', 'Right Ascension']
    # assume same increment in image axes
    assert abs(c.get_increment()[2][0]) == abs(c.get_increment()[2][1])
    cen_y, cen_x = img.topixel([1,1,c_coord[1]*np.pi/180., c_coord[0]*np.pi/180.])[2:]
    max_y, max_x = img.shape()[2:]
    max_dist = max(max(cen_x, max_x - cen_x), max(cen_y, max_y - cen_y))
    max_dist = max_dist * abs(c.get_increment()[2][0])*180/np.pi*3600 / pixsize
    if max_dist > 6400: return 6400
    # multiply distance *2 (so to have the image size) and add 100% to be conservative
    max_dist = (max_dist*2)*2
    goodvalues = np.array([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512, 256, 128])
    shape = min(goodvalues[np.where(goodvalues>=max_dist)])
    del img
    return shape


class Scheduler():
    def __init__(self, qsub = None, max_threads = None, max_processors = None, log_dir = 'logs', dry = False):
        """
        qsub: if true call a shell script which call qsub and then wait 
        for the process to finish before returning
        max_threads: max number of parallel processes
        dry: don't schedule job
        max_processors: max number of processors in a node (ignored if qsub=False)
        """
        self.cluster = get_cluster()
        self.qsub = qsub
        # if qsub/max_thread/max_processors not set, guess from the cluster
        # if they are set, double check number are reasonable
        if self.qsub == None:
            if self.cluster == 'Hamburg': self.qsub = True
            else: self.qsub = False
        else:
            if (self.qsub == False and self.cluster == 'Hamburg') or \
               (self.qsub == True and (self.cluster == 'Leiden' or self.cluster == 'CEP3')):
                logging.critical('Qsub set to %s and cluster is %s.' % (str(qsub), self.cluster))
                sys.exit(1)

        if max_threads == None:
            if self.cluster == 'Hamburg': self.max_threads = 64
            elif self.cluster == 'Leiden': self.max_threads = 32
            elif self.cluster == 'CEP3': self.max_threads = 40
            else: self.max_threads = 12
        else:
            self.max_threads = max_threads

        if max_processors == None:
            if self.cluster == 'Hamburg': self.max_processors = 6
            elif self.cluster == 'Leiden': self.max_processors = 32
            elif self.cluster == 'CEP3': self.max_processors = 40
            else: self.max_processors = 12
        else:
            self.max_processors = max_processors

        self.dry = dry
        logging.info("Scheduler initialized for cluster "+self.cluster+" (Nproc: "+str(self.max_threads)+", multinode: "+str(self.qsub)+", max_processors: "+str(self.max_processors)+").")

        self.action_list = []
        self.log_list = [] # list of 2-lenght tuple of the type: (log filename, type of action)

        if not os.path.isdir(log_dir):
            logging.info('Creating log dir "'+log_dir+'".')
            os.makedirs(log_dir)
        self.log_dir = log_dir

    def add(self, cmd='', log='', log_append = False, cmd_type='', processors = None):
        """
        Add cmd to the scheduler list
        cmd: the command to run
        log: log file name that can be checked at the end
        log_append: if true append, otherwise replace
        cmd_type: can be a list of known command types as "BBS", "NDPPP"...
        processors: number of processors to use, can be "max" to automatically use max number of processors per node
        """
        if log != '': log = self.log_dir+'/'+log
        if log != '' and not log_append: cmd += ' > '+log+' 2>&1'
        if log != '' and log_append: cmd += ' >> '+log+' 2>&1'

        if processors != None and processors == 'max': processors = self.max_processors

        if self.qsub:
            # if number of processors not specified, try to find automatically
            if processors == None:
                processors = 1 # default use single CPU
                if "calibrate-stand-alone" == cmd[:21]: processors = 1
                if "NDPPP" == cmd[:5]: processors = 1
                if "wsclean" == cmd[:7]: processors = self.max_processors
                if "awimager" == cmd[:8]: processors = self.max_processors
            if processors > self.max_processors: processors = self.max_processors
            self.action_list.append([str(processors),'\''+cmd+'\''])
        else:
            self.action_list.append(cmd)

        if log != '':
            self.log_list.append((log,cmd_type))

    def add_casa(self, cmd = '', params = {}, wkd = None, log = '', log_append = False, processors = None):
        """
        Run a casa command pickling the parameters passed in params
        NOTE: running casa commands in parallel is a problem for the log file, better avoid
        alternatively all used MS and CASA must be in a separate working dir

        wkd = working dir (logs and pickle are in the pipeline dir)
        """

        if processors != None and processors == 'max': processors = self.max_processors
        if processors == None: processors=self.max_processors # default use entire node

        # since CASA can run in another dir, be sure log and pickle are in the pipeline working dir
        if log != '': log = os.getcwd()+'/'+self.log_dir+'/'+log
        pfile = os.getcwd()+'/casaparams_'+str(random.randint(0,1e9))+'.pickle'
        pickle.dump( params, open( pfile, "wb" ) )

        # exec in the script dir?
        if wkd == None: casacmd = 'casa --nogui --log2term --nologger -c '+cmd+' '+pfile
        elif os.path.isdir(wkd):
            casacmd = 'cd '+wkd+'; casa --nogui --log2term --nologger -c '+cmd+' '+pfile
        else:
            logging.error('Cannot find CASA working dir: '+wkd)
            sys.exit(1)

        if self.qsub:
            if log != '' and not log_append: casacmd = str(processors)+' \''+casacmd+' > '+log+' 2>&1'
            elif log != '' and log_append: casacmd = str(processors)+' \''+casacmd+' >> '+log+' 2>&1'
            else: casacmd = str(processors)+' \''+casacmd

            # clean up casa remnants in Hamburg cluster
            if self.cluster == 'Hamburg':
                self.action_list.append(casacmd+'; killall -9 -r dbus-daemon Xvfb python casa\*\'')
                if processors != self.max_processors:
                    logging.error('To clean annoying CASA remnants no more than 1 CASA per node is allowed.')
                    sys.exit(1)
            else:
                self.action_list.append(casacmd+'\'')
        else:
            if log != '' and not log_append: self.action_list.append(casacmd+' > '+log+' 2>&1')
            elif log != '' and log_append: self.action_list.append(casacmd+' >> '+log+' 2>&1')
            else: self.action_list.append(casacmd)

        if log != '':
            self.log_list.append((log,'CASA'))

    def run(self, check=False, max_threads=None):
        """
        If check=True then a check is done on every log in the log_list
        if max_thread != None, then it overrides the global values, useful for special commands that need a lower number of threads
        """
        def worker(queue):
            for cmd in iter(queue.get, None):
                if self.qsub and self.cluster == 'Hamburg':
                    # run in priority nodes
                    cmd = 'salloc --job-name LBApipe --reservation=important_science --time=24:00:00 --nodes=1 --tasks-per-node='+cmd[0]+\
                            ' /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env \''+cmd[1]+'\''
                    # run on all cluster
                    #cmd = 'salloc --job-name LBApipe --time=24:00:00 --nodes=1 --tasks-per-node='+cmd[0]+\
                    #        ' /usr/bin/srun --ntasks=1 --nodes=1 --preserve-env \''+cmd[1]+'\''
                #elif self.qsub and self.cluster == 'Leiden': cmd = 'qsub_waiter_lei.sh '+cmd[0]+' '+cmd[1]+' >> qsub.log'
                subprocess.call(cmd, shell=True)
    
        if max_threads == None: max_threads = self.max_threads
        q = Queue()
        threads = [Thread(target=worker, args=(q,)) for _ in range(max_threads)]
    
        for i, t in enumerate(threads): # start workers
            t.daemon = True
            t.start()
    
        for action in self.action_list:
            if self.dry: continue # don't schedule if dry run
            q.put_nowait(action)
        for _ in threads: q.put(None) # signal no more commands
        for t in threads: t.join()

        # check outcomes on logs
        if check:
            for log, cmd_type in self.log_list:
                self.check_run(log, cmd_type)

        # reset list of commands
        self.action_list = []
        self.log_list = []

    def check_run(self, log='', cmd_type=''):
        """
        Produce a warning if a command didn't close the log properly i.e. it crashed
        NOTE: grep, -L inverse match, -l return only filename
        """

        if not os.path.exists(log):
            logging.warning('No log file found to check results: '+log)
            return 1

        if cmd_type == 'BBS':
            out = subprocess.check_output('grep -L success '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('BBS run problem on:\n'+out)
                return 1

        elif cmd_type == 'NDPPP':
            out = subprocess.check_output('grep -L "Finishing processing" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -l "**** uncaught exception ****" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('NDPPP run problem on:\n'+out)
                return 1

        elif cmd_type == 'CASA':
            out = subprocess.check_output('grep -l "[a-z]Error" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -l "An error occurred running" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -l "\*\*\* Error \*\*\*" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('CASA run problem on:\n'+out)
                return 1

        elif cmd_type == 'wsclean':
            out = subprocess.check_output('grep -L "Cleaning up temporary files..." '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('WSClean run problem on:\n'+out)
                return 1

        elif cmd_type == 'python':
            out = subprocess.check_output('grep -l "[a-z]*Error" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            out += subprocess.check_output('grep -l "[a-z]*Critical" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('Python run problem on:\n'+out)
                return 1

        elif cmd_type == 'general':
            out = subprocess.check_output('grep -l -i "error" '+log+' ; exit 0', shell=True, stderr=subprocess.STDOUT)
            if out != '':
                logging.error('Run problem on:\n'+out)
                return 1

        else:
            logging.warning('Unknown command type for log checking: "'+cmd_type+'"')
            return 1

        return 0
