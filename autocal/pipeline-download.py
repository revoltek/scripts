#!/usr/bin/python
# download from LTA using WGET

download_file = 'html.txt'
#download_file = None # just renaming
rename = True

###################################################

import sys, os, re, glob, time
import numpy as np
import pyrap.tables as pt
from autocal.lib_pipeline import *
from astropy.time import Time

def nu2num(nu):
    """
    Get the SB number from the freq
    """
    nu_clk = 200. # 160 or 200 MHz, clock freq
    n = 1 # nyquist zone (1 for LBA, 2 for HBA low, 3 for HBA mid-high)

    if nu_clk == 200:
        SBband = 195312.5/1e6
    elif nu_clk == 160:
        SBband = 156250.0/1e6

    return np.int(np.floor((1024./nu_clk) * (nu - (n-1) * nu_clk/2.)))


set_logger('pipeline-download.logging')
check_rm('logs')
s = Scheduler(dry=False, max_threads = 4) # set here max number of threads here

if not download_file is None:
    df = open(download_file,'r')

    logging.info('Downloading...')
    downloaded = glob.glob('*MS')
    # add renamed files
    if os.path.exists('renamed.txt'):
        with open('renamed.txt','r') as flog:
            downloaded += [line.rstrip('\n') for line in flog]

    for i, line in enumerate(df):
        ms = re.findall(r'L[0-9]*_SB[0-9]*_uv\.dppp\.MS', line)[0]
        if ms in downloaded: continue
        s.add('wget -nv "'+line[:-1]+'" -O - | tar -x', log=str(i)+'.log', cmd_type='general')
    #    print 'wget -nv "'+line[:-1]+'" -O - | tar -x'
        logging.debug('Queue download of: '+ms)
    s.run(check=True)

# rename files using codename and freq
if rename:
    logging.info('Renaming...')
    flog = open('renamed.txt', 'a', 0)
    regex = re.compile(r"^L[0-9]*_")
    regex2 = re.compile(r"_uv\.dppp")
    regex3 = re.compile(r"_SB[0-9]*[_uv]*.MS")
    start = time.time()
    for ms in glob.glob('*MS'):
    
        field = pt.table(ms+'/FIELD', readonly=True, ack=False)
        code = field.getcell('CODE',0)
        field.close()
        newName = regex.sub(code+'_', ms)
        newName = regex2.sub('', newName)

        cycle_obs, sou = code.split('_')
        if not os.path.exists(cycle_obs+'/'+sou): os.makedirs(cycle_obs+'/'+sou)

        # get freq
        spw = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
        freq = spw.getcell('REF_FREQUENCY',0)
        spw.close()

        # get time (saved in ms as MJD in seconds)
        obs = pt.table(ms+'/OBSERVATION', readonly=True, ack=False)
        t = Time(obs.getcell('TIME_RANGE',0)[0]/(24*3600.), format='mjd')
        time = t.iso.replace(':','')[11:15]
        obs.close()

        newName = regex3.sub('_'+time+'_SB'+str(nu2num(freq/1.e6))+'.MS', newName)

        logging.debug('Rename '+ms+' -> '+cycle_obs+'/'+sou+'/'+newName)
        os.system('mv '+ms+' '+cycle_obs+'/'+sou+'/'+newName)
        os.system('fixMS_TabRef.py '+cycle_obs+'/'+sou+'/'+newName)
        flog.write(ms+'\n')
    flog.close()

logging.info("Done.")
