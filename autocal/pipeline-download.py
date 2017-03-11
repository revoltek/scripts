#!/usr/bin/python
# download from LTA using WGET

rename = False
fix_tables = True

###################################################

import sys, os, re, glob, time
import numpy as np
import pyrap.tables as pt
from autocal.lib_pipeline import *
from astropy.time import Time

if os.path.exists('html.txt'):
    download_file = 'html.txt'
else:
    download_file = None # just renaming

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

mss = glob.glob('*MS')

if fix_tables:

    logging.info('Fix MS table...')
    for ms in mss:
        os.system('fixMS_TabRef.py '+ms)

    # only ms created in range (2/2013->2/2014)
    with pt.table(mss[0]+'/OBSERVATION', readonly=True, ack=False) as obs:
        t = Time(obs.getcell('TIME_RANGE',0)[0]/(24*3600.), format='mjd')
        time = np.int(t.iso.replace('-','')[0:8])
    if time > 20130200 and time < 20140300:
        logging.info('Fix beam table...')
        for ms in mss:
            s.add('/home/fdg/scripts/fixinfo/fixbeaminfo '+ms, log=ms+'_fixbeam.log')
        s.run(check=False)


# rename files using codename and freq
if rename:
    logging.info('Renaming...')
    flog = open('renamed.txt', 'a', 0)
    start = time.time()
    for ms in mss:

        with pt.table(ms+'/FIELD', readonly=True, ack=False) as t:
            code = t.getcell('CODE',0)
        if code == '':
            with pt.table(ms+'/OBSERVATION', readonly=True, ack=False) as t:
                code = t.getcell('LOFAR_TARGET',0)[0]
        
        code = code.lower()

        # get freq
        with pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False) as t:
            freq = t.getcell('REF_FREQUENCY',0)

        # get time (saved in ms as MJD in seconds)
        with pt.table(ms+'/OBSERVATION', readonly=True, ack=False) as t:
            time = Time(t.getcell('TIME_RANGE',0)[0]/(24*3600.), format='mjd')
            time = time.iso.replace('-','').replace(' ','').replace(':','')[0:12]

        # make name
        pattern = re.compile("^c[0-9][0-9]-.*$")
        # is survey?
        if pattern.match(code):
            cycle_obs, sou = code.split('_')
            if not os.path.exists(cycle_obs+'/'+sou): os.makedirs(cycle_obs+'/'+sou)
            newName = cycle_obs+'/'+sou+'/'+sou+'_t'+time+'_SB'+str(nu2num(freq/1.e6))+'.MS'
        else:
            newName = code+'_t'+time+'_SB'+str(nu2num(freq/1.e6))+'.MS'
            if newName == ms: continue
    
        logging.debug('Rename '+ms+' -> '+newName)
        os.system('mv '+ms+' '+newName)

        flog.write(ms+'\n')
    flog.close()

logging.info("Done.")
