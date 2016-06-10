#!/usr/bin/python
# download from LTA using WGET

download_file = 'html.txt'
rename = False

###################################################

import sys, os, re, glob, time
import numpy as np
import pyrap.tables as pt
from lib_pipeline import *

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


set_logger()
check_rm('logs')
s = Scheduler(dry=False, max_threads = 2) # set here max number of threads here

df = open(download_file,'r')

logging.info('Downloading...')
downloaded = glob.glob('*MS')
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
    regex = re.compile(r"^L[0-9]*_")
    regex2 = re.compile(r"_uv\.dppp")
    regex3 = re.compile(r"_SB[0-9]*.MS")
    start = time.time()
    for ms in glob.glob('*MS'):
    
        field = pt.table(ms+'/FIELD', readonly=True, ack=False)
        code = field.getcell('CODE',0)
        field.close()
        newName = regex.sub(code+'_', ms)
        newName = regex2.sub('', newName)

        sou = code.split('_')[1]
        if not os.path.exists(sou): os.makedirs(sou)

        spw = pt.table(ms+'/SPECTRAL_WINDOW', readonly=True, ack=False)
        freq = spw.getcell('REF_FREQUENCY',0)
        spw.close()
        newName = regex3.sub('_SB'+str(nu2num(freq/1.e6))+'.MS', newName)

        logging.debug('Rename '+ms+' -> '+sou+'/'+newName)
        os.system('mv '+ms+' '+sou+'/'+newName)

logging.info("Done.")
