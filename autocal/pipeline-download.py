#!/usr/bin/python
# download from LTA using WGET

download_file = 'cyg1.txt'

###################################################

import sys, os, re, glob
import numpy as np
from lib_pipeline import *

set_logger()
check_rm('logs')
s = Scheduler(dry=False, max_threads = 5) # set here max number of threads here

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

logging.info("Done.")
