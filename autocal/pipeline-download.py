#!/usr/bin/python
# download from LTA using WGET

download_file = '3c295.txt'

###################################################

import sys, os, glob
import numpy as np
from lib_pipeline import *

set_logger()
check_rm('logs')
s = Scheduler(dry=False, max_threads = 5) # set here max number of threads here

df = open(download_file,'r')

logging.info('Downloading...')
for i, line in enumerate(df):
    s.add('wget -nv '+line[:-1]+' -O - | tar -x', log=str(i)+'.log', cmd_type='general')
    print 'wget -nv '+line[:-1]+' -O - | tar -x'
s.run(check=True)

logging.info("Done.")
