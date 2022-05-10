#!/usr/bin/env python3

#[structure]
#operation=STRUCTURE
#soltab=sol000/phaseOrig000
#refAnt=CS002LBA
#ant.regexp=CS*
#doUnwrap=True


import os,sys,glob,pickle,subprocess
import numpy as np
from losoto import h5parm
from LiLF import lib_ms

iono_h5parms = glob.glob('*3c*/cal-iono.h5')
#iono_h5parms = glob.glob('id790068_-_c14-o166.3_3c196/cal-iono.h5')
all_data = []
for iono_h5parm in iono_h5parms:
    print('Working on %s' % iono_h5parm)
    cmd = 'losoto %s losoto-structure.parset' % iono_h5parm
    output = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.decode('utf-8')
    beta = output.split('\n')[2].split('=')[-2].split(' ')[0]
    rdiff = output.split('\n')[2].split('=')[-1].split(' ')[0]
    h5 = h5parm.openSoltab(iono_h5parm,address='sol000/tec000')
    data = h5.getValues()
    data[1]['time'] = np.mean(data[1]['time'])
    caldir = iono_h5parm.split('/')[0]
    msname = sorted(glob.glob(caldir+'/*MS'))[0]
    ms = lib_ms.MS(msname)
    data[1]['elev'] = ms.elev
    data[1]['beta'] = beta
    data[1]['rdiff'] = rdiff
    all_data.append(data)

pickle.dump(all_data, open('iono.pickle','wb'))
