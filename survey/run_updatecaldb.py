#!/usr/bin/env python3

import os, sys, glob
from LiLF.surveys_db import SurveysDB

dir_run = '/data/LBAsurvey'

def calibrator_tables_available(sdb, obsid):
    """
    check if calibrator data exist in the database
    """
    sdb.execute('SELECT * FROM observations WHERE id=%f' % obsid)
    r = sdb.cur.fetchall()
    if len(r) != 0 and r[0]['location'] != '': return True
    else: return False

# copy solutions in the repository
with SurveysDB(survey='lba',readonly=False) as sdb:
    for cal in glob.glob(dir_run+'/id*'):
        print('Working on %s...' % cal)
        obsid = int(cal.split('_-_')[0].split('/')[-1][2:])
        if not calibrator_tables_available(sdb, obsid):
            os.chdir(cal)
            # be sure is completed
            if not os.path.exists('cal-iono.h5'):
                print('Incomplete!')
                continue
    
            cal = os.path.basename(cal)
            print('Copy: cal*h5 -> herts:/beegfs/lofar/lba/calibration_solutions/%s' % cal)
            os.system('ssh herts "rm -rf /beegfs/lofar/lba/calibration_solutions/%s"' % cal)
            os.system('ssh herts "mkdir /beegfs/lofar/lba/calibration_solutions/%s"' % cal)
            os.system('scp -q cal-pa.h5 cal-amp.h5 cal-iono.h5 herts:/beegfs/lofar/lba/calibration_solutions/%s' % cal)
            os.system('scp -q -r plots* herts:/beegfs/lofar/lba/calibration_solutions/%s' % cal)
        
            # update the db
            sdb.execute('INSERT INTO observations (id,location,calibratordata) VALUES \
            (%i,"herts","%s")' % (obsid, "/beegfs/lofar/lba/calibration_solutions/"+cal))

