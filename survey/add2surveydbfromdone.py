#!/usr/bin/env python3

import os, sys, glob, pickle, datetime
from LiLF.surveys_db import SurveysDB

done = done = [d for d in glob.glob('/homes/fdg/storage/surveytgts/done/*') if os.path.isdir(d) and not (d.endswith('error') or d.endswith('tests'))]

with SurveysDB(survey='lba',readonly=False) as sdb:
    for d in done:
        target = os.path.basename(d)
        sdb.execute('SELECT status FROM fields WHERE id="%s"' % target)
        r = sdb.cur.fetchall()
        if r[0] == 'Done':
            print('Target %s already in survey db as Done.' % target)
            continue

        print('Adding %s to survey db as Done.' % target)
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        sdb.execute('UPDATE fields SET status="%s", end_date="%s" WHERE id="%s"' % ('Done',timestamp, target))
        sdb.execute('UPDATE fields SET username="%s", clustername="%s", nodename="%s" WHERE id="%s"' % \
                    ('fdg', 'pleiadi', 'unknown', target))
        with open(d+'/quality.pickle', 'rb') as f:
            qdict = pickle.load(f)
            sdb.execute('UPDATE fields SET noise="%s", nvss_ratio="%s", nvss_match="%s", flag_frac="%s" WHERE id="%s"' \
                % (qdict['ddserial_c0_rms'],qdict['nvss_ratio'], qdict['nvss_match'], qdict['flag_frac'],  target))
            

