#!/usr/bin/env python3
# This script query the LTA and populate the field and field_obs table
# At the same time set to "Observed" all fields that have at least 3 observed hours

import os, sys, argparse, re
from LiLF.surveys_db import SurveysDB
from astropy.table import Table
import numpy as np

gridfile = 'allsky-grid.fits'

parser = argparse.ArgumentParser(description='Stage and download MS from the LOFAR LTA.')
parser.add_argument('--gridfile', '-g', dest="gridfile", help='The gridfile as created with update_allsky-grid.py', default=gridfile)
parser.add_argument('--skip', '-s', action="store_true", help='Skip observations already present in field_obs, \
        this is faster but might miss some target to update as "Observed" in the field table.')
parser.add_argument('--updatedb', '-u', action="store_true", help='Update the databse using the gridfile.')
parser.add_argument('--reset', '-r', dest="reset", help='If "all" reset the db to "Not started" for all fields. If a field is specified it only reset it to "Observed".', default=None)
parser.add_argument('--incompletereset', '-i', action="store_true", help='Reset the fields that are not "Done"/"Not started" to "Observed".')
args = parser.parse_args()

if args.reset is not None:
    with SurveysDB(survey='lba',readonly=False) as sdb:  
        if args.reset == 'all':
            print("WARNING: RESET ALL POINTINGS to \"Not started\"")
            input("Press Enter to continue...")
            sdb.execute('UPDATE fields SET status="Not started"')
            sdb.execute('DELETE from field_obs')
        else:
            print("WARNING: reset pointing %s to \"Observed\"" % args.reset)
            input("Press Enter to continue...")
            sdb.execute('UPDATE fields SET status="Observed" where id="%s"' % args.reset)
        sys.exit()

if args.incompletereset:
    with SurveysDB(survey='lba',readonly=False) as sdb:  
        print("WARNING: RESET INCOMPLETE POINTINGS to \"Observed\"")
        input("Press Enter to continue...")
        sdb.execute('UPDATE fields SET status="Observed" where status!="Done" and status!="Not started"')
        sys.exit()

if args.updatedb:
    skip_obs = args.skip

    # get obs_id already done
    with SurveysDB(survey='lba',readonly=True) as sdb:
        sdb.execute('select obs_id from field_obs')
        obs_to_skip = [ x['obs_id'] for x in sdb.cur.fetchall() ]
    print('The following obs are already in the DB:', obs_to_skip)
    
    grid = Table.read('allsky-grid.fits')
    
    with SurveysDB(survey='lba',readonly=False) as sdb:
        for field in grid:
            field_id = field['name']
            nobs = field['hrs']
            for obs_id, cycle in zip(field['obsid'],field['cycle']):
                if obs_id != 0 and cycle != b'bad' and cycle != b'bug' and not obs_id in obs_to_skip:
                    print('Add to the db: %i -> %s (%s)' % (obs_id, field_id, cycle))
                    sdb.execute('INSERT INTO field_obs (obs_id,field_id) VALUES (%i,"%s")' % (obs_id, field_id))
            if nobs >= 3:
                print("%s: set as observed (%i)" % (field_id, nobs))
                sdb.execute('UPDATE fields SET status="Observed" WHERE id="%s"' % (field_id))
                if nobs > 7:
                    sdb.execute('UPDATE fields SET priority=3 WHERE id="%s"' % (field_id))
                else:
                    if field['distAteam'] < 15:
                        sdb.execute('UPDATE fields SET priority=2 WHERE id="%s"' % (field_id))
                    else:
                        sdb.execute('UPDATE fields SET priority=1 WHERE id="%s"' % (field_id))
    sys.exit()

# default: show the db
with SurveysDB(survey='lba',readonly=True) as sdb:
    sdb.execute('SELECT id,status,priority FROM fields WHERE status="Observed" order by priority desc')
    r = sdb.cur.fetchall()
    sdb.execute('SELECT field_id FROM field_obs')
    all_fields = [x['field_id'] for x in sdb.cur.fetchall()]
    for i, entry in enumerate(r):
        hrs = sum(np.array(all_fields) == entry['id'])
        print('%03i) ID: %s - %i hrs (%s - priority: %i)' % (i, entry['id'], hrs, entry['status'], entry['priority']))
    print("############################")
    sdb.execute('SELECT id,status,clustername,nodename FROM fields WHERE status!="Observed" and status!="Not started"')
    r = sdb.cur.fetchall()
    for i, entry in enumerate(r):
        print('%03i) ID: %s (%s - %s: %s)' % (i, entry['id'], entry['status'], entry['clustername'], entry['nodename']))
