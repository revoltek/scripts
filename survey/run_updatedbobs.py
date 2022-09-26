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
    grid.convert_bytestring_to_unicode()
    
    with SurveysDB(survey='lba',readonly=False) as sdb:
        for field in grid:
            field_id = field['name']
            nobs = field['hrs']
            for obs_id, cycle, antset in zip(field['obsid'],field['cycle'],field['antset']):
                if obs_id != 0 and cycle != 'bad' and cycle != 'bug' and not obs_id in obs_to_skip:
                    # select sparse
                    if 'Sparse' in antset:
                        print('Add to the db: %i -> %s (%s)' % (obs_id, field_id+'s', cycle))
                        sdb.execute('INSERT INTO field_obs (obs_id,field_id) VALUES (%i,"%s")' % (obs_id, field_id+'s'))
                    # select outer
                    elif 'Outer' in antset:
                        print('Add to the db: %i -> %s (%s)' % (obs_id, field_id+'o', cycle))
                        sdb.execute('INSERT INTO field_obs (obs_id,field_id) VALUES (%i,"%s")' % (obs_id, field_id+'o'))
            nobs_s = np.sum(field['antset'][(field['cycle'] != 'bad') & (field['cycle'] != 'bug')] == 'LBA Sparse Even')
            nobs_o = np.sum(field['antset'][(field['cycle'] != 'bad') & (field['cycle'] != 'bug')] == 'LBA Outer')
            if nobs_s >= 3 or nobs_o >= 3:
                if nobs > 7:
                    priority = 3
                else:
                    if field['distAteam'] < 15: priority = 2
                    else: priority = 1
                if nobs_s >= 3:
                    print("%s: set as observed (%i hr - priority: %i)" % (field_id+'s', nobs_s, priority))
                    sdb.execute('UPDATE fields SET status="Observed" WHERE id="%s"' % (field_id+'s'))
                    sdb.execute('UPDATE fields SET priority=%i WHERE id="%s"' % (priority,field_id+'s'))
                if nobs_o >= 3:
                    print("%s: set as observed (%i hr - priority: %i)" % (field_id+'o', nobs_o, priority))
                    sdb.execute('UPDATE fields SET status="Observed" WHERE id="%s"' % (field_id+'o'))
                    sdb.execute('UPDATE fields SET priority=%i WHERE id="%s"' % (priority,field_id+'o'))
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
    sdb.execute('SELECT id,status,clustername,nodename,noise,nvss_ratio FROM fields WHERE status!="Observed" and status!="Not started"')
    r = sdb.cur.fetchall()
    for i, entry in enumerate(r):
        if entry['status'] == 'Done':
            print('%03i) ID: %s (%s - %s: %s) Noise: %.2f mJy, NVSSratio: %.2f' \
                    % (i, entry['id'], entry['status'], entry['clustername'], entry['nodename'],entry['noise']*1e3,entry['nvss_ratio']))
        else:
            print('%03i) ID: %s (%s - %s: %s)' % (i, entry['id'], entry['status'], entry['clustername'], entry['nodename']))
