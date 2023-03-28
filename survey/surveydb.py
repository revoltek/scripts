#!/usr/bin/env python3
# This script uses the data from LBA (in the gridfile) and populate the field and field_obs table
# At the same time set to "Observed" all fields that have at least 3 observed hours

import os, sys, argparse, re
from LiLF.surveys_db import SurveysDB
from astropy.table import Table
import numpy as np

gridfile = 'allsky-grid.fits'

parser = argparse.ArgumentParser(description='Stage and download MS from the LOFAR LTA.')
parser.add_argument('--gridfile', '-g', dest="gridfile", help='The gridfile as created with update_allsky-grid.py', default=gridfile)
parser.add_argument('--updatedb', '-u', action="store_true", help='Update the databse using the gridfile.')
parser.add_argument('--reset', '-r', dest="reset", help='If "all" reset the db to "Not started" for all fields. If a field is specified it only reset it to "Observed".', default=None)
parser.add_argument('--incompletereset', '-i', action="store_true", help='Reset the fields that are not "Done"/"Not started" to "Observed".')
parser.add_argument('--sethighpriority', '-p', dest="sethighpriority", help='Give the pointing name so to set its priority to 0 (maximum).', default=None)
parser.add_argument('--caldir', '-c', dest="caldir", help='Update calibrators if given. It needs cal dir, e.g. "/homes/fdg/storage/surveycals/done"', default=None)
parser.add_argument('--show', '-s', dest="show", help="If 'done' shows completed runs; if 'running' shows ongoing/failed runs; if 'all' shows all, incuding runs that have not started yet.")
args = parser.parse_args()

if args.sethighpriority is not None:
    with SurveysDB(survey='lba',readonly=False) as sdb:  
        print("INFO: set %s to maximum priority" % args.sethighpriority)
        sdb.execute('UPDATE fields SET priority=%i WHERE id="%s"' % (0,args.sethighpriority))
        sys.exit()

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
    with SurveysDB(survey='lba',readonly=True) as sdb:
        sdb.execute('SELECT obs_id FROM field_obs')
        obs_to_skip = [ x['obs_id'] for x in sdb.cur.fetchall() ]
    print('The following obs are already in the DB:', obs_to_skip)
    
    grid = Table.read('allsky-grid.fits')
    grid.convert_bytestring_to_unicode()
    
    with SurveysDB(survey='lba',readonly=False) as sdb:
        for field in grid:
            field_id = field['name']
            nobs = field['hrs']

            # first fill the observations db with good/bad obs
            for obs_id, status in zip(field['obsid'],field['cycle']):
                if obs_id == 0: continue
                if status != 'bad' and status != 'bug': status = 'good'
                sdb.execute('SELECT status FROM observations WHERE id="%i"' % obs_id)
                status_old = sdb.cur.fetchall()
                print('Add to the observations db: %i (%s - was: %s)' % (obs_id, status, status_old))

                if status == 'bad' and len(status_old) == 0:
                    # track bad obs
                    sdb.execute('INSERT INTO observations (id, status) VALUES (%i, "bad")' % obs_id)
                elif status == 'good' and len(status_old) == 0:
                    # add good obs
                    sdb.execute('INSERT INTO observations (id, status) VALUES (%i, "good")' % obs_id)
                elif status != status_old:
                    sdb.execute('UPDATE observations SET status="%s" WHERE id="%i"' % (status, obs_id))
                else:
                    pass # already correct in the db

            # now fille the filed_obs considering only good data
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
                    sdb.execute('UPDATE fields SET status="Observed", priority=%i WHERE id="%s"' % (priority,field_id+'s'))
                if nobs_o >= 3:
                    print("%s: set as observed (%i hr - priority: %i)" % (field_id+'o', nobs_o, priority))
                    sdb.execute('UPDATE fields SET status="Observed", priority=%i WHERE id="%s"' % (priority,field_id+'o'))
    sys.exit()

if args.caldir is not None:
    cals = glob.glob(args.caldir+'/id*')
    # copy solutions in the repository
    with SurveysDB(survey='lba',readonly=False) as sdb:
        for cal in cals:
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
                sdb.execute('UPDATE observations SET location="herts", calibratordata="%s" \
                            WHERE id=%i' % ("/beegfs/lofar/lba/calibration_solutions/"+cal, obsid))
    sys.exit()

# default: show the db
if args.show is not None:
    with SurveysDB(survey='lba',readonly=True) as sdb:
        if args.show == 'all':
            sdb.execute('SELECT id,status,priority FROM fields WHERE status="Observed" order by priority desc')
            r = sdb.cur.fetchall()
            sdb.execute('SELECT field_id FROM field_obs')
            all_fields = [x['field_id'] for x in sdb.cur.fetchall()]
            for i, entry in enumerate(r):
                hrs = sum(np.array(all_fields) == entry['id'])
                print('%03i) ID: %s - %i hrs (%s - priority: %i)' % (i, entry['id'], hrs, entry['status'], entry['priority']))
            print("############################")
        elif args.show == 'all' or args.show == 'done':
            sdb.execute('SELECT id,status,clustername,nodename,noise,nvss_ratio,nvss_match,flag_frac FROM fields WHERE status="Done"')
            r = sdb.cur.fetchall()
            for i, entry in enumerate(r):
                print('%03i) ID: %s (%s - %s: %s) Noise: %.2f mJy, NVSSratio: %.2f (matches: %i) - flags: %.1f%%' \
                            % (i, entry['id'], entry['status'], entry['clustername'], entry['nodename'],entry['noise']*1e3,entry['nvss_ratio'],entry['nvss_match'],entry['flag_frac']*100))
        elif args.show == 'all' or args.show == 'running':
            sdb.execute('SELECT id,status,clustername,nodename,noise,nvss_ratio,nvss_match,flag_frac FROM fields WHERE status!="Observed" and status!="Not started" and status!="Done"')
            r = sdb.cur.fetchall()
            for i, entry in enumerate(r):
                print('%03i) ID: %s (%s - %s: %s)' % (i, entry['id'], entry['status'], entry['clustername'], entry['nodename']))
        else:
            print('With "show" use: all, running, or done.')
    sys.exit()

print('No action selected...')
