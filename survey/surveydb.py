#!/usr/bin/env python3
# This script uses the data from LBA (in the gridfile) and populate the field and field_obs table

import sys, argparse, glob
from LiLF.surveys_db import SurveysDB
from astropy.table import Table
import numpy as np
import gspread
from oauth2client.service_account import ServiceAccountCredentials

keyjson = '/homes/fdg/.ssh/lolss-459012-2b1a95f95b82.json'
gridfile = 'allsky-grid.fits'
caldirroot = ('/iranet/groups/ulu/fdg/surveycals/done/')
tgtdirroot = ('/iranet/groups/ulu/fdg/surveytgts/download*/mss/')

parser = argparse.ArgumentParser(description='Stage and download MS from the LOFAR LTA.')
parser.add_argument('--gridfile', '-g', dest="gridfile", help='The gridfile as created with update_allsky-grid.py', default=gridfile)
parser.add_argument('--updatedb', '-u', action="store_true", help='Update the databse using the gridfile.')
parser.add_argument('--reset', '-r', dest="reset", help='If "all" reset the db to "Not observed" for all observed fields. If a field is specified it only reset it to "Observed" (or the value of --status).', default=None)
parser.add_argument('--status', '-t', dest="status", help='To use with --reset to define the finale status, default: "Observed".', default="Observed")
parser.add_argument('--incompletereset', '-i', action="store_true", help='Reset the fields that are not "Done"/"Observed"/"Not observed" to "Downloaded".')
parser.add_argument('--sethighpriority', '-p', dest="sethighpriority", help='Give the pointing name so to set its priority to 0 (maximum).', default=None)
parser.add_argument('--show', '-s', dest="show", help="If 'done' shows completed runs; if 'running' shows ongoing/failed runs; if 'all' shows all, incuding runs that have Observed yet.")
parser.add_argument('--google', '-o', dest="google", help="Update google spreadsheet with the current status of the fields.", action="store_true")
args = parser.parse_args()

if args.google:
    print("Udating google spreadsheet...")
    scope = ["https://spreadsheets.google.com/feeds", "https://www.googleapis.com/auth/drive"]
    creds = ServiceAccountCredentials.from_json_keyfile_name(keyjson, scope)
    client = gspread.authorize(creds)
    sheet = client.open_by_url("https://docs.google.com/spreadsheets/d/1nPjEC78_XOr40FVmTvQ14Pey03su2RQpdOIO6yZOV68/edit?usp=sharing").sheet1  # First worksheet
    id_sheet = sheet.col_values(1) # this columns contains the id of the fields
    # Store all updates in memory
    updates = {}

    with SurveysDB(survey='lba',readonly=True) as sdb:
        sdb.execute('SELECT id,status,clustername,nodename,noise,nvss_ratio,nvss_match,flag_frac,end_date FROM fields')
        r = sdb.cur.fetchall()
        for i, entry in enumerate(r):
            if entry['id'] in id_sheet:
                if entry['noise'] is None: entry['noise'] = 0
                if entry['flag_frac'] is None: entry['flag_frac'] = 0
                if entry['nvss_ratio'] is None: entry['nvss_ratio'] = 0
                # be sure to not list older data
                if entry['status'] == 'Downloaded' or entry['status'] == 'Observed':
                    entry['noise'] = 0
                    entry['flag_frac'] = 0
                    entry['nvss_ratio'] = 0
                    entry['end_date'] = None
                if entry['status'] != 'Done' and entry['status'] != 'Downloaded' and entry['status'] != 'Observed':
                    entry['status'] += ' ('+str(entry['nodename'])+')'

                row_index = id_sheet.index(entry['id']) + 1  # +1 because gspread is 1-indexed
                updates[row_index] = [
                None,  # Column A (ID) — leave unchanged
                None,  # Column B (notes) - not used
                None,  # Column C (ra)
                None,  # Column D (dec)
                entry['status'],                     # Column E
                entry['noise'] * 1e3,                # Column F
                entry['nvss_ratio'],                 # Column G
                entry['flag_frac'] * 100,            # Column H
                str(entry['end_date']),              # Column I
                None, None, None, None, None, None, None, None, None, None, None, None # Columns J to U - leave unchanged
                ]
    
    min_row = min(updates.keys())
    max_row = max(updates.keys())
    data_block = []

    for i in range(min_row, max_row + 1):
        if i in updates:
            data_block.append(updates[i])
        else:
            data_block.append([None] * 21)  # Preserve row alignment

    # Update in one batch call (columns E to I = cols 5–9)
    range_start = f"E{min_row}"
    range_end = f"I{max_row}"
    update_range = f"{range_start}:{range_end}"
    sheet.update(update_range, [row[4:9] for row in data_block])

    print("✅ Spreadsheet updated successfully.")

    sys.exit()            

if args.sethighpriority is not None:
    with SurveysDB(survey='lba',readonly=False) as sdb:  
        print("INFO: set %s to maximum priority" % args.sethighpriority)
        sdb.execute('UPDATE fields SET priority=%i WHERE id="%s"' % (0,args.sethighpriority))
        sys.exit()

if args.reset is not None:
    with SurveysDB(survey='lba',readonly=False) as sdb:  
        if args.reset == 'all':
            print(f"WARNING: RESET ALL POINTINGS to \"Not observed\"")
            input("Press Enter to continue...")
            sdb.execute(f'UPDATE fields SET status="Not observed"')
            sdb.execute('DELETE from field_obs')
        else:
            print("WARNING: reset pointing %s to \"%s\"" % (args.reset, args.status))
            input("Press Enter to continue...")
            sdb.execute('UPDATE fields SET status="%s" where id="%s"' % (args.status, args.reset))
    sys.exit()

if args.incompletereset:
    with SurveysDB(survey='lba',readonly=False) as sdb:  
        print("WARNING: RESET INCOMPLETE POINTINGS to \"Downloaded\"")
        input("Press Enter to continue...")
        sdb.execute('UPDATE fields SET status="Downloaded" where status!="Done" and status!="Observed" and status!="Not observed"')
        sys.exit()

if args.updatedb:
    with SurveysDB(survey='lba',readonly=True) as sdb:
        sdb.execute('SELECT obs_id FROM field_obs')
        obs_to_skip = [ x['obs_id'] for x in sdb.cur.fetchall() ]
    print('The following obs are already in the DB:', obs_to_skip)
    
    grid = Table.read('allsky-grid.fits')
    grid.convert_bytestring_to_unicode()
    grid = grid[np.random.permutation(len(grid))] # shuffle the table


    with SurveysDB(survey='lba',readonly=False) as sdb:
        all_fields = {}
        for field in grid:
            field_id = field['name']
            nobs = field['hrs']

            # first fill the "observations" table with good/bad obs
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
                    continue # already correct in the db, skip re-adding the calibrator

                # and add calibrator data location
                caldir = glob.glob(f'{caldirroot}/id{obs_id}*')
                if len(caldir) == 0:
                    print(f'WARNING: missing cal data for obs_id: {obs_id}')
                else:
                    sdb.execute('UPDATE observations SET calibratordata="%s" WHERE id="%i"' % (caldir, obs_id))

            # now fill the "filed_obs" table considering only good data and separating outer and sparse
            for obs_id, cycle, antset in zip(field['obsid'],field['cycle'],field['antset']):
                if obs_id != 0 and cycle != 'bad' and cycle != 'bug' and not obs_id in obs_to_skip:
                    # select sparse
                    if 'Sparse' in antset:
                        field_id_antset = field_id+'s'
                    # select outer
                    elif 'Outer' in antset:
                        field_id_antset = field_id+'o'

                    print('Add to the field_obs db: %i -> %s (%s)' % (obs_id, field_id_antset, cycle))
                    sdb.execute('INSERT INTO field_obs (obs_id,field_id) VALUES (%i,"%s")' % (obs_id, field_id_antset))
                    
                    if not field_id_antset in all_fields.keys(): all_fields[field_id_antset] = [obs_id] 
                    else: all_fields[field_id_antset].append(obs_id)

        # finally fill the "fields" table
        for field_id in all_fields.keys():
            # check what obs_id is ready
            missing=0
            for obs_id in all_fields[field_id]:
                #print(f"checking: {tgtdirroot}/id{obs_id}_-_{field_id[:-1]}")
                if len(glob.glob(f'{tgtdirroot}/id{obs_id}_-_{field_id[:-1]}')) == 0:
                    missing+=1
                    print('WARNING: %s: missing %i' % (field_id, obs_id))

            hrs = len(all_fields[field_id])
            if missing > 0:
                print('WARNING: %s: %i hrs missing out of %i' % (field_id, missing, hrs))

            hrs -= missing
            if hrs > 0:
                if missing == 0 and (grid[grid['name'] == field_id[:-1]]['distAteam'] > 15): priority = 3
                elif missing > 0 and (grid[grid['name'] == field_id[:-1]]['distAteam'] > 15): priority = 2
                else: priority = 1

                print("%s: set as Downloaded (%i hr - priority: %i)" % (field_id, hrs, priority))
                sdb.execute('UPDATE fields SET status="Downloaded", priority=%i WHERE id="%s"' % (priority,field_id))
            else:
                print("%s: set as Observed (%i hr - priority: 0)" % (field_id, hrs))
                sdb.execute('UPDATE fields SET status="Observed", priority=0 WHERE id="%s"' % (field_id))

    sys.exit()

# default: show the db
if args.show is not None:
    with SurveysDB(survey='lba',readonly=True) as sdb:
        if args.show == 'all':
            print("Observed / Downloaded:")
            sdb.execute('SELECT id,status,priority FROM fields WHERE status="Observed" or status="Downloaded" order by priority desc')
            r = sdb.cur.fetchall()
            #sdb.execute('SELECT field_id FROM field_obs')
            sdb.execute(f'SELECT field_id FROM field_obs t1 JOIN observations t2 ON t1.obs_id=t2.id WHERE t2.status="good"')
            r2 = sdb.cur.fetchall()
            #print(r2)
            all_fields = [x['field_id'] for x in r2]
            for i, entry in enumerate(r):
                # TODO: only select good
                hrs = sum(np.array(all_fields) == entry['id'])
                print('%i) ID: %s - %i hrs (%s) - priority: %s)' \
                        % (i, entry['id'], hrs, entry['status'], entry['priority']))
            print("############################")

        if args.show == 'all' or args.show == 'done':
            print("Done:")
            sdb.execute('SELECT id,status,clustername,nodename,noise,nvss_ratio,nvss_match,flag_frac,end_date FROM fields WHERE status="Done" order by end_date asc')
            r = sdb.cur.fetchall()
            for i, entry in enumerate(r):
                print('%03i) ID: %s (%s - %s: %s; %s) Noise: %.2f mJy, NVSSratio: %.2f (matches: %i) - flags: %.1f%%' \
                            % (i, entry['id'], entry['status'], entry['clustername'], entry['nodename'],entry['end_date'],entry['noise']*1e3,entry['nvss_ratio'],entry['nvss_match'],entry['flag_frac']*100))
            print("############################")      

        if args.show == 'all' or args.show == 'running':
            print("Running:")
            sdb.execute('SELECT id,status,clustername,nodename FROM fields WHERE status!="Not Observed" and status!="Downloaded" and status!="Observed" and status!="Done"')
            r = sdb.cur.fetchall()
            for i, entry in enumerate(r):
                print('%03i) ID: %s (%s - %s: %s)' % (i, entry['id'], entry['status'], entry['clustername'], entry['nodename']))
            print("############################")

        if not args.show in ['all', 'running', 'done']:
            print('With "show" use: all, running, or done.')
    sys.exit()

print('No action selected...')
