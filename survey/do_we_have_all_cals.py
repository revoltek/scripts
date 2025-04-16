#!/usr/bin/env python3

import os, sys, glob
from astropy.table import Table

# tgt-ids
grid = Table.read('/homes/fdg/storage/allsky-grid.fits')
tgt_ids = []
for field in grid['obsid']:
    for obs_id in field:
        if not obs_id in tgt_ids and obs_id != 0:
            tgt_ids.append(obs_id)

# cal-ids
cal_dirs = glob.glob('/homes/fdg/storage/surveycals/download/mss/id*')
print("There are %i cal dirs" % len(cal_dirs))
cal_ids = [int(cal_dir.split('id')[-1].split('_-_')[0]) for cal_dir in cal_dirs]
print("There are %i cal ids" % len(cal_ids))

missing_ids = []
for tgt_id in tgt_ids:
    if not tgt_id in cal_ids:
        missing_ids.append(tgt_id)

print('Missing these tgt-ids from calibrators:')
print(sorted(set(missing_ids)))

missing_ids = []
for cal_id in cal_ids:
    if not cal_id in tgt_ids:
        missing_ids.append(cal_id)

print('Missing these cal-ids from targets:')
print(sorted(set(missing_ids)))


