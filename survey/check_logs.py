#!/usr/bin/env python3

import os, sys, glob

cal_dirs = sorted(glob.glob('id*'))

for cal_dir in cal_dirs:
    try:
        logger = sorted(glob.glob(cal_dir+'/*logger'))[-1]
    except:
        print('%s: no logger' % cal_dir)
        continue

    with open(logger, "r") as f:
        last_line = f.readlines()[-1]
        if not "Done" in last_line:
            print("%s: incomplete calibration - %s " % (cal_dir, last_line[:-1]))

    concatlog = glob.glob(cal_dir+'/logs_pipeline-cal*/concat.log')[0]
    with open(concatlog, "r") as f:
        for ln in f:
            if "ntimes" in ln:
                ntimes = int(ln.split(":")[-1].replace(' ',''))
                if ntimes != 899:
                    print('%s: strange ntimes: %i' % (cal_dir, ntimes))
                break

