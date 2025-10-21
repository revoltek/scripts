#!/usr/bin/env python3

# run a series of scripts in pleiadi

import os, sys, socket
import string, random, glob, time
import logging
logging.basicConfig(level=logging.DEBUG)

singularity_img = f'/homes/fdg/storage/pill.simg'
singularity_cmd = f'singularity exec --cleanenv --pwd /local/work/fdg/surveyrun --env PYTHONPATH=\$PYTHONPATH:/homes/fdg/storage/LiLF/:/homes/fdg/storage/scripts/,PATH=\$PATH:/homes/fdg/storage/LiLF/scripts/ --pid --writable-tmpfs -B/homes/fdg,/local/work/fdg,/iranet/groups/ulu/fdg {singularity_img}'

dir_run = "/homes/fdg/data/surveyrun"
run_only = 1000 # limit run to this number of objects

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

# go in the run dir
for i in range(run_only):
    os.system(f'rm -r {dir_run} 2> /dev/null')
    os.system(f'rm -r /dev/shm/* 2> /dev/null')
    os.makedirs(dir_run)
    os.chdir(dir_run)
    logfile = "/homes/fdg/storage/run/sbatch_files/pill_"+id_generator()+"-"+socket.gethostname()+'.log'
    config_path = f"{dir_run}/lilf.config"
    with open(config_path, "w") as f:
        f.write("[PiLL]\n")
        f.write("minmaxhrs = 5,999\n")
        f.write(f"logfile = {logfile}\n")
        f.write("[LOFAR_timesplit]\n")
        f.write("ateam_clip = [CygA]\n")
    logging.info(f'Starting run num: {i}; log: {logfile}')
    os.system(f'{singularity_cmd} /homes/fdg/storage/LiLF/pipelines/PiLL-survey.py > {logfile} 2>&1')
    os.system(f'rm -r {dir_run} 2> /dev/null')
    os.system(f'rm -r /dev/shm/* 2> /dev/null')
