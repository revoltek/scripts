#!/usr/bin/env python3

# run a series of scripts in pleiadi

import os, sys
import string, random, glob, time
import logging
logging.basicConfig(level=logging.DEBUG)

singularity_img = '/homes/fdg/storage/pill-20220720.simg'
singularity_cmd = 'singularity run --pid --writable-tmpfs --containall --cleanenv -B/homes/fdg:/homes/fdg,/local/work/fdg:/local/work/fdg,/iranet/lofarfs2/lofar2/fdg:/iranet/lofarfs2/lofar2/fdg '+singularity_img

dir_storage_cals = '/iranet/lofarfs2/lofar2/fdg/surveycals'
dir_storage_tgts = '/iranet/lofarfs2/lofar2/fdg/surveytgts'

dir_run = "/homes/fdg/storage/run"

# go in the run dir
os.chdir(dir_run)
os.system('rm sbatch_files/*')

class Schedule():

    def __init__(self, name):

        def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
            return ''.join(random.choice(chars) for _ in range(size))

        self.name = name
        self.file_sbatch = "sbatch_files/"+name+'_'+id_generator()

    def submit(self):
        logging.info('Scheduling: %s' % self.file_sbatch)
        os.system("sbatch -p lofar %s" % (self.file_sbatch))

    def prepare_sbatch(self, mode):
    
        if mode == 'cal':
            dir_orig = dir_storage_cals+'/download/mss/'+self.name
            dir_dest = dir_storage_cals+'/done/'+self.name

            content = """#!/bin/bash
                         #SBATCH --nodes=1
                         ### number of hyperthreading threads
                         #SBATCH --ntasks-per-core=1
                         ### number of MPI tasks per node
                         #SBATCH --ntasks-per-node=1
                         ### number of openmp threads
                         #SBATCH --cpus-per-task=36
                         #SBATCH --time=4:00:00
                         #SBATCH -o %s
                         rm -r /local/work/fdg/*
                         mkdir /local/work/fdg/cal
                         echo -e "[LOFAR_cal]\ndata_dir=%s\n" > /local/work/fdg/lilf.config
                         singularity exec --cleanenv --pwd /local/work/fdg/cal --env PYTHONPATH=\$PYTHONPATH:/homes/fdg/storage/LiLF/:/homes/fdg/storage/scripts/,PATH=\$PATH:/homes/fdg/storage/LiLF/scripts/ --pid --writable-tmpfs --containall -B/homes/fdg,/local/work/fdg,/iranet/lofarfs2/lofar2/fdg %s /homes/fdg/storage/LiLF/pipelines/LOFAR_cal.py
                         mv /local/work/fdg/cal/plots* /local/work/fdg/cal/cal*h5 /local/work/fdg/cal/*logger /local/work/fdg/cal/logs %s
                         rm -r /local/work/fdg/*
                         """ % (self.file_sbatch+".log", dir_orig, singularity_img, dir_dest)
            content = ''.join(line.lstrip(' \t') for line in content.splitlines(True)) # remove spaces

        elif mode == 'self':
            pass

        elif mode == 'dd':
            pass

        with open(self.file_sbatch, 'w') as f:
            f.write(content)


# cals
all_cals = sorted(glob.glob(dir_storage_cals+'/download/mss/id*'))
logging.info('Setting up %i jobs.' % len(all_cals))

for i, dir_orig in enumerate(all_cals):
    dir_dest = dir_orig.replace('download/mss','done')
    if not os.path.exists(dir_dest):
        os.makedirs(dir_dest)

    # skip if already done
    if len(glob.glob(dir_dest+'/*h5')) == 4:
        continue

    # skip if not present
    if len(glob.glob(dir_orig+'/*MS')) == 0:
        continue

    c = Schedule(name=dir_orig.split('/')[-1])
    c.prepare_sbatch('cal')
    c.submit()

    # separate initial calls so initial cp is diluted
    if i < 50:
        time.sleep(60)

# tgts self

# tgts dd-cal
