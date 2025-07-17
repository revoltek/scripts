#!/usr/bin/env python3

# run a series of scripts in pleiadi

import os, sys
import string, random, glob, time
import logging
logging.basicConfig(level=logging.DEBUG)

singularity_img = '/homes/fdg/storage/pill.simg'
singularity_cmd = 'singularity exec --cleanenv --pwd /local/work/fdg --env PYTHONPATH=\$PYTHONPATH:/homes/fdg/storage/LiLF/:/homes/fdg/storage/scripts/,PATH=\$PATH:/homes/fdg/storage/LiLF/scripts/ --pid --writable-tmpfs -B/homes/fdg,/local/work/fdg,/iranet/groups/ulu/fdg '+singularity_img

dir_storage_cals = '/iranet/groups/ulu/fdg/surveycals'
dir_storage_tgts = '/iranet/groups/ulu/fdg/surveytgts'

dir_run = "/homes/fdg/storage/run"
run_only = 500 # limit run to this number of objects

# go in the run dir
os.chdir(dir_run)
#os.system('rm sbatch_files/*')

class Scheduler():

    def __init__(self, name):

        def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
            return ''.join(random.choice(chars) for _ in range(size))

        self.name = name
        self.file_sbatch = "sbatch_files/"+name+'_'+id_generator()
        self.file_log = self.file_sbatch

    def submit(self):
        logging.info('Scheduling: %s' % self.file_sbatch)
        os.system("sbatch -p lofar %s" % (self.file_sbatch))

    def prepare_sbatch(self, mode):
    
        if mode == 'cal':
            dir_orig = dir_storage_cals+'/download/mss/'+self.name
            dir_dest = dir_storage_cals+'/done/'+self.name

            content = f"""#!/bin/bash
                         #SBATCH --nodes=1
                         #SBATCH --partition=lofar
                         ### number of hyperthreading threads
                         #SBATCH --ntasks-per-core=1
                         ### number of MPI tasks per node
                         #SBATCH --ntasks-per-node=1
                         ### number of openmp threads
                         #SBATCH --cpus-per-task=36
                         #SBATCH --time=10:00:00
                         #SBATCH -o {self.file_log}-%N.log
                         #SBATCH --job-name={self.name}
                         rm -r /local/work/fdg/*
                         mkdir -p /local/work/fdg/data-bkp
                         cp {dir_orig}/* /local/work/fdg/data-bkp
                         {singularity_cmd} /homes/fdg/storage/LiLF/pipelines/LOFAR_cal.py
                         mv /local/work/fdg/cal-pa.h5 /local/work/fdg/cal-bp.h5 /local/work/fdg/cal-fr.h5 /local/work/fdg/cal-iono-cs.h5 /local/work/fdg/cal-iono.h5 {dir_dest}
                         mv /local/work/fdg/plots* {dir_dest}
                         mv /local/work/fdg/pipeline-cal_*logger /local/work/fdg/logs_pipeline-cal_* /local/work/fdg/pipeline-cal.walker {dir_dest}
                         rm -r /local/work/fdg/*
                         """
            content = ''.join(line.lstrip(' \t') for line in content.splitlines(True)) # remove spaces

        elif mode == 'pill':

            content = f"""#!/bin/bash
                         #SBATCH --nodes=1
                         ### number of hyperthreading threads
                         #SBATCH --partition=lofar
                         #SBATCH --ntasks-per-core=1
                         ### number of MPI tasks per node
                         #SBATCH --ntasks-per-node=1
                         ### number of openmp threads
                         #SBATCH --cpus-per-task=36
                         #SBATCH --time=240:00:00
                         #SBATCH -o {self.file_log}-%N.log
                         #SBATCH --job-name={self.name}
                         rm -r /local/work/fdg/*
                         rm -r /dev/shm/*
                         mkdir -p /local/work/fdg/
                         echo "[PiLL]" > /local/work/fdg/lilf.config
                         echo "minmaxhrs = 1,3" >> /local/work/fdg/lilf.config
                         echo "[LOFAR_timesplit]" >> /local/work/fdg/lilf.config
                         echo "ateam_clip = [CygA]" >> /local/work/fdg/lilf.config
                         {singularity_cmd} /homes/fdg/storage/LiLF/pipelines/PiLL-survey.py
                         rm -r /local/work/fdg/*
                         rm -r /dev/shm/*
                         """
            content = ''.join(line.lstrip(' \t') for line in content.splitlines(True)) # remove spaces

 
        with open(self.file_sbatch, 'w') as f:
            f.write(content)

###
# cals
#all_cals = sorted(glob.glob(dir_storage_cals+'/download/mss/id*'))
#logging.info('Setting up %i jobs.' % len(all_cals))
#
#i=0
#for dir_orig in all_cals:
#    dir_dest = dir_orig.replace('download/mss','done')
#
#    if not os.path.exists(dir_dest):
#        os.makedirs(dir_dest)
#
#    # skip if already done
#    if len(glob.glob(dir_dest+'/*h5')) == 5 or len(glob.glob(dir_dest+'/*h5')) == 0:
#    #if len(glob.glob(dir_dest+'/*logger')) == 1:
#        continue
#
#    # skip if not present
#    if len(glob.glob(dir_orig+'/*MS')) == 0 and len(glob.glob(dir_orig+'/*tar')) == 0:
#        continue
#
#    c = Scheduler(name=dir_orig.split('/')[-1])
#    c.prepare_sbatch('cal')
#    c.submit()
#
#    # separate initial calls so initial cp is diluted
#    if i < 34:
#        time.sleep(120)
#    i+=1
#    
#    if i == run_only:
#        sys.exit()

###
# tgts
for i in range(run_only):
    c = Scheduler('pill')
    c.prepare_sbatch('pill')
    c.submit()

    # separate initial calls so the stagings+downloads are diluted
    if i < 36:
        time.sleep(600)
