#!/usr/bin/env python3

# run a series of scripts in pleiadi

import os, sys
import string, random, glob, time
import logging
logging.basicConfig(level=logging.DEBUG)

singularity_img = '/homes/fdg/storage/pill.simg'
singularity_cmd = 'singularity exec --cleanenv --pwd /local/work/fdg --env PYTHONPATH=\$PYTHONPATH:/homes/fdg/storage/LiLF/:/homes/fdg/storage/scripts/,PATH=\$PATH:/homes/fdg/storage/LiLF/scripts/ --pid --writable-tmpfs -B/homes/fdg,/local/work/fdg,/iranet/groups/ulu/fdg '+singularity_img

dir_run = "/homes/fdg/data/surveyrun"
run_only = 1 # limit run to this number of objects

# go in the run dir
os.makedirs(dir_run)
os.chdir(dir_run)

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
for i in range(run_only):

