#!/usr/bin/python3
# make_catalogues.py *mosaic.fits [run into mosaics dir]

import bdsf
import os, glob, logging, argparse
from multiprocessing import Process, Queue
import lolss_config as cfg

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

restfrq = cfg.RESTFRQ
dir_mosaics = cfg.DIR_MOSAICS
dir_catalogues = cfg.DIR_CATALOGUES

parser = argparse.ArgumentParser(description='Make catalogues from mosaics')
parser.add_argument('-j', '--njobs', dest='njobs', type=int, default=8,
                    help='Number of parallel bdsf processes (default: 8)')
args = parser.parse_args()


def process_one(infile):
    catprefix = dir_catalogues + os.path.basename(infile).replace('.fits', '')

    if os.path.exists(catprefix + '.cat.fits'):
        logger.info('Skipping %s (catalogue already exists)', infile)
        return

    logger.info('Working on %s', infile)
    img = bdsf.process_image(infile, **cfg.bdsf_kwargs(frequency=restfrq))
    img.write_catalog(outfile=catprefix + '.cat.fits', catalog_type='srl', format='fits', correct_proj='True', clobber=True)
    img.write_catalog(outfile=catprefix + '.gaus.fits', catalog_type='gaul', format='fits', correct_proj='True', clobber=True)
    img.export_image(outfile=catprefix + '.rms.fits', img_type='rms', img_format='fits', clobber=True)
    img.export_image(outfile=catprefix + '.mean.fits', img_type='mean', img_format='fits', clobber=True)
    img.export_image(outfile=catprefix + '.resid.fits', img_type='gaus_resid', img_format='fits', clobber=True)
    img.export_image(outfile=catprefix + '.pybdsmmask.fits', img_type='island_mask', img_format='fits', clobber=True)
    img.write_catalog(outfile=catprefix + '.cat.reg', catalog_type='srl', format='ds9', correct_proj='True', clobber=True)


infiles = sorted(glob.glob(dir_mosaics + '*.fits'))
os.makedirs(dir_catalogues, exist_ok=True)
logger.info('All files: %s', infiles)

def _worker(q):
    while True:
        infile = q.get()
        if infile is None:
            break
        process_one(infile)

task_queue = Queue()
for infile in infiles:
    task_queue.put(infile)
for _ in range(args.njobs):
    task_queue.put(None)

workers = [Process(target=_worker, args=(task_queue,)) for _ in range(args.njobs)]
for w in workers:
    w.start()
for w in workers:
    w.join()
