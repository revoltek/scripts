#!/usr/bin/python3
# make_catalogues.py *mosaic.fits [run into mosaics dir]

import bdsf
import os, glob, logging, argparse
from multiprocessing import Process, Queue

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

restfrq = 54e6
dir_mosaics = "/homes/fdg/storage/surveytgts/mosaics/healpix/"
dir_catalogues = "/homes/fdg/storage/surveytgts/catalogues/"

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
    img = bdsf.process_image(infile, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=50, rms_box_bright=(50,10), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, atrous_do=False, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq)
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
