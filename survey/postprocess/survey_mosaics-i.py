#!/usr/bin/env python3

# For each HEALPix pixel defined in pointing_sets_file, collect the associated
# survey pointings and produce a Stokes-I mosaic using mosaic.py with beam
# correction. An empty FITS template is created first to define the output WCS
# and pixel grid; mosaics that are already complete (all pointings recorded in
# the FITS HISTORY header) are skipped.
#
# The template and the partial output live in <dir-mosaics>/tmp/ and the mosaic
# is moved to its final name only if mosaic.py exits successfully. A final
# HPxxxx-mosaicI.fits therefore always is a complete mosaic, and its HISTORY
# lists exactly the pointings that went into it.

import os, sys, logging, argparse
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from astropy import units as u
from astropy_healpix import HEALPix
import numpy as np
import lolss_config as cfg

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Mosaic survey pointings onto HEALPix tiles.')
parser.add_argument('--dir-mosaics',        default=cfg.DIR_MOSAICS,                    help='Output directory for mosaics')
parser.add_argument('--dir-done',           default=cfg.DIR_DONE,                       help='Directory containing completed pointings')
parser.add_argument('--file-image',         default=cfg.FILE_IMAGE,                                 help='MFS image filename within each pointing directory')
parser.add_argument('--file-beam',          default=cfg.FILE_BEAM,                                          help='Primary beam filename within each pointing directory')
parser.add_argument('--pointing-sets-file', default=cfg.FILE_POINTING_SETS, help='File mapping HEALPix pixels to pointing names')
parser.add_argument('--only-pointing',      default=None, type=int, nargs='+', metavar='HP',                    help='Only mosaic these HEALPix pixel IDs (default: all)')
parser.add_argument('--cellsize',           default=4, type=int,                                                 help='Output pixel size in arcsec (default: 4)')
parser.add_argument('--resolution',         default=15, type=int,                                                help='Output beam FWHM in arcsec (default: 15)')
parser.add_argument('--sbatch',             action='store_true',                                                 help='Submit mosaic jobs via sbatch instead of running locally')
parser.add_argument('-j', '--njobs',        default=1, type=int,                                                 help='Number of mosaics to make at the same time when running locally (default: 1)')
parser.add_argument('--retry',              action='store_true',                                                 help='With --sbatch: resubmit also the mosaics whose previous job never finished (killed or still running)')
args = parser.parse_args()

dir_mosaics        = args.dir_mosaics
dir_done           = args.dir_done
file_MFS_i         = args.file_image
file_beam          = args.file_beam
pointing_sets_file = args.pointing_sets_file
only_this_pointing = args.only_pointing
CELLSIZE           = args.cellsize
RESOLUTION         = args.resolution
use_sbatch         = args.sbatch

hp = HEALPix(nside=cfg.HEALPIX_NSIDE, frame='icrs')
logger.info(f'Total HEALPix pixels on sky: {hp.npix}')
area = hp.pixel_area.to(u.deg**2).value
logger.info(f'Area of one HEALPix pixel: {area:.4f} sq. deg')

singularity_img = '/homes/fdg/storage/pill-20250605.simg'
singularity_cmd = 'singularity exec --cleanenv --pwd /local/work/fdg --env PYTHONPATH=\$PYTHONPATH:/homes/fdg/storage/LiLF/:/homes/fdg/storage/scripts/,PATH=\$PATH:/homes/fdg/storage/LiLF/scripts/ --pid --writable-tmpfs -B/homes/fdg,/local/work/fdg,/iranet/groups/ulu/fdg '+singularity_img

class Pointing():

    def __init__(self, healpix, pointings_name_nearby):
        self.healpix = healpix
        self.pointings_name_nearby = pointings_name_nearby
        self.pointing_files = [dir_done+p+'/'+file_MFS_i for p in pointings_name_nearby]
        self.beam_files = [dir_done+p+'/'+file_beam for p in pointings_name_nearby]
        self.output_file = dir_mosaics+f'HP{healpix:04d}'+'-mosaicI.fits'
        self.log = dir_mosaics+'logs/'+f'HP{healpix:04d}'+'-mosaicI.log'
        self.batch_script = dir_mosaics+'batch/'+f'HP{healpix:04d}'+'-mosaicI.sh'
        # work files: never in dir_mosaics itself, so that every *.fits there is a finished mosaic
        self.template_file = dir_mosaics+'tmp/'+f'HP{healpix:04d}'+'-mosaicI.template.fits'
        self.partial_file = dir_mosaics+'tmp/'+f'HP{healpix:04d}'+'-mosaicI.partial.fits'
        self.failed_file = dir_mosaics+'tmp/'+f'HP{healpix:04d}'+'-mosaicI.failed'

    def missing_inputs(self):
        """Return the list of input images/beams that are not on disk."""
        return [f for f in self.pointing_files + self.beam_files if not os.path.exists(f)]

    def is_pending(self):
        """True if a previous sbatch job was submitted and neither finished nor reported a failure."""
        return os.path.exists(self.template_file) and not os.path.exists(self.failed_file)

    def mosaic(self):
        """
        Run (or submit) mosaic.py. The output is written to a partial file and renamed to the
        final mosaic only if mosaic.py succeeds.
        Returns True if the mosaic was made (local run) or submitted (sbatch), False on failure.
        """
        os.makedirs(os.path.dirname(self.log), exist_ok=True)
        for _stale in (self.partial_file, self.failed_file):
            if os.path.exists(_stale):
                os.remove(_stale)
        cmd = 'mosaic.py --images %s --beams %s --beamcorr --beamcut %s --find_noise --header %s --output %s' % \
                (' '.join(self.pointing_files), ' '.join(self.beam_files), cfg.BEAM_THRESHOLD, self.template_file, self.partial_file)
        finalize = f'mv {self.partial_file} {self.output_file} && rm -f {self.template_file}'
        if use_sbatch:
            os.makedirs(os.path.dirname(self.batch_script), exist_ok=True)
            cmd = singularity_cmd + ' ' + cmd.replace('mosaic.py', '/homes/fdg/storage/scripts/mosaic.py')
            logger.info(cmd)
            with open(self.batch_script, 'w') as _f:
                _f.write('#!/bin/bash\n')
                _f.write(f'#SBATCH --nodes=1\n')
                _f.write(f'#SBATCH --partition=lofar\n')
                _f.write(f'#SBATCH --ntasks-per-core=1\n')
                _f.write(f'#SBATCH --ntasks-per-node=1\n')
                _f.write(f'#SBATCH --cpus-per-task=36\n')
                _f.write(f'#SBATCH --time=240:00:00\n')
                _f.write(f'\nif {cmd} > {self.log} 2>&1; then\n')
                _f.write(f'    {finalize}\n')
                _f.write(f'else\n')
                _f.write(f'    rm -f {self.partial_file}\n')
                _f.write(f'    touch {self.failed_file}\n')
                _f.write(f'    exit 1\n')
                _f.write(f'fi\n')
            if os.system(f'sbatch {self.batch_script}') != 0:
                logger.error(f'sbatch submission failed for {self.batch_script}')
                os.remove(self.template_file)
                return False
            return True
        else:
            logger.info(cmd)
            if os.system(f'{cmd} > {self.log} 2>&1') != 0:
                logger.error(f'mosaic.py failed for HP{self.healpix:04d}, see {self.log}')
                if os.path.exists(self.partial_file):
                    os.remove(self.partial_file)
                return False
            if os.system(finalize) != 0:
                logger.error(f'Cannot move {self.partial_file} to {self.output_file}')
                return False
            return True

    def make_empty_mosaic(self):
        centre = hp.healpix_to_skycoord(self.healpix)
        corners = hp.boundaries_skycoord(self.healpix, step=1)
        dra, ddec = centre.spherical_offsets_to(corners)
        dra = dra.value
        ddec = ddec.value
        maxdist=np.max(np.abs([dra,ddec]))
        size = maxdist*1.15
        himsize = int(size/(CELLSIZE / 3600.0))
        logger.info(f'Creating empty template for {self.output_file} - size: {2*size:.2f} deg, himsize: {2*himsize} pixels')

        header=fits.Header()
        header['SIMPLE']=True
        header['BITPIX']=-32
        header['NAXIS']=2
        header['WCSAXES']=2
        header['NAXIS1']=2*himsize
        header['NAXIS2']=2*himsize
        header['CTYPE1']='RA---SIN'
        header['CTYPE2']='DEC--SIN'
        header['CUNIT1']='deg'
        header['CUNIT2']='deg'
        header['CRPIX1']=himsize
        header['CRPIX2']=himsize
        header['CRVAL1']=centre.ra.value
        header['CRVAL2']=centre.dec.value
        header['CDELT1']=-CELLSIZE / 3600.0
        header['CDELT2']=CELLSIZE / 3600.0
        header['RADESYS']='ICRS'
        header['EQUINOX']=2000.0
        header['LONPOLE']=180.0
        header['LATPOLE']=header['CRVAL2']
        header['BMAJ']=RESOLUTION / 3600.0
        header['BMIN']=RESOLUTION / 3600.0
        header['BPA']=0
        header['TELESCOP']='LOFAR'
        header['RESTFRQ']=cfg.RESTFRQ
        header['OBSERVER']='LoLSS'
        header['BUNIT']='JY/BEAM'
        header['BSCALE']=1.0
        header['BZERO']=0
        header['BTYPE']='Intensity'
        header['OBJECT']='HP-'+str(self.healpix)
        for p in self.pointings_name_nearby:
            header['HISTORY'] = 'POINTING ADDED: ' + p
        data = np.zeros((2*himsize,2*himsize), dtype=np.float32)
        os.makedirs(os.path.dirname(self.template_file), exist_ok=True)
        fits.writeto(self.template_file, header=header, data=data, overwrite=True)

    def check_todo(self):
        if not os.path.exists(self.output_file):
            return True
        # check if the mosaic is complete by looking at the header
        with fits.open(self.output_file) as _f:
            _hdr = _f[0].header
        _prefix = 'POINTING ADDED: '
        _in_mosaic = {str(_h)[len(_prefix):] for _h in _hdr.get('HISTORY', []) if str(_h).startswith(_prefix)}
        _missing = sorted(set(self.pointings_name_nearby) - _in_mosaic)
        _extra = sorted(_in_mosaic - set(self.pointings_name_nearby))
        if _missing:
            logger.info(f'Pointing(s) {_missing} not found in HISTORY of {self.output_file}, needs to be remade')
            return True
        if _extra:
            logger.info(f'Pointing(s) {_extra} in HISTORY of {self.output_file} are not in the pointing set anymore, needs to be remade')
            return True
        return False

# build grouped_pointings dict from pointing_sets.txt
# format: HP<healpix_pixel_id>  <pointing1> <pointing2> ...
grouped_pointings: dict[int, list[str]] = {}
with open(pointing_sets_file) as _f:
    for _line in _f:
        _line = _line.strip()
        if not _line or _line.startswith('#'):
            continue
        _parts = _line.split()
        grouped_pointings[int(_parts[0].lstrip('HP'))] = _parts[1:]
logger.info(f'Loaded {len(grouped_pointings)} coarse HEALPix pixels from {pointing_sets_file}')

# list of pointings to mosaic
failed = []
todo = []
for healpix, close_pointings in grouped_pointings.items():
    # skip if not in the list of pointings to mosaic
    if only_this_pointing is not None and healpix not in only_this_pointing:
        logger.info(f'Skipping pointing {healpix} (not in {only_this_pointing})')
        continue

    p = Pointing(healpix, close_pointings)
    if p.check_todo():
        logger.info('Pointing: %s %s', healpix, close_pointings)
        _missing_inputs = p.missing_inputs()
        if _missing_inputs:
            logger.error(f'HP{healpix:04d}: missing input file(s), not mosaicing: {_missing_inputs}')
            failed.append(healpix)
            continue
        if use_sbatch and p.is_pending() and not args.retry:
            logger.warning(f'HP{healpix:04d}: a previous job has not finished (found {p.template_file}); '
                           'skipping, use --retry to resubmit')
            continue
        todo.append(p)
    else:
        logger.info(f'Skipping pointing {healpix} (already complete)')


def run(p):
    """Make (or submit) one mosaic, return True on success."""
    try:
        p.make_empty_mosaic()
        return p.mosaic()
    except Exception:
        logger.exception(f'HP{p.healpix:04d}: unexpected error')
        return False


logger.info(f'{len(todo)} mosaic(s) to make')
if use_sbatch or args.njobs <= 1:
    results = [run(p) for p in todo]
else:
    # mosaic.py runs as a separate process: threads are enough to keep njobs of them going
    with ThreadPoolExecutor(max_workers=args.njobs) as pool:
        results = list(pool.map(run, todo))
failed += [p.healpix for p, ok in zip(todo, results) if not ok]

if failed:
    logger.error(f'{len(failed)} mosaic(s) failed: ' + ' '.join(f'HP{h:04d}' for h in failed))
    sys.exit(1)
