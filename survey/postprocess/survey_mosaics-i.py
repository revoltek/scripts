#!/usr/bin/env python3

# For each HEALPix pixel defined in pointing_sets_file, collect the associated
# survey pointings and produce a Stokes-I mosaic using mosaic.py with beam
# correction. An empty FITS template is created first to define the output WCS
# and pixel grid; mosaics that are already complete (all pointings recorded in
# the FITS HISTORY header) are skipped.
#
# The template and the partial output live in <dir-mosaics>/tmp/ and the mosaic
# is moved to its final name only if mosaic.py exits successfully. A final
# HPxxxx-<suffix>.fits therefore always is a complete mosaic, and its HISTORY
# lists exactly the pointings that went into it.
#
# The images are not given to mosaic.py directly: first they are copied into
# <dir-tmp>/HPxxxx-<suffix>/, and the copies are corrected with the results of
# survey_astrometry.py (astrometry_pointings.fits):
#   1. CRVAL1/CRVAL2 are shifted by the global offset of the pointing (dCRVAL1/2),
#      unless the offset is not reliable (pointing flagged FEW or SCATTER);
#   2. the facets flagged with --blank-flags (default: OFFSET|SCATTER) are set to nan;
#   3. the pointings flagged SCATTER (facets shifted in different directions) are not used
#      at all, unless --keep-scatter.
# The images in dir-done are only read. A pointing that needs no correction is
# symlinked instead of copied. The work directory is removed when mosaic.py ends.

import os, sys, shutil, logging, argparse
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u
from astropy_healpix import HEALPix
import numpy as np
import lolss_config as cfg
from survey_astrometry import read_facets, FLAG_FEW, FLAG_OFFSET, FLAG_SCATTER, PFLAG_BAD, PFLAG_SCATTER

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
parser.add_argument('--suffix',             default='mosaicI',                                                   help='Suffix of the output names: HPxxxx-<suffix>.fits (default: mosaicI)')
parser.add_argument('--resolution',         default=15, type=int,                                                help='Output beam FWHM in arcsec (default: 15)')
parser.add_argument('--sbatch',             action='store_true',                                                 help='Submit mosaic jobs via sbatch instead of running locally')
parser.add_argument('-j', '--njobs',        default=1, type=int,                                                 help='Number of mosaics to make at the same time when running locally (default: 1)')
parser.add_argument('--retry',              action='store_true',                                                 help='With --sbatch: resubmit also the mosaics whose previous job never finished (killed or still running)')
parser.add_argument('--dir-tmp',            default=None,                                                        help='Where the corrected copies of the images are made (default: <dir-mosaics>/tmp/)')
parser.add_argument('--astrometry',         default=cfg.DIR_ANALYSIS+'astrometry/astrometry_pointings.fits',     help='Output of survey_astrometry.py (default: %(default)s)')
parser.add_argument('--no-astrometry',      action='store_true',                                                 help='Do not correct the astrometry nor blank facets')
parser.add_argument('--no-shift',           action='store_true',                                                 help='Do not shift the pointings (facets are still blanked)')
parser.add_argument('--blank-flags',        default=FLAG_OFFSET|FLAG_SCATTER, type=int,                          help='Blank the facets with any of these flag bits (2 OFFSET, 4 SCATTER; 0: none; default: %(default)s)')
parser.add_argument('--keep-scatter',       action='store_true',                                                 help='Use also the pointings flagged SCATTER by survey_astrometry.py (default: leave them out)')
parser.add_argument('--file-facets',        default='facetsS-c0.reg',                                            help='Facet regions within each pointing directory (default: %(default)s)')
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
dir_tmp            = args.dir_tmp if args.dir_tmp else dir_mosaics+'tmp/'


class Astrometry():
    """Corrections of each pointing from astrometry_pointings.fits."""

    def __init__(self, path, blank_flags, shift=True, exclude_scatter=True):
        self.path = os.path.abspath(path)
        self.settings = f'shift={"yes" if shift else "no"} blank_flags={blank_flags} exclude_scatter={"yes" if exclude_scatter else "no"}'
        with fits.open(path) as _h:
            ptab = Table(_h['POINTINGS'].data)
            ftab = Table(_h['FACETS'].data)
        # global shift (deg, to add to CRVAL1/2); none if the offset is not reliable
        self.shift = {}
        for r in ptab:
            if shift and not r["Flag"] & PFLAG_BAD:
                self.shift[str(r['Pointing']).strip()] = (float(r['dCRVAL1']), float(r['dCRVAL2']))
        self.known = {str(p).strip() for p in ptab['Pointing']}
        # pointings left out of the mosaics
        self.excluded = {str(r['Pointing']).strip() for r in ptab if r['Flag'] & PFLAG_SCATTER} if exclude_scatter else set()
        # facets to blank: {pointing: [(index in the region file, name, flag)]}
        self.bad = {}
        if blank_flags:
            sel = (ftab['Facet'] >= 0) & ((ftab['Flag'] & blank_flags) != 0) & ((ftab['Flag'] & FLAG_FEW) == 0)
            for r in ftab[sel]:
                self.bad.setdefault(str(r['Pointing']).strip(), []).append((int(r['Facet']), str(r['Name']).strip(), int(r['Flag'])))
        logger.info(f'Astrometry from {self.path}: {len(self.known)} pointings, {len(self.excluded)} excluded, {len(self.shift)} shifted, '
                    f'{sum(len(v) for v in self.bad.values())} facets to blank in {len(self.bad)} pointings')

    def history(self):
        """Line written in the HISTORY of the mosaics, to know which correction they have."""
        return f'ASTROMETRY: {self.path} {self.settings}'

    def needs_copy(self, pointing):
        return pointing in self.shift or pointing in self.bad

    def correct(self, pointing, src, dst, regfile):
        """Write in dst the image src with the shift and the blanked facets of pointing."""
        with fits.open(src, memmap=False) as _h:
            hdr = _h[0].header
            data = _h[0].data
        if pointing in self.bad:
            # facet polygons in pixels of the original (not shifted) image
            facets = read_facets(regfile, WCS(hdr).celestial)
            ny, nx = data.shape[-2:]
            for idx, name, flag in self.bad[pointing]:
                if idx >= len(facets) or facets[idx][0] != name:
                    raise RuntimeError(f'{pointing}: facet {idx} {name} of {self.path} not found in {regfile}')
                path = facets[idx][1]
                x0, y0 = np.floor(path.vertices.min(axis=0)).astype(int)
                x1, y1 = np.ceil(path.vertices.max(axis=0)).astype(int)
                x0, y0, x1, y1 = max(x0, 0), max(y0, 0), min(x1+1, nx), min(y1+1, ny)
                if x0 >= x1 or y0 >= y1: continue
                yy, xx = np.mgrid[y0:y1, x0:x1]
                inside = path.contains_points(np.column_stack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
                data[..., y0:y1, x0:x1][..., inside] = np.nan
                hdr['HISTORY'] = f'Facet {name} (flag {flag}) blanked by survey_mosaics-i.py'
                logger.debug(f'{pointing}: blanked facet {name} ({inside.sum()} pixels)')
        if pointing in self.shift:
            dcrval1, dcrval2 = self.shift[pointing]
            hdr['CRVAL1'] += dcrval1
            hdr['CRVAL2'] += dcrval2
            hdr['HISTORY'] = f'Shifted by dRA {dcrval1*3600*np.cos(np.radians(hdr["CRVAL2"])):.2f}" dDEC {dcrval2*3600:.2f}" by survey_mosaics-i.py'
        fits.writeto(dst, data=data, header=hdr, overwrite=True)


astrometry = None
if not args.no_astrometry:
    if not os.path.exists(args.astrometry):
        logger.error(f'{args.astrometry} not found: run survey_astrometry.py, or use --no-astrometry')
        sys.exit(1)
    astrometry = Astrometry(args.astrometry, args.blank_flags, shift=not args.no_shift, exclude_scatter=not args.keep_scatter)

hp = HEALPix(nside=cfg.HEALPIX_NSIDE, frame='icrs')
logger.info(f'Total HEALPix pixels on sky: {hp.npix}')
area = hp.pixel_area.to(u.deg**2).value
logger.info(f'Area of one HEALPix pixel: {area:.4f} sq. deg')

singularity_img = '/homes/fdg/storage/pill-20250605.simg'
singularity_cmd = r'singularity exec --cleanenv --pwd /local/work/fdg --env PYTHONPATH=\$PYTHONPATH:/homes/fdg/storage/LiLF/:/homes/fdg/storage/scripts/,PATH=\$PATH:/homes/fdg/storage/LiLF/scripts/ --pid --writable-tmpfs -B/homes/fdg,/local/work/fdg,/iranet/groups/ulu/fdg '+singularity_img

class Pointing():

    def __init__(self, healpix, pointings_name_nearby):
        self.healpix = healpix
        self.pointings_name_nearby = pointings_name_nearby
        self.pointing_files = [dir_done+p+'/'+file_MFS_i for p in pointings_name_nearby]
        self.beam_files = [dir_done+p+'/'+file_beam for p in pointings_name_nearby]
        self.output_file = dir_mosaics+f'HP{healpix:04d}-{args.suffix}.fits'
        self.log = dir_mosaics+'logs/'+f'HP{healpix:04d}-{args.suffix}.log'
        self.batch_script = dir_mosaics+'batch/'+f'HP{healpix:04d}-{args.suffix}.sh'
        # work files: never in dir_mosaics itself, so that every *.fits there is a finished mosaic
        self.template_file = dir_mosaics+'tmp/'+f'HP{healpix:04d}-{args.suffix}.template.fits'
        self.partial_file = dir_mosaics+'tmp/'+f'HP{healpix:04d}-{args.suffix}.partial.fits'
        self.failed_file = dir_mosaics+'tmp/'+f'HP{healpix:04d}-{args.suffix}.failed'
        # corrected copies of the images, given to mosaic.py
        self.work_dir = dir_tmp+f'HP{healpix:04d}-{args.suffix}/'
        self.work_files = [self.work_dir+p+'-'+file_MFS_i for p in pointings_name_nearby]
        self.facet_files = [dir_done+p+'/'+args.file_facets for p in pointings_name_nearby]

    def missing_inputs(self):
        """Return the list of input images/beams (and needed facet files) that are not on disk."""
        files = self.pointing_files + self.beam_files
        if astrometry is not None:
            files += [f for p, f in zip(self.pointings_name_nearby, self.facet_files) if p in astrometry.bad]
        return [f for f in files if not os.path.exists(f)]

    def prepare_inputs(self):
        """Copy (with the astrometric corrections) or link the images in the work directory."""
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)
        os.makedirs(self.work_dir)
        for p, src, dst, reg in zip(self.pointings_name_nearby, self.pointing_files, self.work_files, self.facet_files):
            if astrometry is not None and p not in astrometry.known:
                logger.warning(f'HP{self.healpix:04d}: {p} not in {astrometry.path}, used without corrections')
            if astrometry is not None and astrometry.needs_copy(p):
                astrometry.correct(p, src, dst, reg)
            else:
                os.symlink(os.path.abspath(src), dst)

    def clean(self):
        shutil.rmtree(self.work_dir, ignore_errors=True)

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
                (' '.join(self.work_files), ' '.join(self.beam_files), cfg.BEAM_THRESHOLD, self.template_file, self.partial_file)
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
                _f.write(f'    rm -rf {self.work_dir}\n')
                _f.write(f'    {finalize}\n')
                _f.write(f'else\n')
                _f.write(f'    rm -rf {self.work_dir}\n')
                _f.write(f'    rm -f {self.partial_file}\n')
                _f.write(f'    touch {self.failed_file}\n')
                _f.write(f'    exit 1\n')
                _f.write(f'fi\n')
            if os.system(f'sbatch {self.batch_script}') != 0:
                logger.error(f'sbatch submission failed for {self.batch_script}')
                os.remove(self.template_file)
                self.clean()
                return False
            return True
        else:
            logger.info(cmd)
            _ok = os.system(f'{cmd} > {self.log} 2>&1') == 0
            self.clean()
            if not _ok:
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
        header['HISTORY'] = astrometry.history() if astrometry is not None else 'ASTROMETRY: none'
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
        # HISTORY cards are split at 72 characters: compare the joined text
        _expected = astrometry.history() if astrometry is not None else 'ASTROMETRY: none'
        if _expected.replace(' ', '') not in ''.join(str(_h) for _h in _hdr.get('HISTORY', [])).replace(' ', ''):
            logger.info(f'{self.output_file} was made with a different astrometric correction, needs to be remade')
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

    if astrometry is not None:
        _excluded = [c for c in close_pointings if c in astrometry.excluded]
        if _excluded:
            logger.warning(f'HP{healpix:04d}: pointing(s) {_excluded} flagged SCATTER, left out')
            close_pointings = [c for c in close_pointings if c not in astrometry.excluded]
        if not close_pointings:
            logger.warning(f'HP{healpix:04d}: no pointing left, not mosaiced')
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
        p.prepare_inputs()
        return p.mosaic()
    except Exception:
        logger.exception(f'HP{p.healpix:04d}: unexpected error')
        p.clean()
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
