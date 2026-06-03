#!/usr/bin/env python3

# For each HEALPix pixel defined in pointing_sets_file, collect the associated
# survey pointings and produce a Stokes-I mosaic using mosaic.py with beam
# correction. An empty FITS template is created first to define the output WCS
# and pixel grid; mosaics that are already complete (all pointings recorded in
# the FITS HISTORY header) are skipped.

import os, logging, argparse
from astropy.io import fits
from astropy import units as u
from astropy_healpix import HEALPix
import numpy as np
from LiLF.surveys_db import SurveysDB

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Mosaic survey pointings onto HEALPix tiles.')
parser.add_argument('--dir-mosaics',        default='/homes/fdg/storage/surveytgts/mosaics/healpix/',                    help='Output directory for mosaics')
parser.add_argument('--dir-done',           default='/homes/fdg/storage/surveytgts/done/',                       help='Directory containing completed pointings')
parser.add_argument('--file-image',         default='wideDDS-c0-MFS-image.fits',                                 help='MFS image filename within each pointing directory')
parser.add_argument('--nchans',             default=None, type=int,                                              help='Number of channels (default: None = MFS only)')
parser.add_argument('--file-beam',          default='primarybeam.fits',                                          help='Primary beam filename within each pointing directory')
parser.add_argument('--pointing-sets-file', default='/homes/fdg/storage/surveytgts/analysis/pointing_sets.txt', help='File mapping HEALPix pixels to pointing names')
parser.add_argument('--only-pointing',      default=None, type=int, nargs='+', metavar='HP',                    help='Only mosaic these HEALPix pixel IDs (default: all)')
parser.add_argument('--cellsize',           default=4, type=int,                                                 help='Output pixel size in arcsec (default: 4)')
parser.add_argument('--resolution',         default=15, type=int,                                                help='Output beam FWHM in arcsec (default: 15)')
parser.add_argument('--sbatch',             action='store_true',                                                 help='Submit mosaic jobs via sbatch instead of running locally')
args = parser.parse_args()

dir_mosaics        = args.dir_mosaics
dir_done           = args.dir_done
file_MFS_i         = args.file_image
nchans             = args.nchans
file_beam          = args.file_beam
pointing_sets_file = args.pointing_sets_file
only_this_pointing = args.only_pointing
CELLSIZE           = args.cellsize
RESOLUTION         = args.resolution
use_sbatch         = args.sbatch

hp = HEALPix(nside=16, frame='icrs')
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

    def mosaic(self):
        os.makedirs(os.path.dirname(self.log), exist_ok=True)
        cmd = 'mosaic.py --images %s --beams %s --beamcorr --beamcut 0.4 --find_noise --header %s --output %s > %s 2>&1' % \
                (' '.join(self.pointing_files), ' '.join(self.beam_files), self.output_file, self.output_file, self.log)
        logger.info(cmd)
        if use_sbatch:
            os.makedirs(os.path.dirname(self.batch_script), exist_ok=True)
            cmd = singularity_cmd + ' ' + cmd.replace('mosaic.py', '/homes/fdg/storage/scripts/mosaic.py')
            with open(self.batch_script, 'w') as _f:
                _f.write('#!/bin/bash\n')
                _f.write(f'#SBATCH --nodes=1\n')
                _f.write(f'#SBATCH --partition=lofar\n')
                _f.write(f'#SBATCH --ntasks-per-core=1\n')
                _f.write(f'#SBATCH --ntasks-per-node=1\n')
                _f.write(f'#SBATCH --cpus-per-task=36\n')
                _f.write(f'#SBATCH --time=240:00:00\n')
                _f.write(f'\n{cmd}\n')
            os.system(f'sbatch {self.batch_script}')
        else:
            os.system(cmd)

    def make_empty_mosaic(self):
        centre = hp.healpix_to_skycoord(self.healpix)
        corners = hp.boundaries_skycoord(self.healpix, step=1)
        dra, ddec = centre.spherical_offsets_to(corners)
        dra = dra.value
        ddec = ddec.value
        maxdist=np.max(np.abs([dra,ddec]))
        size = maxdist*1.15
        himsize = int(size/(CELLSIZE / 3600.0))
        logger.info(f'Creating empty mosaic for {self.output_file} - size: {2*size:.2f} deg, himsize: {2*himsize} pixels')

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
        header['RESTFRQ']=54e6
        header['OBSERVER']='LoLSS'
        header['BUNIT']='JY/BEAM'
        header['BSCALE']=1.0
        header['BZERO']=0
        header['BTYPE']='Intensity'
        header['OBJECT']='HP-'+str(self.healpix)
        for p in self.pointings_name_nearby:
            header['HISTORY'] = 'POINTING ADDED: ' + p
        data = np.zeros((2*himsize,2*himsize), dtype=np.float32)
        fits.writeto(self.output_file, header=header, data=data, overwrite=True)

    def check_todo(self):
        if not os.path.exists(self.output_file):
            return True
        # check if the mosaic is complete by looking at the header
        with fits.open(self.output_file) as _f:
            _hdr = _f[0].header
        for p in self.pointings_name_nearby:
            if 'POINTING ADDED: ' + p not in _hdr['HISTORY']:
                logger.info(f'Pointing {p} not found in HISTORY of {self.output_file}, needs to be remade')
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
for healpix, close_pointings in grouped_pointings.items():
    # skip if not in the list of pointings to mosaic
    if only_this_pointing is not None and healpix not in only_this_pointing: 
        logger.info(f'Skipping pointing {healpix} (not in {only_this_pointing})')
        continue

    p = Pointing(healpix, close_pointings)
    if p.check_todo():
        logger.info('Pointing: %s %s', healpix, close_pointings)
        coord = hp.healpix_to_skycoord(healpix)
        ra = coord.ra.deg
        dec = coord.dec.deg
        p.make_empty_mosaic()
        p.mosaic()
    else:
        logger.info(f'Skipping pointing {healpix} (already complete)')