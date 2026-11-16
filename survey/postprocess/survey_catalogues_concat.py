#!/usr/bin/env python3

# Concatenate the per-HEALPix PyBDSF catalogues into the survey catalogues.
#
# Each mosaic is larger than its HEALPix pixel, so every source is extracted in
# more than one mosaic. A source is taken from the mosaic of the pixel that
# contains its position; Gaussians follow their parent source (they are not
# selected by their own position) and carry its Source_Name.
#
# Mosaics of adjacent pixels are made from different sets of pointings, so the
# position of a source close to a pixel border can fall on different sides in
# the two mosaics. After concatenation the borders are therefore cleaned up:
#   - duplicates: sources from different mosaics closer than --dedup-radius are
#     reduced to the one with the lowest Isl_rms;
#   - losses: a source that each mosaic assigns to the other pixel is recovered.

import numpy as np
from astropy.io import fits
from astropy.table import vstack
from astropy.coordinates import search_around_sky
import glob
import subprocess
import argparse
from astropy.table import Table
from tqdm import tqdm
import os
from astropy_healpix import HEALPix
import astropy.units as u
from astropy.coordinates import SkyCoord

def run(s):
    print(s)
    subprocess.run(s, shell=True, check=True)

parser = argparse.ArgumentParser(description='Concatenate survey catalogues')
parser.add_argument('--catdir', default=os.path.expanduser('~/storage/surveytgts/catalogues/'),
                    help='Directory containing the per-HEALPix catalogue FITS files (default: ~/storage/surveytgts/catalogues/)')
parser.add_argument('--force', action='store_true', default=False,
                    help='Re-process files even if filtered output already exists')
parser.add_argument('--output', default='LoLSS_DR2',
                    help='Base name for the output FITS files (default: LoLSS_DR2)')
parser.add_argument('--outdir', default=None,
                    help='Directory for the filtered and concatenated catalogues (default: same as --catdir)')
parser.add_argument('--dedup-radius', type=float, default=7.5,
                    help='Sources from different mosaics closer than this are the same source [arcsec] (default: 7.5, half a beam)')
parser.add_argument('--border-margin', type=float, default=30.,
                    help='Sources assigned to another pixel are kept as border candidates if closer than this to a pixel edge [arcsec] (default: 30)')
args = parser.parse_args()

CATALOGUE_DIR = args.catdir
force = args.force
    
hp = HEALPix(nside=16)
print(hp.npix,'total healpix pixels on sky')
area=hp.pixel_area.to(u.deg**2).value
print('area of one healpix is',area,'sq. deg')

OUTDIR = os.path.abspath(args.outdir if args.outdir is not None else CATALOGUE_DIR)
os.makedirs(OUTDIR, exist_ok=True)
os.chdir(CATALOGUE_DIR)

CATVER = 2  # format of the filtered files: bump to force their regeneration
ID_FACTOR = 10**7  # (Mosaic_id, Source_id) -> unique integer


def uid(t):
    """Unique identifier of the parent source: mosaic of origin + PyBDSF Source_id."""
    return np.array(t['Mosaic_id'], dtype=np.int64) * ID_FACTOR + np.array(t['Source_id'], dtype=np.int64)


def source_names(ra, dec):
    sc = SkyCoord(ra, dec, frame='icrs')
    strings = sc.to_string(style='hmsdms', sep='', precision=2)
    return [('LoLJ'+s).replace(' ', '')[:-1] for s in strings]


def near_border(ra, dec, pix, margin):
    """True for positions closer than margin [arcsec] to an edge of their HEALPix pixel."""
    sc = SkyCoord(ra, dec, frame='icrs')
    near = np.zeros(len(sc), dtype=bool)
    for pa in np.arange(0, 360, 45):
        off = sc.directional_offset_by(pa*u.deg, margin*u.arcsec)
        near |= hp.lonlat_to_healpix(off.ra, off.dec) != pix
    return near


def clean_columns(t):
    """Remove the columns that are not released and convert to the catalogue units."""
    for k in t.colnames:
        if 'img_plane' in k or k.startswith('Resid') or k.endswith('_max') or k in ['Flag_beam', 'Isl_mean', 'Isl_id']:
            del t[k]
    for k in ['E_RA', 'E_DEC', 'Maj', 'E_Maj', 'Min', 'E_Min', 'DC_Maj', 'E_DC_Maj', 'DC_Min', 'E_DC_Min']:
        t[k].convert_unit_to(u.arcsec)
    for k in ['Peak_flux', 'E_Peak_flux', 'Isl_rms']:
        t[k].convert_unit_to(u.mJy/u.beam)
    for k in ['Total_flux', 'E_Total_flux', 'Isl_Total_flux', 'E_Isl_Total_flux']:
        t[k].convert_unit_to(u.mJy)
    t.meta = {'CATVER': CATVER, 'INIMAGE': t.meta.get('INIMAGE', '')}
    return t


def is_current(filename):
    """True if a filtered file exists and was written in the current format."""
    if not os.path.isfile(filename):
        return False
    return fits.getheader(filename, 1).get('CATVER') == CATVER


def process_pixel(pix):
    """
    Split the source and Gaussian catalogues of mosaic HP<pix> in
      own    : sources whose position is in this HEALPix pixel
      border : sources assigned to another pixel but within the margin of the edge
    and write them as HP<pix>-mosaicI.{cat,gaus}-{filtered,border}.fits in OUTDIR.
    Gaussians are selected through their parent source and get its Source_Name.
    """
    stem = f'HP{pix:04d}-mosaicI'
    out = {(cattype, part): os.path.join(OUTDIR, f'{stem}.{cattype}-{part}.fits')
           for cattype in ('cat', 'gaus') for part in ('filtered', 'border')}
    if not force and is_current(out['cat', 'filtered']) and is_current(out['gaus', 'filtered']):
        return
    for outname in out.values():
        if os.path.isfile(outname):
            os.remove(outname)

    cat = Table.read(f'{stem}.cat.fits')
    gaus = Table.read(f'{stem}.gaus.fits')

    cat['HEALPIX'] = hp.lonlat_to_healpix(cat['RA'].data*u.deg, cat['DEC'].data*u.deg)
    cat['Mosaic_id'] = np.int32(pix)
    own = np.array(cat['HEALPIX'] == pix)
    border = ~own & near_border(cat['RA'].data*u.deg, cat['DEC'].data*u.deg,
                                np.array(cat['HEALPIX']), args.border_margin)
    if not own.any():
        print(f'*** no source in its own pixel for {stem} ***')
    cat = clean_columns(cat[own | border])
    cat.add_column(source_names(cat['RA'], cat['DEC']), name='Source_Name', index=0)
    own = np.array(cat['HEALPIX'] == pix)

    # Gaussians inherit name and pixel from their parent source
    parent = {sid: i for i, sid in enumerate(cat['Source_id'])}
    gaus = gaus[np.isin(gaus['Source_id'], cat['Source_id'])]
    gaus = clean_columns(gaus)
    idx = np.array([parent[sid] for sid in gaus['Source_id']], dtype=int)
    gaus.add_column(cat['Source_Name'][idx], name='Source_Name', index=0)
    gaus['HEALPIX'] = cat['HEALPIX'][idx]
    gaus['Mosaic_id'] = np.int32(pix)
    gaus_own = np.array(gaus['HEALPIX'] == pix)

    for part, sel_cat, sel_gaus in (('filtered', own, gaus_own), ('border', ~own, ~gaus_own)):
        if sel_cat.any():
            cat[sel_cat].write(out['cat', part], overwrite=True)
            gaus[sel_gaus].write(out['gaus', part], overwrite=True)


def read_all(cattype, part):
    files = sorted(glob.glob(os.path.join(OUTDIR, f'HP*-mosaicI.{cattype}-{part}.fits')))
    with open(os.path.join(OUTDIR, f'{cattype}-catlist.txt' if part == 'filtered' else f'{cattype}-borderlist.txt'), 'w') as flist:
        flist.write('\n'.join(files) + '\n')
    if not files:
        return None
    return vstack([Table.read(f) for f in tqdm(files, desc=f'{cattype}-{part}')], metadata_conflicts='silent')


def pairs_within(coords1, coords2, radius):
    """Indices (i1, i2) of all the pairs closer than radius [arcsec]."""
    if len(coords1) == 0 or len(coords2) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    i1, i2, _, _ = search_around_sky(coords1, coords2, radius*u.arcsec)
    return i1, i2


def recover_lost(cat, border, radius):
    """
    Return the rows of border to add to cat: sources without a counterpart in cat
    that the mosaic of their own pixel also detected but assigned to another pixel.
    Duplicates among the recovered rows are left to remove_duplicates().
    """
    if border is None or len(border) == 0:
        return np.zeros(0, dtype=bool)
    cb = SkyCoord(border['RA'], border['DEC'], unit='deg')
    i1, _ = pairs_within(cb, SkyCoord(cat['RA'], cat['DEC'], unit='deg'), radius)
    unmatched = np.ones(len(border), dtype=bool)
    unmatched[i1] = False
    # seen by the mosaic of the pixel it lies in?
    i1, i2 = pairs_within(cb, cb, radius)
    seen_at_home = np.zeros(len(border), dtype=bool)
    home = np.array(border['Mosaic_id'])[i2] == np.array(border['HEALPIX'])[i1]
    seen_at_home[i1[home]] = True
    return unmatched & seen_at_home


def remove_duplicates(cat, radius):
    """
    Return a boolean mask of the rows to keep. Two sources from different mosaics
    closer than radius are the same source: the one with the lowest Isl_rms is kept.
    """
    coords = SkyCoord(cat['RA'], cat['DEC'], unit='deg')
    i1, i2 = pairs_within(coords, coords, radius)
    mosaic = np.array(cat['Mosaic_id'])
    sel = (i1 != i2) & (mosaic[i1] != mosaic[i2])
    neighbours = {}
    for a, b in zip(i1[sel], i2[sel]):
        neighbours.setdefault(int(a), set()).add(int(b))
    keep = np.ones(len(cat), dtype=bool)
    rms = np.array(cat['Isl_rms'])
    kept = set()
    # best sources first: a source is dropped if a better duplicate was already kept
    for a in sorted(neighbours, key=lambda i: (rms[i], mosaic[i])):
        if neighbours[a] & kept:
            keep[a] = False
        else:
            kept.add(a)
    return keep


pixels = sorted(int(os.path.basename(f)[2:6]) for f in glob.glob('HP*-mosaicI.cat.fits'))
for pix in tqdm(pixels, desc='filtering'):
    process_pixel(pix)

cat, gaus = read_all('cat', 'filtered'), read_all('gaus', 'filtered')
border_cat, border_gaus = read_all('cat', 'border'), read_all('gaus', 'border')
if cat is None:
    raise SystemExit(f'No filtered catalogue found in {OUTDIR}')
print(f'Concatenated: {len(cat)} sources, {len(gaus)} Gaussians')

# recover the sources that no mosaic assigned to its own pixel
recover = recover_lost(cat, border_cat, args.dedup_radius)
log = []
if recover.any():
    recovered = border_cat[recover]
    cat = vstack([cat, recovered], metadata_conflicts='silent')
    gaus = vstack([gaus, border_gaus[np.isin(uid(border_gaus), uid(recovered))]], metadata_conflicts='silent')
recovered_uid = uid(border_cat[recover]) if recover.any() else np.zeros(0, dtype=np.int64)

# remove the sources extracted twice
keep = remove_duplicates(cat, args.dedup_radius)
is_recovered = np.isin(uid(cat), recovered_uid)
for action, sel in (('removed', ~keep & ~is_recovered), ('recovered', keep & is_recovered)):
    t = cat['Source_Name', 'RA', 'DEC', 'Total_flux', 'Isl_rms', 'HEALPIX', 'Mosaic_id'][sel]
    t['Action'] = action
    log.append(t)
print(f'Border clean-up within {args.dedup_radius}": {np.sum(~keep & ~is_recovered)} duplicates removed, '
      f'{np.sum(keep & is_recovered)} lost sources recovered')
cat = cat[keep]
gaus = gaus[np.isin(uid(gaus), uid(cat))]

names, counts = np.unique(cat['Source_Name'], return_counts=True)
if (counts > 1).any():
    print(f'WARNING: {np.sum(counts > 1)} Source_Name are not unique, e.g. {list(names[counts > 1][:5])}')

log = vstack(log, metadata_conflicts='silent')
log.sort('RA')
log.write(os.path.join(OUTDIR, f'{args.output}-border-log.fits'), overwrite=True)

for cattype, t in (('cat', cat), ('gaus', gaus)):
    del t['Source_id'], t['Mosaic_id']
    t.sort('RA')
    t.meta = {'EXTNAME': f'{args.output}_{cattype}'}
    outname = os.path.join(OUTDIR, f'{args.output}-{cattype}.fits')
    t.write(outname, overwrite=True)
    print(f'Written {outname}: {len(t)} rows')
