#!/usr/bin/env python3
"""
survey_astrometry.py

Astrometry of the LoLSS pointings before mosaicing: offsets of every pointing
and of every facet of every pointing against LoTSS-DR3, measured on the
single-pointing 15" images (wideDDS-c0-MFS-image.fits).

Outputs (in --outdir, default analysis/astrometry/):
  catalogues/<pointing>.srl.fits  PyBDSF source list of each pointing
                                  (made once; --redo-bdsf to make it again)
  pointings/<pointing>.fits       LoLSS-DR3 matches with their facet (cache;
                                  --redo to match again)
  astrometry_pointings.fits       HDU POINTINGS: global offset of each pointing
                                  HDU FACETS:    offset of each facet, flags
  astrometry_pointings.pdf        offsets of the pointings
  astrometry_facets.pdf           offsets of the facets
The corrections are applied by survey_mosaics-i.py.

Method
------
PyBDSF (survey parameters of lolss_config.py) is run on the flat-noise image
of each pointing, through a symbolic link in catalogues/, so that nothing is
written in the pointing directory. The LoLSS sources used have S_Code in
--lolss-codes (default S and M), S/N >= --snr-min, Maj < --lolss-maj-max beams,
no other LoLSS source within --lolss-iso-radius and primary beam >= BEAM_THRESHOLD.
Reference sources are LoTSS-DR3 sources that are compact (Maj < --ref-maj-max)
and isolated (no DR3 source brighter than --ref-iso-frac times their flux within
--ref-iso-radius). LoLSS and DR3 sources are matched as mutual nearest
neighbours within --match-radius. Each match is assigned to the facet of
facetsS-c0.reg that contains the LoLSS position.

All offsets are reference - LoLSS, in arcsec on the sky: dRA is the RA
difference times cos(Dec). They are the shift to apply to the LoLSS positions;
the corresponding change of CRVAL1 is dRA / cos(Dec).

Only sources with S/N >= --snr-cut are used. Medians are weighted with
1/(PyBDSF position error^2 + --sys-err^2). In each pointing, iteratively:
  - pointing offset = median offset of the sources in the facets not flagged bad;
  - single-source scatter = robust scatter of the sources around the offset of
    their own facet, from the facets with >= --facet-nmin-scatter sources (so
    that differences between facets do not count);
  - every facet gets its median offset and the residual (facet - pointing offset);
    the error of the residual comes from the single-source scatter / sqrt(N),
    so that also facets with few sources can be judged;
  - a facet is flagged OFFSET if the residual exceeds --facet-max-offset and is
    significant at --facet-nsigma, SCATTER if it has at least --facet-nmin-scatter
    sources and their robust scatter exceeds --facet-max-scatter, FEW if it has
    fewer than --facet-nmin sources (FEW facets are not judged, but their sources
    do count for the pointing offset).
A pointing is flagged FEW if fewer than --pointing-nmin sources remain,
SCATTER if the robust scatter of their offsets exceeds --pointing-max-scatter
(only if not FEW; e.g. facets shifted in different directions, so that the
median is not meaningful; FEW pointings are not corrected and SCATTER pointings
are left out of the mosaics by survey_mosaics-i.py), and
OUTLIER if its offset differs from the survey median by more than
--pointing-max-offset with a significance of --pointing-nsigma.

PyBDSF is the slow step (minutes per pointing): use -j. The match is redone
automatically when a catalogue is newer than its match file. The analysis is
done again at every run from the cache, so
the thresholds can be changed cheaply with --analyse-only.
"""

import os, re, argparse, glob, time, warnings
import multiprocessing
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

import lolss_analysis_common as common
import lolss_config as cfg

from astropy.wcs import FITSFixedWarning
warnings.simplefilter('ignore', FITSFixedWarning)

# facet flags (bit mask)
FLAG_FEW     = 1   # too few sources to judge the facet
FLAG_OFFSET  = 2   # facet offset differs from the pointing offset
FLAG_SCATTER = 4   # facet positions are scattered (distorted facet)
FLAG_OUTSIDE = 8   # pseudo-facet of the sources outside every polygon
# pointing flags (bit mask)
PFLAG_FEW     = 1  # too few sources for a reliable offset
PFLAG_OUTLIER = 2  # offset differs from the survey median
PFLAG_NOFACETS = 4 # no facetsS-c0.reg: facets not checked
PFLAG_SCATTER = 8  # source offsets too scattered: the median offset is not meaningful
PFLAG_BAD = PFLAG_FEW | PFLAG_SCATTER  # offset not reliable: not used, not corrected


############################################################################
# pointings: offsets of pointings and facets against LoTSS-DR3

def unit_vectors(ra, dec):
    ra, dec = np.radians(np.asarray(ra, dtype=float)), np.radians(np.asarray(dec, dtype=float))
    return np.column_stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])


def chord(arcsec):
    """Chord on the unit sphere of an angle in arcsec (for KD-tree queries)."""
    return 2*np.sin(np.radians(np.asarray(arcsec, dtype=float)/3600.)/2)


def mad_sigma(x):
    x = np.asarray(x, dtype=float)
    if len(x) == 0:
        return np.nan
    return 1.4826*np.median(np.abs(x - np.median(x)))


def sky_offset(ra_ref, dec_ref, ra, dec):
    """Reference - measured, in arcsec on the sky (RA difference times cos(Dec))."""
    dra = (np.asarray(ra_ref) - np.asarray(ra) + 180.) % 360. - 180.
    return dra*np.cos(np.radians(dec_ref))*3600., (np.asarray(dec_ref) - np.asarray(dec))*3600.


def load_dr3(path_cache, path_cat):
    """Positions, fluxes (mJy) and sizes (arcsec) of the DR3 sources, sorted by Dec (cached in an npz)."""
    if not os.path.exists(path_cache):
        print('Reading %s (once, cached in %s)...' % (path_cat, path_cache), flush=True)
        with fits.open(path_cat, memmap=True) as h:
            d = h[1].data
            dec = np.array(d['DEC'], dtype=np.float64)
            order = np.argsort(dec)
            out = dict(dec=dec[order], ra=np.array(d['RA'], dtype=np.float64)[order],
                       flux=np.array(d['Total_flux'], dtype=np.float32)[order],
                       maj=np.array(d['Maj'], dtype=np.float32)[order])
        tmp = path_cache + '.tmp%i.npz' % os.getpid()
        np.savez(tmp, **out)
        os.replace(tmp, path_cache)
    with np.load(path_cache) as c:
        return {k: c[k] for k in ('ra', 'dec', 'flux', 'maj')}


def read_facets(regfile, wcs):
    """Facet polygons of a ds9 region file: list of (name, matplotlib Path in pixels)."""
    from matplotlib.path import Path
    facets = []
    with open(regfile) as f:
        for line in f:
            m = re.match(r'\s*polygon\(([^)]*)\)', line)
            if m:
                v = np.array(m.group(1).split(','), dtype=float).reshape(-1, 2)
                x, y = wcs.all_world2pix(v[:, 0], v[:, 1], 0)
                facets.append(['facet%03i' % len(facets), Path(np.column_stack([x, y]))])
                continue
            m = re.search(r'text\s*=\s*[{"]([^}"]*)[}"]', line)
            if m and line.lstrip().startswith('point') and facets:
                facets[-1][0] = m.group(1)   # the name of the calibrator follows its polygon
    return facets


# DR3 catalogue, set in the parent before the pool is created (shared by fork)
DR3 = None


def catalogue_path(pointing, args):
    return os.path.join(args.outdir, 'catalogues', pointing + '.srl.fits')


def find_sources(pointing, args):
    """
    PyBDSF source list of one pointing in <outdir>/catalogues/<pointing>.srl.fits.
    PyBDSF runs on a symbolic link placed in that directory, so that its log is
    written there and nothing is ever written in the pointing directory.
    """
    import bdsf
    catdir = os.path.join(args.outdir, 'catalogues', '')
    out = catalogue_path(pointing, args)
    link = catdir + pointing + '-image.fits'
    if os.path.lexists(link):
        os.remove(link)
    os.symlink(os.path.abspath(os.path.join(args.dir_done, pointing, args.file_image)), link)
    t0 = time.time()
    kwargs = cfg.bdsf_kwargs(frequency=cfg.RESTFRQ, quiet=True)
    if args.bdsf_ncores:
        kwargs['ncores'] = args.bdsf_ncores
    img = bdsf.process_image(link, **kwargs)
    tmp = catdir + pointing + '.srl.tmp%i.fits' % os.getpid()
    if img.nsrc > 0:
        img.write_catalog(outfile=tmp, catalog_type='srl', format='fits', clobber=True)
    else:
        Table(names=['RA', 'DEC'], dtype=[float, float]).write(tmp, overwrite=True)
    os.replace(tmp, out)
    os.remove(link)
    print('%s: PyBDSF found %i sources (%.0f s)' % (pointing, img.nsrc, time.time()-t0), flush=True)
    return out


def match_pointing(pointing, args):
    """Match the PyBDSF sources of one pointing to DR3; write and return the cache file name (None if not possible)."""
    from astropy.wcs import WCS
    from scipy.ndimage import map_coordinates
    from scipy.spatial import cKDTree

    d = os.path.join(args.dir_done, pointing, '')
    file_img, file_pb, file_reg = d + args.file_image, d + cfg.FILE_BEAM, d + args.file_facets
    out = os.path.join(args.outdir, 'pointings', pointing + '.fits')
    file_cat = catalogue_path(pointing, args)

    hdr = fits.getheader(file_img)
    wcs = WCS(hdr).celestial
    ra0, dec0 = hdr['CRVAL1'] % 360, hdr['CRVAL2']
    beam = hdr.get('BMAJ', args.beam/3600.)*3600.
    with fits.open(file_pb) as h:
        pb = np.array(np.squeeze(h[0].data), dtype=np.float32)
        pb_wcs = WCS(h[0].header).celestial
        pb_pix = abs(h[0].header['CDELT2'])
    yy, xx = np.indices(pb.shape)
    cy, cx = pb_wcs.wcs.crpix[1]-1, pb_wcs.wcs.crpix[0]-1
    radius = (np.hypot(xx-cx, yy-cy)[pb >= args.beamcut].max() + 2)*pb_pix if (pb >= args.beamcut).any() else 0
    del yy, xx

    def pb_at(ra, dec):
        x, y = pb_wcs.all_world2pix(ra, dec, 0)
        return map_coordinates(pb, [y, x], order=1, mode='constant', cval=0.)

    # LoLSS: compact, isolated sources above the beam threshold
    cat = Table.read(file_cat)
    n_all = len(cat)
    if n_all:
        lra, ldec = np.array(cat['RA'], dtype=float) % 360, np.array(cat['DEC'], dtype=float)
        lv = unit_vectors(lra, ldec)
        iso = np.ones(n_all, dtype=bool)
        if n_all > 1:
            dist, _ = cKDTree(lv).query(lv, k=2)
            iso = dist[:, 1] > chord(args.lolss_iso_radius)
        snr = np.array(cat['Peak_flux'])/np.array(cat['Isl_rms'])
        scode = np.char.strip(np.array(cat['S_Code']).astype(str))
        sel = np.isin(scode, list(args.lolss_codes)) & iso & (snr >= args.snr_min) & (np.array(cat['Maj'])*3600 < args.lolss_maj_max*beam) \
            & (pb_at(lra, ldec) >= args.beamcut)
        cat, lra, ldec, lv, snr = cat[sel], lra[sel], ldec[sel], lv[sel], snr[sel]

    # DR3: compact and isolated sources (the isolation is checked against the whole DR3 catalogue)
    margin = args.ref_iso_radius/3600.
    lo, hi = np.searchsorted(DR3['dec'], [dec0 - radius - margin, dec0 + radius + margin])
    ra, dec, flux, maj = (DR3[k][lo:hi] for k in ('ra', 'dec', 'flux', 'maj'))
    v = unit_vectors(ra, dec)
    keep = v @ unit_vectors([ra0], [dec0])[0] > np.cos(np.radians(radius + margin))
    ra, dec, flux, maj, v = ra[keep], dec[keep], flux[keep], maj[keep], v[keep]
    ref = np.where((maj < args.ref_maj_max) & (flux > args.ref_flux_min))[0]
    if len(ref):
        groups = cKDTree(v).query_ball_point(v[ref], chord(args.ref_iso_radius))
        iso = np.array([np.all(flux[[j for j in g if j != i]] < args.ref_iso_frac*flux[i]) for i, g in zip(ref, groups)])
        ref = ref[iso]

    # mutual nearest neighbours within the match radius
    rows = []
    if len(cat) and len(ref):
        rv = v[ref]
        d12, i12 = cKDTree(rv).query(lv)
        _, i21 = cKDTree(lv).query(rv)
        ok = (i21[i12] == np.arange(len(lv))) & (d12 < chord(args.match_radius))
        for k in np.where(ok)[0]:
            i = ref[i12[k]]
            dra, ddec = sky_offset(ra[i], dec[i], lra[k], ldec[k])
            c = cat[k]
            rows.append((ra[i], dec[i], flux[i], maj[i], lra[k], ldec[k], c['E_RA']*3600, c['E_DEC']*3600, dra, ddec,
                         snr[k], c['Total_flux'], c['Peak_flux'], c['Isl_rms'], c['Maj']*3600, c['Min']*3600,
                         float(pb_at([lra[k]], [ldec[k]])[0])))
    names = ['RA_ref', 'DEC_ref', 'Flux_ref', 'Maj_ref', 'RA', 'DEC', 'E_RA', 'E_DEC', 'dRA', 'dDEC',
             'SNR', 'Total_flux', 'Peak_flux', 'Isl_rms', 'Maj', 'Min', 'PB']
    tab = Table(rows=rows, names=names) if rows else Table(names=names, dtype=[float]*len(names))

    # facet of each source (-1 if outside every polygon)
    facets = read_facets(file_reg, wcs) if os.path.exists(file_reg) else []
    fidx = np.full(len(tab), -1, dtype=np.int16)
    if len(tab) and facets:
        x, y = wcs.all_world2pix(np.array(tab['RA']), np.array(tab['DEC']), 0)
        pts = np.column_stack([x, y])
        for j, (_, path) in enumerate(facets):
            s = (fidx < 0) & path.contains_points(pts)
            fidx[s] = j
    tab['Facet'] = fidx

    # facet table: name and centre (centroid of the polygon vertices)
    ftab = Table(names=['Facet', 'Name', 'RA', 'DEC'], dtype=[np.int16, 'U32', float, float])
    for j, (name, path) in enumerate(facets):
        cxy = path.vertices.mean(axis=0)
        fra, fdec = wcs.all_pix2world(cxy[0], cxy[1], 0)
        ftab.add_row((j, name, float(fra) % 360, float(fdec)))

    hdr_out = fits.Header()
    hdr_out['POINTING'] = pointing
    hdr_out['RA'] = (ra0, 'pointing centre [deg]')
    hdr_out['DEC'] = (dec0, 'pointing centre [deg]')
    hdr_out['NLOLSS'] = (n_all, 'sources in the PyBDSF catalogue')
    hdr_out['NLOLSSOK'] = (len(cat), 'LoLSS sources usable for the match')
    hdr_out['NREF'] = (len(ref), 'DR3 reference sources')
    hdr_out['NFACETS'] = (len(facets), 'facets in ' + args.file_facets)
    hdr_out['IMAGE'] = args.file_image
    hdu_src = fits.BinTableHDU(tab, name='SOURCES')
    hdu_src.header.update(hdr_out)
    hdu_fac = fits.BinTableHDU(ftab, name='FACETS')
    tmp = out + '.tmp%i' % os.getpid()
    fits.HDUList([fits.PrimaryHDU(header=hdr_out), hdu_src, hdu_fac]).writeto(tmp, overwrite=True)
    os.replace(tmp, out)
    print('%s: %i matches (LoLSS usable %i of %i, DR3 references %i), %i facets'
          % (pointing, len(tab), len(cat), n_all, len(ref), len(facets)), flush=True)
    return out


def measure_pointing(pointing, args):
    """Source finding (if needed) and match of one pointing."""
    d = os.path.join(args.dir_done, pointing, '')
    for f in (d + args.file_image, d + cfg.FILE_BEAM):
        if not os.path.exists(f):
            print('%s: missing %s, skipped' % (pointing, f), flush=True)
            return None
    if args.redo_bdsf or not os.path.exists(catalogue_path(pointing, args)):
        find_sources(pointing, args)
    return match_pointing(pointing, args)


def _measure(item):
    pointing, args = item
    try:
        return measure_pointing(pointing, args)
    except Exception as e:
        import traceback
        print('%s: FAILED %s\n%s' % (pointing, e, traceback.format_exc()), flush=True)
        return None


def wmedian(x, w):
    """Weighted median."""
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return x[o][np.searchsorted(c, 0.5*c[-1])]


def robust_offset(dra, ddec, w=None, sigma=None):
    """
    (Weighted) median offset, robust scatter of the sources and error of the median, per axis.
    sigma: scatter of a single source to use for the error (default: the robust scatter of these sources).
    """
    n = len(dra)
    if n == 0:
        return (np.nan,)*6
    if w is None:
        w = np.ones(n)
    sr, sd = mad_sigma(dra), mad_sigma(ddec)
    er, ed = (sr, sd) if sigma is None else sigma
    k = 1.2533/np.sqrt(np.sum(w)**2/np.sum(w**2))   # error of a median, with the effective number of sources
    return wmedian(dra, w), wmedian(ddec, w), sr, sd, k*er, k*ed


def analyse(args):
    """Offsets and flags of pointings and facets from the cache files."""
    files = sorted(glob.glob(os.path.join(args.outdir, 'pointings', '*.fits')))
    if args.only:
        files = [f for f in files if os.path.basename(f)[:-5] in args.only]
    prow, frow = [], []
    for f in files:
        with fits.open(f) as h:
            hdr = h[0].header
            src = Table(h['SOURCES'].data)
            fac = Table(h['FACETS'].data)
        pointing = hdr['POINTING']
        src = src[src['SNR'] >= args.snr_cut]
        dra, ddec, fidx = np.array(src['dRA']), np.array(src['dDEC']), np.array(src['Facet'])
        # weights: PyBDSF position error plus a systematic floor, so that no single bright source dominates
        wgt = 1./((np.array(src['E_RA'])**2 + np.array(src['E_DEC'])**2)/2 + args.sys_err**2)
        has_facets = len(fac) > 0
        fids = list(range(len(fac))) + ([-1] if np.any(fidx < 0) else [])

        bad = set()
        for it in range(5):
            use = ~np.isin(fidx, list(bad))
            p_ra, p_dec, p_sra, p_sdec, p_era, p_edec = robust_offset(dra[use], ddec[use], wgt[use])
            # scatter of a single source around the offset of its own facet (facet-to-facet differences
            # removed), from the facets with enough sources: it gives the error of the facets with few sources
            r_ra, r_dec = [], []
            for j in fids:
                m = (fidx == j) & use
                if j >= 0 and m.sum() >= args.facet_nmin_scatter:
                    corr = np.sqrt(m.sum()/(m.sum() - 1.))
                    r_ra.append((dra[m] - wmedian(dra[m], wgt[m]))*corr); r_dec.append((ddec[m] - wmedian(ddec[m], wgt[m]))*corr)
            if sum(len(r) for r in r_ra) >= 10:
                s_src = (mad_sigma(np.concatenate(r_ra)), mad_sigma(np.concatenate(r_dec)))
            else:
                s_src = (p_sra, p_sdec)
            stats = {}
            newbad = set()
            for j in fids:
                m = fidx == j
                n = int(m.sum())
                f_ra, f_dec, f_sra, f_sdec, f_era, f_edec = robust_offset(dra[m], ddec[m], wgt[m], sigma=s_src)
                res_ra, res_dec = f_ra - p_ra, f_dec - p_dec
                res = np.hypot(res_ra, res_dec)
                e_res = np.hypot(f_era, f_edec)/np.sqrt(2) if n else np.nan
                scatter = np.sqrt((f_sra**2 + f_sdec**2)/2) if n >= args.facet_nmin_scatter else np.nan
                flag = FLAG_OUTSIDE if j < 0 else 0
                if n < args.facet_nmin:
                    flag |= FLAG_FEW
                else:
                    if res > args.facet_max_offset and res > args.facet_nsigma*e_res:
                        flag |= FLAG_OFFSET
                    if scatter > args.facet_max_scatter:   # nan (too few sources to judge) is never flagged
                        flag |= FLAG_SCATTER
                    if flag & (FLAG_OFFSET | FLAG_SCATTER) and j >= 0:
                        newbad.add(j)
                stats[j] = (n, f_ra, f_dec, f_sra, f_sdec, f_era, f_edec, res_ra, res_dec, res, e_res, scatter, flag)
            if newbad == bad:
                break
            bad = newbad

        n_good = int((~np.isin(fidx, list(bad))).sum())
        pflag = 0 if has_facets else PFLAG_NOFACETS
        if n_good < args.pointing_nmin:
            pflag |= PFLAG_FEW
        elif np.sqrt((p_sra**2 + p_sdec**2)/2) > args.pointing_max_scatter:   # judged only with enough sources
            pflag |= PFLAG_SCATTER
        nfac = len(fac)
        n_meas = sum(1 for j in range(nfac) if not stats[j][-1] & FLAG_FEW)
        n_bad = sum(1 for j in range(nfac) if j in bad)
        prow.append((pointing, hdr['RA'], hdr['DEC'], len(src), n_good, nfac, n_meas, n_bad,
                     p_ra, p_dec, p_era, p_edec, p_sra, p_sdec, s_src[0], s_src[1], pflag))
        for j in fids:
            if j >= 0:
                name, fra, fdec = fac['Name'][j], fac['RA'][j], fac['DEC'][j]
            else:
                name, fra, fdec = 'outside', np.nan, np.nan
            dist = np.degrees(np.arccos(np.clip(unit_vectors([fra], [fdec])[0] @ unit_vectors([hdr['RA']], [hdr['DEC']])[0], -1, 1))) if j >= 0 else np.nan
            frow.append((pointing, j, name, fra, fdec, dist) + stats[j])

    ptab = Table(rows=prow, names=['Pointing', 'RA', 'DEC', 'N_src', 'N_src_good', 'N_facets', 'N_facets_measured', 'N_facets_bad',
                                   'dRA', 'dDEC', 'E_dRA', 'E_dDEC', 'Sigma_dRA', 'Sigma_dDEC', 'Sigma_src_dRA', 'Sigma_src_dDEC', 'Flag'])
    ftab = Table(rows=frow, names=['Pointing', 'Facet', 'Name', 'RA', 'DEC', 'Dist', 'N_src',
                                   'dRA', 'dDEC', 'Sigma_dRA', 'Sigma_dDEC', 'E_dRA', 'E_dDEC',
                                   'Res_dRA', 'Res_dDEC', 'Res', 'E_Res', 'Scatter', 'Flag'])

    # pointings that disagree with the rest of the survey
    ok = (ptab['Flag'] & PFLAG_BAD) == 0
    s_ra, s_dec = np.median(ptab['dRA'][ok]), np.median(ptab['dDEC'][ok])
    dev = np.hypot(ptab['dRA'] - s_ra, ptab['dDEC'] - s_dec)
    e_dev = np.hypot(ptab['E_dRA'], ptab['E_dDEC'])/np.sqrt(2)
    out = ok & (dev > args.pointing_max_offset) & (dev > args.pointing_nsigma*e_dev)
    ptab['Flag'][out] |= PFLAG_OUTLIER
    ptab['Offset'] = np.hypot(ptab['dRA'], ptab['dDEC'])
    ptab['dCRVAL1'] = ptab['dRA']/np.cos(np.radians(ptab['DEC']))/3600.   # deg, to add to CRVAL1
    ptab['dCRVAL2'] = ptab['dDEC']/3600.

    for c in ('RA', 'DEC', 'Dist'):
        for t in (ptab, ftab):
            if c in t.colnames: t[c].unit = 'deg'
    for t in (ptab, ftab):
        for c in t.colnames:
            if c.startswith(('dRA', 'dDEC', 'E_', 'Sigma_', 'Res', 'Scatter', 'Offset')):
                t[c].unit = 'arcsec'
    ptab['dCRVAL1'].unit = ptab['dCRVAL2'].unit = 'deg'

    hdr = fits.Header()
    hdr['COMMENT'] = 'Offsets are reference (LoTSS-DR3) - LoLSS, arcsec on the sky (dRA includes cos(Dec)).'
    hdr['COMMENT'] = 'To correct a pointing add dCRVAL1/dCRVAL2 to CRVAL1/CRVAL2 of its image.'
    hdr['COMMENT'] = 'Facet flags: 1 few sources, 2 offset from the pointing, 4 large scatter, 8 outside the polygons.'
    hdr['COMMENT'] = 'Pointing flags: 1 few sources, 2 outlier wrt the survey, 4 no facet file, 8 large scatter.'
    for k, val in (('SNRCUT', args.snr_cut), ('SYSERR', args.sys_err), ('FNMIN', args.facet_nmin), ('FNMINSC', args.facet_nmin_scatter), ('FMAXOFF', args.facet_max_offset),
                   ('FNSIG', args.facet_nsigma), ('FMAXSCAT', args.facet_max_scatter), ('PNMIN', args.pointing_nmin),
                   ('PMAXOFF', args.pointing_max_offset), ('PMAXSCAT', args.pointing_max_scatter), ('PNSIG', args.pointing_nsigma),
                   ('SURVDRA', float(s_ra)), ('SURVDDEC', float(s_dec))):
        hdr[k] = val
    outfile = os.path.join(args.outdir, 'astrometry_pointings.fits')
    hp = fits.BinTableHDU(ptab, name='POINTINGS')
    hp.header.update(hdr)
    hf = fits.BinTableHDU(ftab, name='FACETS')
    hf.header.update(hdr)
    fits.HDUList([fits.PrimaryHDU(header=hdr), hp, hf]).writeto(outfile, overwrite=True)

    # summary
    real = ftab['Facet'] >= 0
    meas = real & ((ftab['Flag'] & FLAG_FEW) == 0)
    print('Pointings: %i, offset survey median dRA %.2f" dDEC %.2f"; few sources: %i; large scatter: %i; outliers: %i'
          % (len(ptab), s_ra, s_dec, np.sum((ptab['Flag'] & PFLAG_FEW) > 0),
             np.sum((ptab['Flag'] & PFLAG_BAD) == PFLAG_SCATTER), np.sum(out)))
    if ok.any():
        print('Pointing offsets |off|: median %.2f", 90%% %.2f", max %.2f"' % tuple(np.percentile(ptab['Offset'][ok], [50, 90, 100])))
    print('Facets: %i, measured (>= %i sources): %i, OFFSET: %i, SCATTER: %i'
          % (real.sum(), args.facet_nmin, meas.sum(), np.sum(meas & ((ftab['Flag'] & FLAG_OFFSET) > 0)),
             np.sum(meas & ((ftab['Flag'] & FLAG_SCATTER) > 0))))
    print('Written %s' % outfile)
    plot_pointings(ptab, args)
    plot_facets(ftab, args)


def plot_pointings(ptab, args):
    ok = (ptab['Flag'] & PFLAG_BAD) == 0
    out = (ptab['Flag'] & PFLAG_OUTLIER) > 0
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.5), gridspec_kw=dict(width_ratios=[1, 1.6, 1]))

    # symmetric log axes: linear within +-LT (where the bulk and the outlier threshold are), log outside
    LT = max(1., args.pointing_max_offset)
    sca = (ptab['Flag'] & PFLAG_BAD) == PFLAG_SCATTER
    s_ra, s_dec = (np.median(ptab['dRA'][ok]), np.median(ptab['dDEC'][ok])) if ok.any() else (0., 0.)
    shown = ok | sca
    lim = max(2., np.nanmax(np.abs(np.r_[ptab['dRA'][shown], ptab['dDEC'][shown]]))*1.3) if shown.any() else 2.

    def symlog_axes(ax, x=True, y=True):
        if x: ax.set_xscale('symlog', linthresh=LT, linscale=1.); ax.set_xlim(-lim, lim)
        if y: ax.set_yscale('symlog', linthresh=LT, linscale=1.); ax.set_ylim(-lim, lim)

    ax = axs[0]
    ax.axvspan(-LT, LT, color='0.94', zorder=0, lw=0)
    ax.axhspan(-LT, LT, color='0.94', zorder=0, lw=0)
    ax.scatter(ptab['dRA'][ok & ~out], ptab['dDEC'][ok & ~out], s=12, color='k', alpha=0.2, lw=0,
               label='pointings (%i)' % np.sum(ok & ~out))
    ax.scatter(ptab['dRA'][out], ptab['dDEC'][out], s=16, color='r', alpha=0.4, lw=0, label='outliers (%i)' % np.sum(out))
    ax.scatter(ptab['dRA'][sca], ptab['dDEC'][sca], s=16, marker='x', color='tab:blue', alpha=0.5, lw=0.8,
               label='large scatter, not corrected (%i)' % np.sum(sca))
    for t in ptab[out][np.argsort(-np.array(ptab['Offset'][out]))][:10]:   # names of the 10 largest outliers
        ax.annotate(t['Pointing'], (t['dRA'], t['dDEC']), fontsize=6, color='r', alpha=0.8)
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot(s_ra + args.pointing_max_offset*np.cos(th), s_dec + args.pointing_max_offset*np.sin(th), 'r--', lw=0.8)
    ax.axhline(0, color='grey', lw=0.5); ax.axvline(0, color='grey', lw=0.5)
    symlog_axes(ax)
    ax.set_box_aspect(1)
    ax.set_xlabel('dRA [arcsec]'); ax.set_ylabel('dDec [arcsec]')
    ax.set_title('Pointing offsets (DR3 - LoLSS)\nlinear within $\\pm$%g" (shaded), log outside' % LT, fontsize=10)
    ax.legend(fontsize=7, loc='upper left')

    ax = axs[1]
    ra = np.array(ptab['RA']); ra = np.where(ra > 180, ra - 360, ra)
    off = np.array(ptab['Offset'])
    sc = ax.scatter(ra[ok], ptab['DEC'][ok], c=off[ok], cmap='viridis', s=12, vmin=0, vmax=np.nanpercentile(off[ok], 98) if ok.any() else 1)
    # arrows drawn in screen units (same length per arcsec in every direction): the axes are not equal,
    # 360 deg of RA are squeezed in the same width as ~70 deg of Dec, so data-unit arrows would hide dRA.
    # Screen x grows to the right, where RA decreases on the inverted axis: hence -dRA
    ARROW_IN = 0.05   # inches per arcsec
    q = ax.quiver(ra[ok], ptab['DEC'][ok], -ptab['dRA'][ok], ptab['dDEC'][ok], angles='uv', scale_units='inches',
                  scale=1./ARROW_IN, width=0.0015, color='k', alpha=0.6)
    ax.quiverkey(q, 0.9, 1.03, 5., '5"', labelpos='E', coordinates='axes', fontproperties=dict(size=8))
    ax.scatter(ra[out], ptab['DEC'][out], s=40, facecolors='none', edgecolors='r')
    few = (ptab['Flag'] & PFLAG_FEW) > 0
    ax.scatter(ra[few], ptab['DEC'][few], s=12, marker='x', color='grey', label='too few sources')
    ax.scatter(ra[sca], ptab['DEC'][sca], s=16, marker='x', color='tab:blue', label='large scatter')
    ax.invert_xaxis()
    ax.set_xlabel('RA [deg]'); ax.set_ylabel('Dec [deg]')
    ax.set_title('Offsets on the sky (arrow = shift to apply, red: outliers)', loc='left')
    plt.colorbar(sc, ax=ax, label='|offset| [arcsec]')
    ax.legend(fontsize=8, loc='lower left')

    ax = axs[2]
    # bins uniform on the symlog axis
    lin = np.arange(-LT, LT + 1e-9, LT/10)
    out_b = np.geomspace(LT, lim, 25)[1:]
    bins = np.r_[-out_b[::-1], lin, out_b]
    ax.hist(np.clip(ptab['dRA'][ok], -lim, lim), bins=bins, histtype='step', label='dRA')
    ax.hist(np.clip(ptab['dDEC'][ok], -lim, lim), bins=bins, histtype='step', label='dDec')
    ax.axvspan(-LT, LT, color='0.94', zorder=0, lw=0)
    symlog_axes(ax, y=False)
    ax.set_yscale('log')
    ax.set_xlabel('offset [arcsec]'); ax.set_ylabel('pointings')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, 'astrometry_pointings.pdf'), dpi=150)
    plt.close(fig)


def plot_facets(ftab, args):
    real = ftab['Facet'] >= 0
    meas = real & ((ftab['Flag'] & FLAG_FEW) == 0)
    boff = meas & ((ftab['Flag'] & FLAG_OFFSET) > 0)
    bsca = meas & ((ftab['Flag'] & FLAG_SCATTER) > 0)
    good = meas & ~boff & ~bsca
    fig, axs = plt.subplots(2, 2, figsize=(12, 11))

    ax = axs[0, 0]
    ax.scatter(ftab['Res_dRA'][good], ftab['Res_dDEC'][good], s=3, color='k', alpha=0.3, label='good (%i)' % good.sum())
    ax.scatter(ftab['Res_dRA'][bsca & ~boff], ftab['Res_dDEC'][bsca & ~boff], s=8, color='orange', label='SCATTER (%i)' % (bsca & ~boff).sum())
    ax.scatter(ftab['Res_dRA'][boff], ftab['Res_dDEC'][boff], s=8, color='r', label='OFFSET (%i)' % boff.sum())
    lim = max(3., np.nanpercentile(np.abs(np.r_[ftab['Res_dRA'][meas], ftab['Res_dDEC'][meas]]), 99.5)*1.1) if meas.any() else 3.
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.add_artist(plt.Circle((0, 0), args.facet_max_offset, fill=False, ls='--', color='r', lw=0.8))
    ax.set_xlabel('facet - pointing dRA [arcsec]'); ax.set_ylabel('facet - pointing dDec [arcsec]')
    ax.set_title('Facet offsets relative to their pointing')
    ax.legend(fontsize=8)

    ax = axs[0, 1]
    ax.scatter(ftab['Dist'][good], ftab['Res'][good], s=3, color='k', alpha=0.3)
    ax.scatter(ftab['Dist'][bsca & ~boff], ftab['Res'][bsca & ~boff], s=8, color='orange')
    ax.scatter(ftab['Dist'][boff], ftab['Res'][boff], s=8, color='r')
    ax.axhline(args.facet_max_offset, color='r', ls='--', lw=0.8)
    ax.set_xlabel('facet distance from the pointing centre [deg]'); ax.set_ylabel('|facet - pointing| [arcsec]')
    ax.set_ylim(0, max(lim, np.nanmax(ftab['Res'][meas])*1.05) if meas.any() else lim)

    ax = axs[1, 0]
    bins = np.linspace(0, max(args.facet_max_scatter*2, np.nanpercentile(ftab['Scatter'][meas], 99.5) if meas.any() else 1), 60)
    ax.hist(ftab['Scatter'][meas], bins=bins, color='k', histtype='step')
    ax.axvline(args.facet_max_scatter, color='orange', ls='--')
    ax.set_xlabel('robust scatter of the source offsets in a facet [arcsec]'); ax.set_ylabel('facets')
    ax.set_yscale('log')

    ax = axs[1, 1]
    fra = np.array(ftab['RA']); fra = np.where(fra > 180, fra - 360, fra)
    ax.scatter(fra[good], ftab['DEC'][good], s=1, color='lightgrey')
    ax.scatter(fra[bsca & ~boff], ftab['DEC'][bsca & ~boff], s=6, color='orange')
    ax.scatter(fra[boff], ftab['DEC'][boff], s=6, color='r')
    ax.invert_xaxis()
    ax.set_xlabel('RA [deg]'); ax.set_ylabel('Dec [deg]')
    ax.set_title('Flagged facets on the sky')
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, 'astrometry_facets.pdf'), dpi=150)
    plt.close(fig)


def run_pointings(args):
    global DR3
    os.makedirs(os.path.join(args.outdir, 'pointings'), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, 'catalogues'), exist_ok=True)
    if not args.analyse_only:
        pointings = args.only if args.only else sorted(p for p in os.listdir(args.dir_done)
                                                         if os.path.isdir(os.path.join(args.dir_done, p)))

        def needs_match(p):
            match, cat = os.path.join(args.outdir, 'pointings', p + '.fits'), catalogue_path(p, args)
            return args.redo or args.redo_bdsf or not os.path.exists(match) or not os.path.exists(cat) \
                or os.path.getmtime(cat) > os.path.getmtime(match)

        todo = [p for p in pointings if needs_match(p)]
        n_bdsf = sum(1 for p in todo if args.redo_bdsf or not os.path.exists(catalogue_path(p, args)))
        print('%i pointings, %i to match, of which %i need PyBDSF' % (len(pointings), len(todo), n_bdsf), flush=True)
        if n_bdsf:
            import bdsf  # noqa: F401  imported before the workers are forked
        if todo:
            DR3 = load_dr3(args.dr3_cache, args.dr3)
            items = [(p, args) for p in todo]
            if args.njobs > 1:
                # ProcessPoolExecutor and not multiprocessing.Pool: PyBDSF starts processes of its own,
                # which daemonic Pool workers cannot do. fork: the workers share the DR3 arrays loaded above
                from concurrent.futures import ProcessPoolExecutor
                with ProcessPoolExecutor(max_workers=args.njobs, mp_context=multiprocessing.get_context('fork')) as ex:
                    res = list(ex.map(_measure, items))
            else:
                res = [_measure(i) for i in items]
            print('%i/%i pointings measured' % (sum(r is not None for r in res), len(todo)), flush=True)
    analyse(args)


def main():
    class _Fmt(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter): pass
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=_Fmt)
    p = parser
    p.add_argument('--outdir', default=os.path.join(cfg.DIR_ANALYSIS, 'astrometry', ''), help='Output directory')
    p.add_argument('--dir-done', default=cfg.DIR_DONE, help='Directory with the pointings (only read)')
    p.add_argument('--file-image', default=cfg.FILE_IMAGE, help='Image of each pointing')
    p.add_argument('--file-facets', default='facetsS-c0.reg', help='Facet regions of each pointing')
    p.add_argument('--dr3', default=common.FILE_LOTSS_CAT, help='LoTSS-DR3 source catalogue')
    p.add_argument('--dr3-cache', default=None, help='Cache of the DR3 columns (default: <outdir>/dr3_cache.npz)')
    p.add_argument('--only', nargs='+', default=None, metavar='POINTING', help='Only these pointings (e.g. P000+31s)')
    p.add_argument('-j', '--njobs', type=int, default=1, help='Pointings measured in parallel')
    p.add_argument('--redo', action='store_true', help='Match again the pointings already in the cache (PyBDSF catalogues are kept)')
    p.add_argument('--redo-bdsf', action='store_true', help='Run PyBDSF again also on the pointings that already have a catalogue')
    p.add_argument('--bdsf-ncores', type=int, default=1, help='Cores used by each PyBDSF run (0: PyBDSF default, all of them)')
    p.add_argument('--analyse-only', action='store_true', help='Only redo the analysis and the plots from the cache')
    p.add_argument('--beam', type=float, default=cfg.RESOLUTION, help='Beam FWHM [arcsec] if the image has no BMAJ')
    p.add_argument('--beamcut', type=float, default=cfg.BEAM_THRESHOLD, help='Primary beam level below which sources are not used')
    g = p.add_argument_group('selection and match (need --redo to change)')
    g.add_argument('--ref-maj-max', type=float, default=20., help='Max DR3 major axis [arcsec]')
    g.add_argument('--ref-flux-min', type=float, default=3., help='Min DR3 total flux [mJy]')
    g.add_argument('--ref-iso-radius', type=float, default=45., help='DR3 isolation radius [arcsec]')
    g.add_argument('--ref-iso-frac', type=float, default=0.2, help='No DR3 source brighter than this fraction of the reference within the isolation radius')
    g.add_argument('--lolss-iso-radius', type=float, default=45., help='No other LoLSS source within this radius [arcsec]')
    g.add_argument('--lolss-codes', default='SM', help='PyBDSF S_Code of the LoLSS sources used (S single, M multi-Gaussian, C single in an island with others)')
    g.add_argument('--lolss-maj-max', type=float, default=2.5, help='Max LoLSS major axis, in beams')
    g.add_argument('--snr-min', type=float, default=5., help='Min S/N (Peak_flux/Isl_rms) of a LoLSS source')
    g.add_argument('--match-radius', type=float, default=15., help='Max LoLSS-DR3 distance [arcsec]')
    g = p.add_argument_group('flags (analysis: can be changed with --analyse-only)')
    g.add_argument('--snr-cut', type=float, default=5., help='Min S/N of the sources used in the analysis')
    g.add_argument('--sys-err', type=float, default=1., help='Systematic error added in quadrature to the PyBDSF position error in the weights [arcsec]')
    g.add_argument('--facet-nmin', type=int, default=3, help='Min sources to judge the offset of a facet')
    g.add_argument('--facet-nmin-scatter', type=int, default=5, help='Min sources to judge the scatter of a facet')
    g.add_argument('--facet-max-offset', type=float, default=5., help='Facet flagged OFFSET above this |facet - pointing| [arcsec]...')
    g.add_argument('--facet-nsigma', type=float, default=3., help='...if also significant at this level')
    g.add_argument('--facet-max-scatter', type=float, default=5., help='Facet flagged SCATTER above this robust scatter [arcsec]')
    g.add_argument('--pointing-nmin', type=int, default=20, help='Min sources for a pointing offset')
    g.add_argument('--pointing-max-scatter', type=float, default=4., help='Pointing flagged SCATTER (offset not used) above this robust scatter of its sources [arcsec]')
    g.add_argument('--pointing-max-offset', type=float, default=1., help='Pointing flagged OUTLIER above this distance from the survey median [arcsec]...')
    g.add_argument('--pointing-nsigma', type=float, default=5., help='...if also significant at this level')

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    if args.dr3_cache is None:
        args.dr3_cache = os.path.join(args.outdir, 'dr3_cache.npz')
    run_pointings(args)


if __name__ == '__main__':
    main()
