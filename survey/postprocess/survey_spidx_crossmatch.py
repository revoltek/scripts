#!/usr/bin/env python3
"""
survey_spidx_crossmatch.py  (from imagemakers/crossmatch.ipynb)

Cross-match the LoLSS catalogue with LoTSS and NVSS and compute the spectral
indices, summing the components of the higher-resolution catalogue that fall
within the matching radius of each lower-resolution source:

  cat_lolss-lotss.fits : LoLSS sources + summed LoTSS flux within 15" -> spidxLoLSSLoTSS
  cat_lolss-nvss.fits  : NVSS sources + summed LoLSS flux within 45"  -> spidxLoLSSNVSS

Outputs in --outdir.
"""

import os, argparse
import numpy as np
from astropy.table import Table

import lolss_analysis_common as common
from lib_analysis_lolss import cross_match_multicomp

FREQ_LOLSS, FREQ_LOTSS, FREQ_NVSS = 54., 144., 1400.  # MHz


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cat', default=common.FILE_LOLSS_CAT, help='LoLSS source catalogue (default: %(default)s)')
    parser.add_argument('--lotss', default=common.FILE_LOTSS_CAT, help='LoTSS source catalogue (default: %(default)s)')
    parser.add_argument('--nvss', default=common.FILE_NVSS_CAT, help='NVSS catalogue (default: %(default)s)')
    parser.add_argument('--compactness-fit', default=common.FILE_COMPACTNESS_FIT,
                        help='Output of survey_compactness.py (default: %(default)s)')
    parser.add_argument('--radius-lotss', type=float, default=15., help='LoLSS-LoTSS matching radius [arcsec] (default: 15)')
    parser.add_argument('--radius-nvss', type=float, default=45., help='NVSS-LoLSS matching radius [arcsec] (default: 45)')
    parser.add_argument('--outdir', default=common.DIR_EXTCAT, help='Output directory (default: %(default)s)')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cat_lolss = Table.read(args.cat)
    cat_lolss.rename_column('Total_flux', 'flux_lolss')
    cat_lolss.rename_column('Peak_flux', 'peak_lolss')
    cat_lolss.rename_column('E_Total_flux', 'flux_lolss_e')

    # LoLSS x LoTSS
    cat_lotss = Table.read(args.lotss)
    cat_lotss.rename_column('Total_flux', 'flux_lotss')
    cat_lotss.rename_column('Peak_flux', 'peak_lotss')
    cat_lotss.rename_column('E_Total_flux', 'flux_lotss_e')
    cat_lotss['is_pnt'] = common.is_compact(cat_lotss['flux_lotss'], cat_lotss['peak_lotss'],
                                            cat_lotss['flux_lotss_e'], common.COMPACT_LOTSS)
    print('Is pnt (LoTSS) %f%%' % (100*np.sum(cat_lotss['is_pnt'])/len(cat_lotss)))

    cat_lolss_lotss = cross_match_multicomp(cat_lolss, cat_lotss, dist=args.radius_lotss, flux_col='flux_lotss',
                                            columns_to_keep=['is_pnt'])
    del cat_lotss
    print('Match LoTSS: %i' % np.sum(np.isfinite(cat_lolss_lotss['flux_lotss_sum'])))
    print('un-Match LoTSS: %i' % np.sum(np.isnan(cat_lolss_lotss['flux_lotss_sum'])))
    cat_lolss_lotss['spidxLoLSSLoTSS'] = np.log10(cat_lolss_lotss['flux_lotss_sum']/cat_lolss_lotss['flux_lolss']) \
        / np.log10(FREQ_LOTSS/FREQ_LOLSS)
    cat_lolss_lotss.write(os.path.join(args.outdir, 'cat_lolss-lotss.fits'), overwrite=True)

    # NVSS x LoLSS
    cat_lolss['is_pnt'] = common.is_compact(cat_lolss['flux_lolss'], cat_lolss['peak_lolss'], cat_lolss['flux_lolss_e'],
                                            common.load_compact_lolss(args.compactness_fit))
    print('Is pnt (LoLSS) %f%%' % (100*np.sum(cat_lolss['is_pnt'])/len(cat_lolss)))

    cat_nvss = Table.read(args.nvss)
    cat_nvss.rename_column('Total_flux', 'flux_nvss')
    cat_nvss['flux_nvss'] = cat_nvss['flux_nvss']/1e-3  # Jy -> mJy

    cat_lolss_nvss = cross_match_multicomp(cat_nvss, cat_lolss, dist=args.radius_nvss, flux_col='flux_lolss',
                                           columns_to_keep=['is_pnt'])
    print('Match NVSS: %i' % np.sum(np.isfinite(cat_lolss_nvss['flux_lolss_sum'])))
    print('un-Match NVSS: %i' % np.sum(np.isnan(cat_lolss_nvss['flux_lolss_sum'])))
    cat_lolss_nvss['spidxLoLSSNVSS'] = np.log10(cat_lolss_nvss['flux_nvss']/cat_lolss_nvss['flux_lolss_sum']) \
        / np.log10(FREQ_NVSS/FREQ_LOLSS)
    cat_lolss_nvss.write(os.path.join(args.outdir, 'cat_lolss-nvss.fits'), overwrite=True)


if __name__ == '__main__':
    main()
