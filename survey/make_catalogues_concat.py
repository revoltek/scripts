#!/usr/bin/env python3

import os,sys, glob, time
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astropy.coordinates import match_coordinates_sky

# Code to concatenate catalogs made from mosaics, remove repeated sources and to manipulate catalog entries to the final catalog format.

## To do
# -- Need to finalise the column choices
# -- We probably want to check the masks for each image use those mask islands to determine which mosaic to select the source from as perhaps there  are very extneded sources approximately midway between pointings (especially for gaussian components list)
# -- We probably want an entry in the catalog to say if the source was in mask used for the final deconvolution
# -- Needs speeding up -- each pointing separately
# -- update the errors that are put on the RA and DEC (could use astrometry error maps?)
# -- update the source extension definition

deg2arcsec = 3600

def replace_to(t, col, dtype):
    float_col = t[col].astype(dtype)
    t.replace_column(col, float_col)

def concat_catalogs(cats, outconcatcat):
    # Use the first catalogue as a dummy and then just update the entries

    concat_table = vstack([ Table.read(c) for c in cats ], join_type='exact')
    concat_table.sort('Source_name')
    concat_table.write(outconcatcat, overwrite=True)


def filter_catalogs(t_full, t_this, srl_outname, gaus_outname):

    cat = Table.read(t_this['catalogue'])
    cat_gaus = Table.read(t_this['catalogue_gaus'])

    print('Matching (%s)...' % t_this['fieldname'])
    idx_match, sep, _ = match_coordinates_sky(SkyCoord(cat['RA'], cat['DEC']), SkyCoord(t_full['RA'], t_full['DEC']))
    # find all idx_match with the same index of the current fieldname and keep only those entry of the catalogue
    initial_len = len(cat)
    cat = cat[idx_match == list(t_full['fieldname']).index(t_this['fieldname'])]
    print("Cat: %i -> %i" % (initial_len, len(cat)))
    # select only gaussians that are part of selected sources
    initial_len = len(cat_gaus)
    idx = [i for i in range(len(cat_gaus)) if cat_gaus['Source_id'][i] in cat['Source_id']]
    cat_gaus = cat_gaus[idx]
    print("Gaus: %i -> %i" % (initial_len, len(cat_gaus)))

    # name
    sc = SkyCoord(cat['RA'], cat['DEC'], frame='icrs')
    identityRA = sc.to_string(style='hmsdms',sep='',precision=1)
    identityDEC = sc.to_string(style='hmsdms',sep='',precision=0)
    identity = [ 'LOL1J'+i[:8]+j[7:] for (i,j) in zip(identityRA,identityDEC) ]
    cat['Source_name'] = identity

    identity = []
    for idx in cat_gaus['Source_id']:
        identity.append( cat['Source_name'][ cat['Source_id'] == idx ][0] )
    cat_gaus['Source_name'] = identity

    # add mosaic id
    cat['Mosaic_id'] = t_this['fieldname']
    cat_gaus['Mosaic_id'] = t_this['fieldname']

    cat.remove_columns(['Source_id', 'Isl_id',
                  'RA_max', 'E_RA_max', 'DEC_max', 'E_DEC_max',
                  'Maj_img_plane','E_Maj_img_plane',
                  'Min_img_plane','E_Min_img_plane',
                  'PA_img_plane','E_PA_img_plane',
                  'DC_Maj_img_plane','E_DC_Maj_img_plane',
                  'DC_Min_img_plane','E_DC_Min_img_plane',
                  'DC_PA_img_plane','E_DC_PA_img_plane',
                  'Isl_Total_flux','E_Isl_Total_flux',
                  'Isl_mean','Resid_Isl_rms','Resid_Isl_mean'])

    cat_gaus.remove_columns(['Source_id', 'Isl_id', 'Wave_id',
                  'Xposn','E_Xposn','Yposn','E_Yposn',
                  'Maj_img_plane','E_Maj_img_plane',
                  'Min_img_plane','E_Min_img_plane',
                  'PA_img_plane','E_PA_img_plane',
                  'DC_Maj_img_plane','E_DC_Maj_img_plane',
                  'DC_Min_img_plane','E_DC_Min_img_plane',
                  'DC_PA_img_plane','E_DC_PA_img_plane',
                  'Isl_Total_flux','E_Isl_Total_flux',
                  'Isl_mean','Resid_Isl_rms','Resid_Isl_mean'])

#    print(cat.columns)
#    print(cat_gaus.columns)

    for t in [cat, cat_gaus]:

        # values in mJy
        t['Peak_flux'] *= 1e3; t['Peak_flux'].unit = 'mJy/beam'
        t['E_Peak_flux'] *= 1e3; t['E_Peak_flux'].unit = 'mJy/beam'
        t['Total_flux'] *= 1e3; t['Total_flux'].unit = 'mJy'
        t['E_Total_flux'] *= 1e3; t['E_Total_flux'].unit = 'mJy'
        t['Isl_rms'] *= 1e3; t['Isl_rms'].unit = 'mJy/beam'
        
        # values in arcsec
        t['E_RA'] *= deg2arcsec; t['E_RA'].unit = 'arcsec'
        t['E_DEC'] *= deg2arcsec; t['E_DEC'].unit = 'arcsec'
        
        t['Maj'] *= deg2arcsec; t['Maj'].unit = 'arcsec'
        t['E_Maj'] *= deg2arcsec; t['E_Maj'].unit = 'arcsec'
        t['Min'] *= deg2arcsec; t['Min'].unit = 'arcsec'
        t['E_Min'] *= deg2arcsec; t['E_Min'].unit = 'arcsec'
        
        t['DC_Maj'] *= deg2arcsec; t['DC_Maj'].unit = 'arcsec'
        t['E_DC_Maj'] *= deg2arcsec; t['E_DC_Maj'].unit = 'arcsec'
        t['DC_Min'] *= deg2arcsec; t['DC_Min'].unit = 'arcsec'
        t['E_DC_Min'] *= deg2arcsec; t['E_DC_Min'].unit = 'arcsec'
        
        for col in ['E_RA','E_DEC','Total_flux','E_Total_flux','Peak_flux','E_Peak_flux','Maj','E_Maj','Min','E_Min','PA','E_PA','DC_Maj','E_DC_Maj','DC_Min','E_DC_Min','DC_PA','E_DC_PA','Isl_rms']:
            replace_to(t, col, 'f4') # to float from double

    # reorder columns
    new_order = ['Source_name','RA','E_RA','DEC','E_DEC','Total_flux','E_Total_flux','Peak_flux','E_Peak_flux','Maj','E_Maj','Min','E_Min','PA','E_PA','DC_Maj','E_DC_Maj','DC_Min','E_DC_Min','DC_PA','E_DC_PA','Isl_rms','S_Code', 'Mosaic_id']
    cat = cat[new_order]
    new_order = ['Source_name','Gaus_id','RA','E_RA','DEC','E_DEC','Total_flux','E_Total_flux','Peak_flux','E_Peak_flux','Maj','E_Maj','Min','E_Min','PA','E_PA','DC_Maj','E_DC_Maj','DC_Min','E_DC_Min','DC_PA','E_DC_PA','Isl_rms','S_Code', 'Mosaic_id']
    cat_gaus = cat_gaus[new_order]

    with open(srl_outname.replace('.fits','.reg'), 'w') as regionfile:
        regionfile.write('# Region file format: DS9 version 4.0\n')
        regionfile.write('global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
        regionfile.write('fk5\n')
        for c in cat:
            if not np.isnan(c['Maj']):
                regionfile.write('ellipse(%s,%s,%s,%s,%s)\n' % (c['RA'],c['DEC'],c['Maj']/3600,c['Min']/3600,c['PA']+90))
            else:
                regionfile.write('box(%s,%s,5.0",5.0",0.0)\n' % (c['RA'],c['DEC']))

    with open(gaus_outname.replace('.fits','.reg'), 'w') as regionfile:
        regionfile.write('# Region file format: DS9 version 4.0\n')
        regionfile.write('global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
        regionfile.write('fk5\n')
        for c in cat_gaus:
            if not np.isnan(c['Maj']):
                regionfile.write('ellipse(%s,%s,%s,%s,%s)\n' % (c['RA'],c['DEC'],c['Maj']/3600,c['Min']/3600,c['PA']+90))
            else:
                regionfile.write('box(%s,%s,5.0",5.0",0.0)\n' % (c['RA'],c['DEC']))

    cat.write(srl_outname, overwrite=True)
    cat_gaus.write(gaus_outname, overwrite=True)

def do_concat(mosdir):

    mosaiccats = sorted(glob.glob('%s/catalogues/*cat.fits' % mosdir))
    gauscats = [x.replace('cat.fits','gaus.fits') for x in mosaiccats]

    # Determine all pointing coordinates
    fieldnames = []
    ras = []
    decs = []
    for mosaiccat, gauscat in zip(mosaiccats, gauscats):
        mosaicfile = mosaiccat.replace('.cat.fits', '.fits').replace('catalogues/', '')
        fieldnames.append(mosaicfile[9:16])
        print('Adding catalogue: %s (field: %s)' % (mosaiccat, fieldnames[-1]))
        with fits.open(mosaicfile) as f:
            ras.append(f[0].header['CRVAL1'])
            decs.append(f[0].header['CRVAL2'])
    
    t_mosaics = Table({'fieldname':fieldnames, 'catalogue':mosaiccats, 'catalogue_gaus':gauscats, 'RA':ras, 'DEC':decs}, units={'RA':u.deg,'DEC':u.deg})
 
    srlcatnames = []
    gauscatnames = []

    for t_mosaic in t_mosaics:
        print('### Working on %s' % t_mosaic['fieldname'])
        srl_outname = t_mosaic['catalogue'].replace('cat.fits','cat-filt.fits')
        gaus_outname = t_mosaic['catalogue_gaus'].replace('gaus.fits','gaus-filt.fits')
        if not os.path.exists(srl_outname) or not os.path.exists(gaus_outname):
            filter_catalogs(t_mosaics, t_mosaic, srl_outname, gaus_outname)
            
        srlcatnames.append(srl_outname)
        gauscatnames.append(gaus_outname)

    print('### Concatenating %s files' % len(srlcatnames))
    concat_catalogs(srlcatnames,'LoLSS_DR1_rolling.srl.fits')
    concat_catalogs(gauscatnames,'LoLSS_DR1_rolling.gaus.fits')

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Concatenate ddf-pipeline mosaic directories')
    parser.add_argument('--mosdir', type=str, help='mosaic directory name')
    args = parser.parse_args()

    do_concat(args.mosdir)

# call as e.g. make_catalogues_concat.py --mosdir=mosaic-i
