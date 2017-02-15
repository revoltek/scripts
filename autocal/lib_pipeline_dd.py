#!/usr/bin/python

import os, sys
import logging
import numpy as np
from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord, match_coordinates_sky
from astropy.io import fits as pyfits
from astropy import wcs
import astropy.units as u
import pyregion
from pyregion.parser_helper import Shape
from autocal.lib_pipeline_img import *
try:
    from scipy.spatial import Voronoi
except:
    logging.error("Load latest scipy with 'use Pythonlibs'")
    sys.exit(1)

def table_to_circ_region(table, outfile, racol='RA', deccol='DEC', sizecol='dd_size', color='red'):
    """
    Get a table with ra, dec, size and generate a circular ds9 region 
    """

    color = [color] * len(table)

    regions = []
    for i, r in enumerate(table):
        s = Shape('circle', None)
        s.coord_format = 'fk5'
        s.coord_list = [ r[racol], r[deccol], r[sizecol] ] # ra, dec, radius
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '5', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        s.comment = 'color={} text="{}"'.format(color[i], str(i))
        regions.append(s)

    regions = pyregion.ShapeList(regions)
    regions.write(outfile)


def make_directions_from_skymodel(filename, outdir='regions/', flux_min_Jy=1.0, size_max_arcmin=3.0,
    directions_separation_max_arcmin=5.0, directions_max_num=20, flux_min_for_merging_Jy=0.2):
    """
    Selects appropriate calibrators from srl file
    Parameters
    ----------
    filename : srl file
        Skymodel made by grouping clean components of dir-independent model
    outdir: string
        Directory where to save the ds9 regions
    flux_min_Jy : float
        Minimum flux density for a calibrator in Jy
    size_max_arcmin : float
        Maximum size for a calibrator in arcmin
    directions_separation_max_arcmin : float
        Maximum separation in arcmin between two calibrators for gouping into a
        single direction
    directions_max_num : int, optional
        Limit total number of directions to this value
    flux_min_for_merging_Jy : float, optional
        Minimum flux density for a source to be considered for merging
    Returns
    -------
    table : astropy table
        Table with direction information
    """
    # open astropy table
    t = Table.read(filename, format='fits')['RA','DEC','Maj','Peak_flux','Total_flux'] # restrict to some cols
    t.rename_column('Maj', 'dd_size') # use Maj as proxy for region size
    logging.info('# sources initial: %i' % len(t))

    # exclude sources that are too extended
    t_large = t[ (t['dd_size'] >= size_max_arcmin*u.arcmin) ]
    t = t[ (t['dd_size'] < size_max_arcmin*u.arcmin) ]
    logging.info('# sources after cut on size: %i' % len(t))
    logging.info('# large sources: %i' % len(t_large))
    if len(t) == 0:
        logging.critical("No sources found that meet the specified max size criterion.")
        sys.exit(1)

    # exclude sources that are too faint
    t_large = t_large[ (t_large['Total_flux'] > flux_min_for_merging_Jy) ]
    t = t[ (t['Total_flux'] > flux_min_for_merging_Jy) ]
    logging.info('# sources after cut min flux for merging: %i' % len(t))
    if len(t) == 0:
        logging.critical("No sources found above %f Jy." % flux_min_for_merging_Jy )
        sys.exit(1)

    t.sort('Total_flux')
    t.reverse()

    # combine nearby sources
    for s in t:
        # if ra/dec changes, continue finding nearby sources until no-sources are found
        updated = True
        while updated:
            dists = SkyCoord(ra=s['RA']*u.degree, dec=s['DEC']*u.degree).separation(SkyCoord(ra=t['RA'], dec=t['DEC']))
            updated = False
            for i, dist in enumerate(dists):
                if dist < directions_separation_max_arcmin*u.arcmin and dist > 0.*u.degree:
                    # if a source is dominant keep that at the center of the patch
                    if t['Peak_flux'][i] > 3*s['Peak_flux']:
                        s['RA'] = t['RA'][i]
                        s['DEC'] = t['DEC'][i]
                        updated = True
                    # other wise weighted mean
                    elif t['Peak_flux'][i] > s['Peak_flux']:
                        s['RA'] = (s['RA']*s['Peak_flux'] + t['RA'][i]*t['Peak_flux'][i])/(s['Peak_flux']+t['Peak_flux'][i])
                        s['DEC'] = (s['DEC']*s['Peak_flux'] + t['DEC'][i]*t['Peak_flux'][i])/(s['Peak_flux']+t['Peak_flux'][i])
                        updated = True

                    s['dd_size'] = max(s['dd_size'], t['dd_size'][i]) + dist.degree

                    s['Total_flux'] += t['Total_flux'][i]
                    s['Peak_flux'] = max(s['Peak_flux'], t['Peak_flux'][i])

                    t.remove_rows(i)
    logging.info('# sources after combining close-by sources: %i' % len(t))

    t.sort('Total_flux')
    t.reverse()

    # Filter patches on total flux density limit
    t = t[ (t['Total_flux'] > flux_min_Jy) ]
    logging.info('# sources after cut min flux: %i' % len(t))
    if len(t) == 0:
        logging.critical("No sources or merged groups found that meet the specified "
            "min total flux density criterion.")
        sys.exit(1)

    # Trim directions list to get directions_max_num of directions
    if directions_max_num is not None:
        t = t[:directions_max_num]
        logging.info('# sources after cut on max directions: %i' % len(t))

    for s in t_large:
        dists = SkyCoord(ra=s['RA']*u.degree, dec=s['DEC']*u.degree).separation(SkyCoord(ra=t['RA'], dec=t['DEC']))
        for i, dist in enumerate(dists):
            if dist < directions_separation_max_arcmin*u.arcmin:
                t['Total_flux'][i] += s['Total_flux']
                t['dd_size'][i] = max(s['dd_size'], t['dd_size'][i]) + dist.degree

    t.sort('Total_flux')
    t.reverse()

    # Writedd global region files
    table_to_circ_region(t, outdir+'/all.reg')

    # save source by source
    for i, s in enumerate(t):
        table_to_circ_region(t[i:i+1], outdir+'/ddcal%02i.reg' % i)

    t['name'] = ['ddcal%02i' % i for i in xrange(len(t))]

    return t


def voronoi_finite_polygons_2d_box(vor, box):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    box.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    box : (2,2) float array
        corners of bounding box
        numpy.array([[x1,y1],[x2,y2]])

    Returns
    -------
    poly : array of M (N,2) arrays
        polygon coordinates for M revised Voronoi regions.

    """
    import matplotlib.transforms as mplTrans
    import matplotlib.path as mplPath

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")
    if box.shape != (2,2):
        raise ValueError("Bounding box should be 2x2 array ((x1,y1),(x2,y2))")

    radius = np.max(box)

    # define the bounding box transform from the box extent - to be used to intersect with the regions
    bbox = mplTrans.Bbox(box)

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge
            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    regions, imvertices = new_regions, np.asarray(new_vertices)
    #return new_regions, np.asarray(new_vertices)

    ## now force them to be in the bounding box
    poly = np.asarray([imvertices[v] for v in regions])

    newpoly = []

    for p in poly:
        polyPath = mplPath.Path(p)
        newpolyPath = polyPath.clip_to_bbox(bbox)
        pp = newpolyPath.vertices.transpose()
        newpoly.append(pp.transpose())

    return np.asarray(newpoly)


def make_tassellation(t, fitsfile, outdir='regions/', pb_cut=5.):
    """
    t : table with directions
    firsfile : model fits file to tassellate
    pb_cut : primary beam diameter cut (~fwhm) in degree,
             sources outside will not have a facet
    """
    logging.debug("Image used for tasselation reference: "+fitsfile)
    fits = pyfits.open(fitsfile)
    hdr, data = flatten(fits)
    w = wcs.WCS(hdr)

    # Add facet size column
    t['facet_size'] = np.zeros(len(t))
    t['facet_size'].unit = u.degree

    # Cut at FWHM: if source outside do not do facet
    ra_c, dec_c = w.all_pix2world(data.shape[0]/2., data.shape[1]/2., 0, ra_dec_order=True)
    x_c = data.shape[0]/2.
    y_c = data.shape[1]/2.
    idx_for_facet = []
    for i, dd in enumerate(t):
        if angsep(dd['RA'], dd['DEC'], ra_c, dec_c) < pb_cut/2.:
            idx_for_facet.append(i)
    #print idx_for_facet

    # convert to pixel space (voronoi must be in eucledian space)
    t['x'], t['y'] = w.all_world2pix(t['RA'],t['DEC'], 0, ra_dec_order=True)
    pixsize = np.abs(hdr['CDELT1'])
    beam_radius_pix = (pb_cut/2.)/pixsize
    x1 = x_c - beam_radius_pix
    y1 = y_c - beam_radius_pix
    x2 = x_c + beam_radius_pix
    y2 = y_c + beam_radius_pix

    # do tasselization
    vor = Voronoi(np.array((t['x'][idx_for_facet], t['y'][idx_for_facet])).transpose())
    box = np.array([[x1,y1],[x2,y2]])
    impoly = voronoi_finite_polygons_2d_box(vor, box)

    # plot tasselization
    import matplotlib.pyplot as pl
    pl.figure(figsize=(8,8))
    c1 = pl.Circle((x_c, y_c), beam_radius_pix, color='g', fill=False)
    ax1 = pl.gca()
    ax1.add_artist(c1)
    ax1.plot(t['x'],t['y'],'*')
    ax1.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1])
    for p in impoly:
        pp = p.transpose()
        ax1.plot(pp[0],pp[1])
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
    pl.savefig('voronoi_facets.png')

    # save regions
    for i, poly in enumerate(impoly):
        ra, dec = w.all_pix2world(poly[:,0],poly[:,1], 0, ra_dec_order=True)
        coords = np.array([ra,dec]).T.flatten()

        s = Shape('Polygon', None)
        s.coord_format = 'fk5'
        s.coord_list = coords # ra, dec, radius
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '5', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        s.comment = 'color=red text=""'

        regions = pyregion.ShapeList([s])
        regions.write(outdir+t['name'][i]+'-facet.reg')

        t['facet_size'][idx_for_facet[i]] = pixsize * np.max( [np.max(poly[:,0]) - np.min(poly[:,0]), np.max(poly[:,1]) - np.min(poly[:,1])] )

    t.remove_columns(['x','y'])

    return t

