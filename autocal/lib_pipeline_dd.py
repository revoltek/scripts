#!/usr/bin/python

import os, sys
import logging
import numpy as np
from astropy.coordinates import Angle, SkyCoord, match_coordinates_sky
import astropy.units as u
from astropy.table import Table

def make_directions_from_skymodel(filename, outdir='regions/' flux_min_Jy=1.0, size_max_arcmin=3.0,
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
    logging.info('# sources initial: %i' % len(t))

    # exclude sources that are too extended
    t_large = t[ (t['Maj'] >= size_max_arcmin*u.arcmin) ]
    t = t[ (t['Maj'] < size_max_arcmin*u.arcmin) ]
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
                        s['Maj'] = max(s['Maj'], t['Maj'][i]) + dist.degree
                        updated = True
                    # other wise weighted mean
                    elif t['Peak_flux'][i] > s['Peak_flux']:
                        s['RA'] = (s['RA']*s['Peak_flux'] + t['RA'][i]*t['Peak_flux'][i])/(s['Peak_flux']+t['Peak_flux'][i])
                        s['DEC'] = (s['DEC']*s['Peak_flux'] + t['DEC'][i]*t['Peak_flux'][i])/(s['Peak_flux']+t['Peak_flux'][i])
                        s['Maj'] = max(s['Maj'], t['Maj'][i]) + dist.degree/2.
                        updated = True
                    else:
                        s['Maj'] = max(s['Maj'], t['Maj'][i]) + dist.degree

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
                t['Maj'][i] = max(s['Maj'], t['Maj'][i]) + dist.degree

    t.sort('Total_flux')
    t.reverse()

    print t

    # Writedd global region files
    table_to_circ_region(t, outdir+'/all.reg')

    # save source by source
    for i, s in enumerate(t):
        table_to_circ_region(t[i:i+1], outdir+'/ddcal%02i.reg' % i)

    t['name'] = ['ddcal%02i' % i for i in xrange(len(t))]
    t[''] = ['ddcal%02i' % i for i in xrange(len(t))]
    return t


def table_to_circ_region(table, outfile, racol='RA', deccol='DEC', sizecol='Maj', color='green'):

    import pyregion
    from pyregion.parser_helper import Shape

    color = [color] * len(table)

    regions = []
    for i, r in enumerate(table):
        s = Shape('circle', None)
        s.coord_format = 'fk5'
        s.coord_list = [ r[racol], r[deccol], r[sizecol] ]
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '5', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        s.comment = 'color={} text="{}"'.format(color[i], str(i))
        regions.append(s)

    regions = pyregion.ShapeList(regions)
    regions.write(outfile)


def thiessen(directions_list, field_ra_deg, field_dec_deg, faceting_radius_deg,
    s=None, check_edges=False, target_ra=None, target_dec=None,
    target_radius_arcmin=None, beam_ratio=None):
    """
    Generates and add thiessen polygons or patches to input directions
    Parameters
    ----------
    directions_list : list of Direction objects
        List of input directions
    field_ra_deg : float
        RA in degrees of field center
    field_dec_deg : float
        Dec in degrees of field center
    faceting_radius_deg : float
        Maximum radius within which faceting will be done. Direction objects
        with postions outside this radius will get small rectangular patches
        instead of thiessen polygons
    s : LSMTool SkyModel object, optional
        Sky model to use to check for source near facet edges
    check_edges : bool, optional
        If True, check whether any know source falls on a facet edge. If sources
        are found that do, the facet is adjusted
    target_ra : str, optional
        RA of target source. E.g., '14h41m01.884'
    target_dec : str, optional
        Dec of target source. E.g., '+35d30m31.52'
    target_radius_arcmin : float, optional
        Radius in arcmin of target source
    beam_ratio : float, optional
        Ratio of semi-major (N-S) axis to semi-minor (E-W) axis for the primary
        beam
    """

    # Select directions inside FOV (here defined as ellipse given by
    # faceting_radius_deg and the mean elevation)
    faceting_radius_pix = faceting_radius_deg / 0.066667 # radius in pixels
    field_x, field_y = radec2xy([field_ra_deg], [field_dec_deg],
        refRA=field_ra_deg, refDec=field_dec_deg)

    fx = []
    fy = []
    for th in range(0, 360, 1):
        fx.append(faceting_radius_pix * np.cos(th * np.pi / 180.0) + field_x[0])
        fy.append(faceting_radius_pix * beam_ratio * np.sin(th * np.pi / 180.0) + field_y[0])
    fov_poly_tuple = tuple([(xp, yp) for xp, yp in zip(fx, fy)])
    fov_poly = Polygon(fx, fy)

    points, _, _ = getxy(directions_list, field_ra_deg, field_dec_deg)
    for x, y, d in zip(points[0], points[1], directions_list):
        dist = fov_poly.is_inside(x, y)
        if dist < 0.0:
            # Source is outside of FOV, so use simple rectangular patches
            d.is_patch = True

    # Now do the faceting (excluding the patches)
    directions_list_thiessen = [d for d in directions_list if not d.is_patch]
    points, _, _ = getxy(directions_list_thiessen, field_ra_deg, field_dec_deg)
    points = points.T

    # Generate array of outer points used to constrain the facets
    nouter = 64
    means = np.ones((nouter, 2)) * points.mean(axis=0)
    offsets = []
    angles = [np.pi/(nouter/2.0)*i for i in range(0, nouter)]
    for ang in angles:
        offsets.append([np.cos(ang), np.sin(ang)])

    # Generate initial facets
    radius = 5.0 * faceting_radius_deg / 0.066667 # radius in pixels
    scale_offsets = radius * np.array(offsets)
    outer_box = means + scale_offsets
    points_all = np.vstack([points, outer_box])
    tri = Delaunay(points_all)
    circumcenters = np.array([_circumcenter(tri.points[t])
                              for t in tri.vertices])
    thiessen_polys = [_thiessen_poly(tri, circumcenters, n)
                      for n in range(len(points_all) - nouter)]

    # Check for vertices that are very close to each other, as this gives problems
    # to the edge adjustment below
    for thiessen_poly in thiessen_polys:
        dup_ind = 0
        for i, (v1, v2) in enumerate(zip(thiessen_poly[:-1], thiessen_poly[1:])):
            if (approx_equal(v1[0], v2[0], rel=1e-6) and
                approx_equal(v1[1], v2[1], rel=1e-6)):
                thiessen_poly.pop(dup_ind)
                dup_ind -= 1
            dup_ind += 1

    # Clip the facets at FOV
    for i, thiessen_poly in enumerate(thiessen_polys):
        polyv = np.vstack(thiessen_poly)
        poly_tuple = tuple([(xp, yp) for xp, yp in zip(polyv[:, 0], polyv[:, 1])])
        p1 = shapely.geometry.Polygon(poly_tuple)
        p2 = shapely.geometry.Polygon(fov_poly_tuple)
        if p1.intersects(p2):
            p1 = p1.intersection(p2)
            xyverts = [np.array([xp, yp]) for xp, yp in
                zip(p1.exterior.coords.xy[0].tolist(),
                p1.exterior.coords.xy[1].tolist())]
            thiessen_polys[i] = xyverts

    # Check for sources near / on facet edges and adjust regions accordingly
    if check_edges:
        log.info('Adjusting facets to avoid sources...')
        RA, Dec = s.getPatchPositions(asArray=True)
        sx, sy = radec2xy(RA, Dec, refRA=field_ra_deg, refDec=field_dec_deg)
        sizes = s.getPatchSizes(units='degree').tolist()
        fluxes_jy = s.getColValues('I', units='Jy', aggregate='sum').tolist()

        if target_ra is not None and target_dec is not None and target_radius_arcmin is not None:
            log.info('Including target ({0}, {1}) in facet adjustment'.format(
                target_ra, target_dec))
            tra = Angle(target_ra).to('deg').value
            tdec = Angle(target_dec).to('deg').value
            tx, ty = radec2xy([tra], [tdec], refRA=field_ra_deg, refDec=field_dec_deg)
            sx.extend(tx)
            sy.extend(ty)
            sizes.append(target_radius_arcmin*2.0/1.2/60.0)
            fluxes_jy.append(0.0)

        # Set minimum size to 2 - 10 * FWHM of resolution of high-res image, scaled by
        # sqrt(flux) to include strong artifacts in the avoidance region
        fwhm = 25.0 / 3600.0 # degrees
        min_sizes = [fwhm*min(10.0, max(2.0, np.sqrt(flux_jy/0.01))) for flux_jy in fluxes_jy]
        sizes = [max(size, min_size) for size, min_size in zip(sizes, min_sizes)]

        # Filter sources to get only those close to a boundary. We need to iterate
        # until no sources are found
        niter = 0
        while niter < 3:
            niter += 1
            ind_near_edge = []
            for i, thiessen_poly in enumerate(thiessen_polys):
                polyv = np.vstack(thiessen_poly)
                poly_tuple = tuple([(x, y) for x, y in zip(polyv[:, 0], polyv[:, 1])])
                poly = Polygon(polyv[:, 0], polyv[:, 1])
                dists = poly.is_inside(sx, sy)
                for j, dist in enumerate(dists):
                    pix_radius = sizes[j] * 1.2 / 2.0 / 0.066667 # radius of source in pixels
                    if abs(dist) < pix_radius and j not in ind_near_edge:
                        ind_near_edge.append(j)
            if len(ind_near_edge) == 0:
                break
            sx_filt = np.array(sx)[ind_near_edge]
            sy_filt = np.array(sy)[ind_near_edge]
            sizes_filt = np.array(sizes)[ind_near_edge]

            # Adjust all facets for each source near a boundary
            for x, y, size in zip(sx_filt, sy_filt, sizes_filt):
                for i, thiessen_poly in enumerate(thiessen_polys):
                    polyv = np.vstack(thiessen_poly)
                    poly_tuple = tuple([(xp, yp) for xp, yp in zip(polyv[:, 0], polyv[:, 1])])
                    poly = Polygon(polyv[:, 0], polyv[:, 1])
                    dist = poly.is_inside(x, y)
                    p1 = shapely.geometry.Polygon(poly_tuple)

                    pix_radius = size * 1.2 / 2.0 / 0.066667 # size of source in pixels
                    if abs(dist) < pix_radius:
                        p2 = shapely.geometry.Point((x, y))
                        p2buf = p2.buffer(pix_radius)
                        if dist < 0.0:
                            # If point is outside, difference the polys
                            p1 = p1.difference(p2buf)
                        else:
                            # If point is inside, union the polys
                            p1 = p1.union(p2buf)
                        try:
                            xyverts = [np.array([xp, yp]) for xp, yp in
                                zip(p1.exterior.coords.xy[0].tolist(),
                                p1.exterior.coords.xy[1].tolist())]
                            thiessen_polys[i] = xyverts
                        except AttributeError:
                            continue

    # Add the final facet and patch info to the directions
    patch_polys = []
    for d in directions_list:
        # Make calibrator patch
        sx, sy = radec2xy([d.ra], [d.dec], refRA=field_ra_deg, refDec=field_dec_deg)

        # Compute size of patch in pixels, with a factor of 0.6 so that
        # sources are not added along the edges (the areas outside of this 60%
        # region are masked during imaging in the make_clean_mask() call with
        # trim_by = 0.4)
        patch_width = d.cal_imsize * 0.6 * d.cellsize_selfcal_deg / 0.066667
        x0 = sx[0] - patch_width / 2.0
        y0 = sy[0] - patch_width / 2.0
        selfcal_poly = [np.array([x0, y0]),
                np.array([x0, y0+patch_width]),
                np.array([x0+patch_width, y0+patch_width]),
                np.array([x0+patch_width, y0])]

        if d.is_patch:
            # For sources beyond max radius, set facet poly to calibrator poly
            # and clip to facet polys and previous patch polys
            #
            # First, make a copy of the calibrator poly so that it is not
            # altered by the clipping
            patch_poly = [np.copy(vert) for vert in selfcal_poly]

            # Now loop over the facets and patches and clip
            for facet_poly in thiessen_polys + patch_polys:
                polyv = np.vstack(facet_poly)
                poly_tuple = tuple([(xp, yp) for xp, yp in zip(polyv[:, 0], polyv[:, 1])])
                p1 = shapely.geometry.Polygon(poly_tuple)
                p2 = shapely.geometry.Polygon(patch_poly)
                if p2.intersects(p1):
                    p2 = p2.difference(p1)
                    try:
                        xyverts = [np.array([xp, yp]) for xp, yp in
                            zip(p2.exterior.coords.xy[0].tolist(),
                            p2.exterior.coords.xy[1].tolist())]
                        patch_poly = xyverts
                    except AttributeError:
                        pass

            add_facet_info(d, selfcal_poly, patch_poly, field_ra_deg, field_dec_deg)
            patch_polys.append(patch_poly)
        else:
            facet_poly = thiessen_polys[directions_list_thiessen.index(d)]
            add_facet_info(d, selfcal_poly, facet_poly, field_ra_deg, field_dec_deg)
