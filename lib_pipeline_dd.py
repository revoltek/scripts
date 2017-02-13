#!/usr/bin/python

import os, sys
import logging
import numpy as np
import lsmtool
import shapely.geometry
from shapely.ops import cascaded_union
from itertools import combinations
from astropy.coordinates import Angle

def make_directions_file_from_skymodel(filename, flux_min_Jy=4.0, size_max_arcmin=2.0,
    directions_separation_max_arcmin=10.0, directions_max_num=None, flux_min_for_merging_Jy=1.0):
    """
    Selects appropriate calibrators from sky model and makes the directions file
    Parameters
    ----------
    s : LSMTool SkyModel object
        Skymodel made by grouping clean components of dir-independent model
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
    directions_file : str
        Filename of resulting Factor-formatted directions file
    """
    s = lsmtool.load(filename)

    # Filter patches on size
    s_large = s.copy()
    sizes = s.getPatchSizes(units='arcmin', weight=True)
    s.select(sizes < size_max_arcmin, aggregate=True, force=True)
    if len(s) == 0:
        logging.critical("No sources found that meet the specified max size criterion.")
        sys.exit(1)
    logging.info('Found {0} sources with sizes below {1} '
        'arcmin'.format(len(s.getPatchNames()), size_max_arcmin))
    s_large.remove(sizes < size_max_arcmin, aggregate=True, force=True)

    # Filter patches on flux density limit for merging
    s.select('I > {0} Jy'.format(flux_min_for_merging_Jy), aggregate='sum', force=True)
    if len(s) == 0:
        logging.critical("No sources found above {} Jy.".format(flux_min_for_merging_Jy))
        sys.exit(1)
    logging.info('Found {0} sources with flux densities above {1} Jy'.format(
        len(s.getPatchNames()), flux_min_for_merging_Jy))
    if len(s_large) > 0:
        s_large.select('I > {0} Jy'.format(flux_min_for_merging_Jy), aggregate='sum', force=True)

    # Look for nearby pairs
    pRA, pDec = s.getPatchPositions(asArray=True)
    for ra, dec in zip(pRA.tolist()[:], pDec.tolist()[:]):
        dist = s.getDistance(ra, dec, byPatch=True, units='arcmin')
        nearby = np.where(dist < directions_separation_max_arcmin)
        if len(nearby[0]) > 1:
            patches = s.getPatchNames()[nearby]
            s.merge(patches.tolist())
    # update patch positions
    s.setPatchPositions(method='mid')
    logging.info('Souece merged in {0} groups within {1} arcmin of each other.'.format(
        len(s.getPatchNames()), directions_separation_max_arcmin))

    # Filter patches on total flux density limit
    s.select('I > {0} Jy'.format(flux_min_Jy), aggregate='sum', force=True)
    if len(s) == 0:
        logging.critical("No sources or merged groups found that meet the specified "
            "min total flux density criterion.")
        sys.exit(1)
    logging.info('Found {0} sources or merged groups with total flux densities above {1} Jy'.format(
        len(s.getPatchNames()), flux_min_Jy))

    # Trim directions list to get directions_max_num of directions
    if directions_max_num is not None:
        dir_fluxes = s.getColValues('I', aggregate='sum')
        dir_fluxes_sorted = dir_fluxes.tolist()
        dir_fluxes_sorted.sort(reverse=True)
        cut_jy = dir_fluxes_sorted[-1]
        while len(dir_fluxes_sorted) > directions_max_num:
            cut_jy = dir_fluxes_sorted.pop() + 0.00001
        s.remove('I < {0} Jy'.format(cut_jy), aggregate='sum')
    logging.info('Kept {0} directions in total'.format(len(s.getPatchNames())))

    # Now, merge calibrator patches with any extended sources that meet the
    # flux density limit for merging
    if len(s_large) > 0:
        logging.info('Merging extended sources within {0} arcmin of calibrators...'.format(
            directions_separation_max_arcmin))
        calibrator_names = s.getPatchNames().tolist()
        # calibrator positions
        pRA, pDec = s.getPatchPositions(asArray=True)
        s.concatenate(s_large)
        for ra, dec in zip(pRA.tolist()[:], pDec.tolist()[:]):
            dist = s.getDistance(ra, dec, byPatch=True, units='arcmin')
            nearby = np.where(dist < directions_separation_max_arcmin)
            if len(nearby[0]) > 1:
                patches = s.getPatchNames()[nearby].tolist()

                # Ensure that calibrator patch is first in list (as merged
                # patch will get its name). If there are two calibrator patches
                # in the list, use the first one
                for calibrator_name in calibrator_names:
                    if calibrator_name in patches:
                        patches.remove(calibrator_name)
                        patches.insert(0, calibrator_name)
                        break
                s.merge(patches)

        # Remove any non-calibrator patches from the merged model
        all_names = s.getPatchNames().tolist()
        calibrator_ind = np.array([all_names.index(calibrator_name) for
            calibrator_name in calibrator_names if calibrator_name in all_names])
        s.select(calibrator_ind, aggregate=True)

    # Logs
    s.info()

    # Write region files, one per patch
    s.setPatchPositions(method='mid')
    for i, patch in enumerate(s.getPatchNames()):
        w = s.getRowValues(patch)
        w.write('regions/src%02i.reg' % i, format='ds9')


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
