#!/usr/bin/python
# Input: 
# CASA region file with one region for each DD cal (every shape is OK)
# Output:
# region files for each facet

ddreg = 'dd.crtf'

import matplotlib.pyplot as pl
import numpy as np
import os, sys, glob, re
import pyrap.images
from lib_pipeline import *
from lib_coordinates_mode import *
set_logger()
try:
    from scipy.spatial import Voronoi
except:
    logging.error("Load latest scipy with 'source /opt/cep/scripts/doPythonlibs'")
    sys.exit(1)

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

    radius = numpy.max(box)

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
            t /= numpy.linalg.norm(t)
            n = numpy.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = numpy.sign(numpy.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = numpy.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = numpy.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = numpy.array(new_region)[numpy.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    regions, imvertices = new_regions, numpy.asarray(new_vertices)
    #return new_regions, numpy.asarray(new_vertices)

    ## now force them to be in the bounding box
    poly = numpy.asarray([imvertices[v] for v in regions])

    newpoly = []

    for p in poly:
        polyPath = mplPath.Path(p)
        newpolyPath = polyPath.clip_to_bbox(bbox)
        pp = newpolyPath.vertices.transpose()
        newpoly.append(pp.transpose())

    return numpy.asarray(newpoly)


# used to move to pixel space (voronoi must be in eucledian space) and back
model = glob.glob('self/models/wide_g*.model*')[0]
logging.debug("Model image is "+model)
img = pyrap.images.image(model)
c = img.coordinates()
# assumes images in a standard casa shape
assert c.get_axes()[2] == ['Declination', 'Right Ascension']
# assume same increment in image axes
assert abs(c.get_increment()[2][0]) == abs(c.get_increment()[2][1])

# find a "center" for each dd cal
logging.info("Find the center of all dd cal regions")
dds = []
i = 0
with open(ddreg) as f:
    for dd in f.readlines():
        if dd[0] == '\n': continue
        if dd[0] == '#': 
            init = dd
            continue

        # return an [(ra,dec), (ra,dec), ..] containing all matches
        c = re.findall("(?<=\[)(\d+:\d+:\d+\.\d+),\s?([+-]\d+\.\d+\.\d+\.\d+)(?=\])", dd)

        # get from baricenter
        if 'poly' == dd[:4]:
            ra = []
            dec = []
            for cra, cdec in c:
                ra.append( hmstora(*cra.split(':')) )
                dec.append( dmstodec(*cdec.split('.', 2)) )
            dds.append([ np.mean([min(ra),max(ra)]) , np.mean([min(dec),max(dec)]) ])

        # get from center coord
        elif 'ellipse' == dd[:7]:
            ra, dec = c[0]
            dds.append([hmstora(*ra.split(':')), dmstodec(*dec.split('.', 2))])

        # get from mean of angles
        elif 'box' == dd[:3]:
            ra, dec = c[0]
            c0 = [hmstora(*ra.split(':')), dmstodec(*dec.split('.', 2))]
            ra, dec = c[1]
            c1 = [hmstora(*ra.split(':')), dmstodec(*dec.split('.', 2))]
            dds.append([ np.mean([c0[0],c1[0]]) , np.mean([c0[1],c1[1]]) ])

        else:
            logging.error('Unsopported shape: '+dd)
            sys.exit(1)

        # save dd region files
        logging.info('writing reg_dd-%03d.crtf' % i)
        with open('reg_dd-%03d.crtf' % i, 'w') as fo:
            fo.write(init)
            fo.write(dd)

        logging.debug(' - centre: RA:'+str(dds[-1][0])+' - DEC:'+str(dds[-1][1]))
        i = i+1

# convert to pixel space (NOTE: [y,x])
dds_p = []
for dd in dds:
    dds_p.append( img.topixel([1, 1, dd[1]*np.pi/180., dd[0]*np.pi/180.])[2:] )
dds_p = np.array(dds_p)

# do tasselization
vor = Voronoi(np.array((dds_p[:,1], dds_p[:,0])).transpose())
x1, y1 = 0, 0
x2, y2 = img.shape()[2:]
box = np.array([[x1,y1],[x2,y2]])
impoly = voronoi_finite_polygons_2d_box(vor, box)

# convert back to ra/dec
impoly_p = np.copy(impoly)
for poly in impoly:
    for angle in poly:
        f, c, angle[0], angle[1] = img.toworld([0,0,angle[0],angle[1]])

# save facet region files
sall = '#CRTFv0 CASA Region Text Format version 0\n'
for i, dd_facet in enumerate(impoly):
    s= '#CRTFv0 CASA Region Text Format version 0\n'
    s1 = 'poly ['
    dirs = ['[{0:.0f}:{1:.0f}:{2:.4f}, {3:+.0f}.{4:.0f}.{5:.4f}]'.format( *( ratohms(p[1]*180/np.pi) + dectodms(p[0]*180/np.pi) ) ) for p in impoly[i]]
    s1 += ', '.join(dirs) 
    s1 += '] coord=J2000, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=magenta, font="DejaVu Sans", fontsize=11, fontstyle=normal, usetex=false\n'
    with open('reg_facet-%03d.crtf' % i, 'w') as fo:
        fo.write(s+s1)
    sall += s1
with open('reg_facet-all.crtf', 'w') as fo:
    fo.write(sall)

# plot tasselization
pl.figure(figsize=(8,8))
ax1 = pl.subplot(111)
ax1.plot(dds_p.transpose()[1],dds_p.transpose()[0],'*')
ax1.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1])
for p in impoly_p:
    pp = p.transpose()
    ax1.plot(pp[0],pp[1])
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
pl.savefig('voronoi_facets.png')
