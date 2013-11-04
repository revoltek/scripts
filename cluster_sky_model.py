#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013 - Reinout van Weeren
# Copyright (C) 2013 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import os, sys
import matplotlib
matplotlib.use('GTK')
import numpy as np
import pylab as pl
import matplotlib, itertools
import pyrap.tables as pt
import lofar.stationresponse as lsr
from lib_coordinates_mode import *
from lib_read_skymodel import read_skymodel

def compute_patch_center(data, beam_ms = None):
    """
    Return the patches names, central (weighted) RA and DEC and total flux
    """
    patches = np.unique(data['Patch'])
    print 'These are', len(patches) ,' input patches'
    ra_patches   = []
    dec_patches  = []
    flux_patches = []

    # get the average time of the obs, OK for 1st order correction
    if beam_ms != None:
        t = pt.table(beam_ms, ack=False)
        tt = t.query('ANTENNA1==0 AND ANTENNA2==1', columns='TIME')
        time = tt.getcol("TIME")
        time = min(time) + ( max(time) - min(time) ) / 2.
        t.close()
        sr = lsr.stationresponse(beam_ms, False, True)

    for patch_id, patch in enumerate(patches):
        components_idx = np.where(data['Patch'] == patch)[0]
        print 'Patch', patch, 'has', len(components_idx), 'components'

        # set to zero
        ra_patch    = 0.
        dec_patch   = 0.
        ra_weights  = 0.
        dec_weights = 0.
        flux_patch  = 0.

        for component_idx in components_idx:

            # conver RA, DEC to degrees for component
            ra_comp   = data['Ra'][component_idx]
            dec_comp  = data['Dec'][component_idx]
            flux_comp = np.float(data['I'][component_idx])

            # beam correction
            if beam_ms != None:
                sr.setDirection(ra_comp*np.pi/180.,dec_comp*np.pi/180.)
                # use station 0 to compute the beam and get mid channel
                beam = sr.evaluateStation(time,0)
                r = abs(beam[int(len(beam)/2.)])
                beam = ( r[0][0] + r[1][1] ) / 2.
                print "Beam:", beam,
                flux_comp *= beam

            # calculate the average weighted patch center, and patch flux
            flux_patch  += flux_comp
            ra_patch    = ra_patch + (flux_comp*ra_comp)
            dec_patch   = dec_patch+ (flux_comp*dec_comp)
            ra_weights  += flux_comp
            dec_weights += flux_comp

        print 'Center RA', ra_patch/ra_weights,
        print 'Center DEC', dec_patch/dec_weights,
        print 'Flux', flux_patch

        ra_patches.append(ra_patch/ra_weights)
        dec_patches.append(dec_patch/dec_weights)
        flux_patches.append(flux_patch)

    return np.array(patches), np.array(ra_patches), np.array(dec_patches), np.array(flux_patches)


def create_clusters(patches, ra_patches, dec_patches, flux_patches, Q, show_plot):
    """
    Clusterize all the patches of the skymodel iteratively around the brightest patches
    """

    # sort the patches by brightest first
    idx = np.argsort(flux_patches)[::-1] # -1 to reverse sort
    flux_patches = flux_patches[idx]
    ra_patches   = ra_patches[idx]
    dec_patches  = dec_patches[idx]
    patches      = patches[idx]

    patch_cluster_id = np.zeros(len(patches[Q:]), dtype=int)

    # find the Q brightest patches
    clusters = np.copy(patches[0:Q])
    clusters_ra = np.copy(ra_patches[0:Q])
    clusters_dec = np.copy(dec_patches[0:Q])
    clusters_flux= np.copy(flux_patches[0:Q])

    # find the closest cluster to a given patch to initialize the clustering
    patch_cluster_id = find_closest_cluster(clusters,clusters_ra,clusters_dec,\
                       patches,ra_patches,dec_patches,Q)

    clusters_ra_new, clusters_dec_new = update_cluster_centers(clusters,clusters_ra,\
                                        clusters_dec,clusters_flux,ra_patches,      \
                                        dec_patches,flux_patches,patch_cluster_id,Q)

    patch_cluster_id_new = find_closest_cluster(clusters,clusters_ra_new,clusters_dec_new,\
                                                patches,ra_patches,dec_patches,Q)

    # Iterate until clusters are properly defined
    count = 1
    while (patch_cluster_id != patch_cluster_id_new).any() == True:

        clusters_ra_new, clusters_dec_new = np.copy(update_cluster_centers(clusters,clusters_ra_new,\
                                                       clusters_dec_new,clusters_flux,ra_patches,      \
                                                       dec_patches,flux_patches,patch_cluster_id,Q))
        patch_cluster_id = np.copy(patch_cluster_id_new)
        patch_cluster_id_new = find_closest_cluster(clusters,clusters_ra_new,clusters_dec_new,\
                               patches,ra_patches,dec_patches,Q)
        print 'Iteration', count
        count += 1

    print "Central patch fluxes:", clusters_flux
    clusters_flux = compute_total_cluster_flux(clusters_flux,flux_patches,patch_cluster_id_new,Q)
    print "Cluster fluxes:", clusters_flux

    if show_plot:
        # check out http://phrogz.net/css/distinct-colors.html
        color=itertools.cycle(["#d96c6c", "#b23000", "#402310", "#ffd9bf", "#998273", "#b25f00", "#ffaa00", "#736039",
                               "#332f26", "#fff240", "#838c00", "#d5d9a3", "#aaff00", "#334d00", "#657356", "#bfffc8",
                               "#00e65c", "#008c38", "#00ffee", "#004d47", "#59adb3", "#00c2f2", "#267399", "#738c99",
                               "#0066ff", "#79aaf2", "#264599", "#0000ff", "#0000cc", "#000033", "#323040", "#220080",
                               "#b866cc", "#9d7ca6", "#ee00ff", "#59164c", "#ff40a6", "#330d21", "#8c4662", "#ff0022",
                               "#73000f", "#ffbfc8"])
        for cluster_id, cluster in enumerate(clusters):
            c = color.next()

            # plot patches
            idx = np.where(patch_cluster_id == cluster_id)
            for el in range(len(idx[0])):
                matplotlib.pyplot.plot(ra_patches[idx[0][el]+Q],dec_patches[idx[0][el]+Q],'o',color=c,markersize=10)

            # plot clusters
            matplotlib.pyplot.plot(clusters_ra_new[cluster_id],clusters_dec_new[cluster_id],'*',color=c,markersize=25)
            matplotlib.pyplot.xlabel('RA [deg]')
            matplotlib.pyplot.ylabel('DEC [deg]')

            matplotlib.pyplot.text(clusters_ra_new[cluster_id]+0.05, clusters_dec_new[cluster_id]+0.05, str(cluster_id))

        # Reverse RA axis because of RA definition
        ax = matplotlib.pyplot.gca()
        ax.invert_xaxis()

        matplotlib.pyplot.show()

    # clusters_ra_new are the centroids, note return new values (they do not match anymore with the defining patch!!)
    return clusters,clusters_ra_new,clusters_dec_new,clusters_flux,patches,ra_patches,dec_patches,patch_cluster_id_new


def write_skymodel(clusters,clusters_ra,clusters_dec,clusters_flux,patches,patch_cluster_id_new,\
                   data,Q,outfilename,cluster_ra_min_max, cluster_dec_min_max) :
    """
    Write the outfile, note that clusters ra and dec are recaluclated as the mean between min and max and not the centroid
    """

    outfile =open(outfilename,'w')
    outfile.write("#(Name, Type, Patch, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, Orientation, ReferenceFrequency='1.49803e+08', SpectralIndex='[]') = format\n")
    outfile.write("\n")
    outfile.write("\n")
    outfile.write("## updated clustered skymodel\n")
    outfile.write("\n")
    outfile.write("\n")

    def ratohms_string(ra):
        rah, ram, ras = ratohms(ra)
        return str(rah) + ':' + str(ram) + ':' +str(round(ras,2))
    def dectodms_string(dec):
        decd, decm, decs = dectodms(dec)
        return str(decd) + '.' + str(decm) + '.' +str(round(decs,2))

    for cluster_id, cluster in enumerate(clusters):

        clustername = 'CLUSTER_' + str(cluster_id)

        # compute RA, DEC patch in sexagesimal format (from degr)
        mid_ra  = np.mean(cluster_ra_min_max[cluster_id,:])
        mid_dec = np.mean(cluster_dec_min_max[cluster_id,:])

        outfile.write(", , " + clustername + ", " + ratohms_string(mid_ra) + ", " + dectodms_string(mid_dec) + "\n")
        outfile.write("# NDPPP patch coordinates: " + ratohms_string(mid_ra) + ", " + dectodms_string(mid_dec) + "\n")
        outfile.write("# Total cluster flux="+str(round(clusters_flux[cluster_id],4))+ " Jy\n")

        dec_size = (np.max(cluster_dec_min_max[cluster_id,:])-np.min(cluster_dec_min_max[cluster_id,:]))*3600.
        ra_size =  (np.max(cluster_ra_min_max[cluster_id,:])-np.min(cluster_ra_min_max[cluster_id,:]))*3600.*np.cos(np.pi*mid_dec/180.)

        mmstr = 'RAsize, DECsize [arcsec] = '+ str([np.int(ra_size),np.int(dec_size)])

        outfile.write("# "+ mmstr + "\n")

        # write paatches defining cluster first
        idx = np.where(data['Patch'] == cluster)[0]
        for el in idx:

            linep1 = data['Name'][el]+', '+data['Type'][el]+', '+clustername+', '+ratohms_string(data['Ra'][el])+', '+dectodms_string(data['Dec'][el])+', '
            linep2 = str(data['I'][el])+', '+str(data['Q'][el]) +', '+str(data['U'][el]) +', '+str(data['V'][el])+', '
            linep3 = str(data['MajorAxis'][el])+', '+str(data['MinorAxis'][el]) +', '+str(data['Orientation'][el]) +', '
            linep4 = str(data['ReferenceFrequency'][el])+', '+ data['SpectralIndex'][el]
            outfile.write(linep1+linep2+linep3+linep4+'\n')

        #write other patches in this cluster
        idx = np.where(patch_cluster_id_new == cluster_id)[0]
        patches_list = patches[idx+Q]
        for patch in patches_list:

            idx = np.where(data['Patch'] == patch)[0]
            for el in (idx):
                linep1 = data['Name'][el]+', '+data['Type'][el]+', '+clustername+', '+ratohms_string(data['Ra'][el])+', '+dectodms_string(data['Dec'][el])+', '
                linep2 = str(data['I'][el])+', '+str(data['Q'][el]) +', '+str(data['U'][el]) +', '+str(data['V'][el])+', '
                linep3 = str(data['MajorAxis'][el])+', '+str(data['MinorAxis'][el]) +', '+str(data['Orientation'][el]) +', '
                linep4 = str(data['ReferenceFrequency'][el])+', '+ data['SpectralIndex'][el]
                outfile.write(linep1+linep2+linep3+linep4+'\n')

        outfile.write("\n")
        outfile.write("\n")

    outfile.close()


def compute_total_cluster_flux(clusters_flux,flux_patches,patch_cluster_id_new,Q):
    """
    Return a list of the total flux in all clusters
    """
    clusters_flux_new = np.copy(clusters_flux) # copy into new array

    for cluster_id in range(len(clusters_flux)):
        idx = np.where(patch_cluster_id_new == cluster_id)[0]

        # add the flux from the mathcing patches
        for el in xrange(len(idx)):
            clusters_flux_new[cluster_id] = clusters_flux_new[cluster_id] + flux_patches[idx[el]+Q]

    return clusters_flux_new


def find_closest_cluster(clusters,clusters_ra,clusters_dec,patches,ra_patches,dec_patches,Q):
    """
    Return an array of dimension "len(patches) - Q" which indicate for each non-center patch the closest cluster
    """
    patch_cluster_id = np.zeros(len(patches[Q:]), dtype=int)
    for (patch_id, patch) in enumerate(patches[Q:]): #select patches larger than Q

        #print 'Finding closest cluster for patch', patch
        angulardistance = 1.0e10

        for (cluster_id,cluster) in enumerate(clusters): # loop over clusters
            ra_patch  = ra_patches[patch_id+Q]
            dec_patch = dec_patches[patch_id+Q]
            adis = abs(angsep2(clusters_ra[cluster_id],clusters_dec[cluster_id],ra_patch,dec_patch))
            if adis < angulardistance:
                angulardistance = adis
                patch_cluster_id[patch_id] = cluster_id
                #print abs(angsep2(clusters_ra[cluster_id],clusters_dec[cluster_id],ra_patch,dec_patch))
                #print 'cloest cluster found till now', cluster
        #print 'Closest cluster to', patch, 'is at', angulardistance ,\
        #      '[cluster is', clusters[patch_cluster_id[patch_id]],']'

    return patch_cluster_id


def update_cluster_centers(clusters,clusters_ra,clusters_dec,clusters_flux,ra_patches,dec_patches,flux_patches,patch_cluster_id,Q):
    """
    Return the new cluster centres give a patch_cluster_id
    """
    clusters_ra_new  = np.zeros(len(clusters_ra))
    clusters_dec_new =  np.zeros(len(clusters_dec))
    for (cluster_id,cluster) in enumerate(clusters):
        cl_ra      = 0.
        cl_dec     = 0.
        weight_dec = 0.
        weight_ra  = 0.
        cl_ra      = clusters_ra[cluster_id]*clusters_flux[cluster_id]
        cl_dec     = clusters_dec[cluster_id]*clusters_flux[cluster_id]
        weight_ra  = clusters_flux[cluster_id]
        weight_dec = clusters_flux[cluster_id]

        #print 'input cluster coordinates', cl_ra/weight_ra, cl_dec/weight_dec
        idx = np.where(patch_cluster_id == cluster_id)[0]
        for el in xrange(len(idx)):
            cl_ra      = cl_ra + ra_patches[idx[el]+Q]*flux_patches[idx[el]+Q]
            cl_dec     = cl_dec + dec_patches[idx[el]+Q]*flux_patches[idx[el]+Q]
            weight_ra  = weight_ra + flux_patches[idx[el]+Q]
            weight_dec = weight_dec+ flux_patches[idx[el]+Q]

        cl_ra = cl_ra/weight_ra
        cl_dec= cl_dec/weight_dec
        clusters_ra_new[cluster_id] = cl_ra
        clusters_dec_new[cluster_id] = cl_dec

    return clusters_ra_new, clusters_dec_new


def compute_cluster_info(clusters,clusters_ra,clusters_dec,patches,ra_patches,dec_patches,patch_cluster_id_new,Q):
    """
    Return min and max clusters ra and dec
    """
    cluster_ra_min_max  = np.zeros((len(clusters),2))
    cluster_dec_min_max = np.zeros((len(clusters),2))

    for cluster_id, cluster in enumerate(clusters):

        idx = np.copy(np.where(patch_cluster_id_new == cluster_id)[0])
        idx = idx + Q

        idx = np.sort(np.append(idx, cluster_id))

        ra_min = np.min(ra_patches[idx])
        ra_max = np.max(ra_patches[idx])
        dec_min= np.min(dec_patches[idx])
        dec_max= np.max(dec_patches[idx])

        cluster_ra_min_max[cluster_id,:] = [ra_min, ra_max]
        cluster_dec_min_max[cluster_id,:] = [dec_min, dec_max]

    return cluster_ra_min_max, cluster_dec_min_max


if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <PyBDSM BBS catalog>')
    parser.add_option('-o', '--outfile', dest='file_clusters', help='Name of the output skymodel file [default catalogue.clusters.skymodel]', \
                      metavar='VAL', default='catalog.clusters.skymodel')
    parser.add_option('-n', '--numclusters', dest='Q', help='Number of clusters [default: 10]', metavar='VAL', default=10, type=int)
    parser.add_option('-p', action='store_true', dest='show_plot', help='Show plot of the cluster disposition', default=False)
    parser.add_option('-b', '--beamMS', help='Give an MS used to find freq/dir/ant and correct for beam effect when clustering [default=None]', default=None)
    (o, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit()

    data = read_skymodel(args[0])
    patches, ra_patches, dec_patches, flux_patches = compute_patch_center(data, o.beamMS)

    clusters,clusters_ra,clusters_dec,clusters_flux,patches,ra_patches,dec_patches,patch_cluster_id = \
               create_clusters(patches, ra_patches, dec_patches, flux_patches, o.Q, o.show_plot)

    cluster_ra_min_max, cluster_dec_min_max = compute_cluster_info(clusters,clusters_ra,clusters_dec,patches,ra_patches, dec_patches,patch_cluster_id, o.Q)

    write_skymodel(clusters,clusters_ra,clusters_dec,clusters_flux,patches,patch_cluster_id,data,o.Q,o.file_clusters, cluster_ra_min_max, cluster_dec_min_max)

    print "Done."
