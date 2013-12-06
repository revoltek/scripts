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

class Patch():
    def __init__(self, name, ra, dec, flux):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.flux = flux

class Cluster():
    def __init__(self, name, patch):
        self.name = name
        self.init_patch = patch
        self.patches = []
        self.update_cluster_coords()

    def add_patch(self, patch):
        """
        Add a patch to this cluster
        """
        self.patches.append(patch)

    def total_cluster_flux(self):
        """
        Return total cluster flux
        """
        cluster_flux = self.init_patch.flux
        for patch in self.patches:
            cluster_flux += patch.flux
    
        return cluster_flux
    
    def update_cluster_coords(self):
        """
        update the self.centroid_ra, self.centroid_dec, self.mid_ra, self.mid_dec
        """
        self.centroid_ra = self.init_patch.ra*self.init_patch.flux
        self.centroid_dec = self.init_patch.dec*self.init_patch.flux
        min_ra = np.inf
        max_ra = -np.inf
        min_dec = np.inf
        max_dec = -np.inf
        
        for patch in self.patches:
            self.centroid_ra += patch.ra*patch.flux
            self.centroid_dec += patch.dec*patch.flux
            if patch.ra < min_ra: min_ra = patch.ra
            if patch.ra > max_ra: max_ra = patch.ra
            if patch.dec < min_dec: min_dec = patch.dec
            if patch.dec > max_dec: max_dec = patch.dec

        self.centroid_ra /= self.total_cluster_flux()
        self.centroid_dec /= self.total_cluster_flux()
        self.mid_ra = min_ra + (max_ra - min_ra)/2.
        self.mid_dec = min_dec + (max_dec - min_dec)/2.

   
def ratohms_string(ra):
    rah, ram, ras = ratohms(ra)
    return str(rah) + ':' + str(ram) + ':' +str(round(ras,2))


def dectodms_string(dec):
    decd, decm, decs = dectodms(dec)
    return str(decd) + '.' + str(decm) + '.' +str(round(decs,2))


def compute_patch_center(data, beam_ms = None):
    """
    Return the patches names, central (weighted) RA and DEC and total flux
    """
    patch_names = np.unique(data['Patch'])
    print 'There are', len(patch_names) ,' input patches'
    patches = []

    # get the average time of the obs, OK for 1st order correction
    if beam_ms != None:
        t = pt.table(beam_ms, ack=False)
        tt = t.query('ANTENNA1==0 AND ANTENNA2==1', columns='TIME')
        time = tt.getcol("TIME")
        time = min(time) + ( max(time) - min(time) ) / 2.
        t.close()
        sr = lsr.stationresponse(beam_ms, False, True)

    for patch_name in patch_names:
        comp_ids = np.where(data['Patch'] == patch_name)[0]
        print 'Patch', patch_name, 'has', len(comp_ids), 'components'

        patch_ra    = 0.
        patch_dec   = 0.
        weights_ra  = 0.
        weights_dec = 0.
        patch_flux  = 0.

        for comp_id in comp_ids:

            # conver RA, DEC to degrees for component
            comp_ra   = data['Ra'][comp_id]
            comp_dec  = data['Dec'][comp_id]
            comp_flux = np.float(data['I'][comp_id])

            # beam correction
            if beam_ms != None:
                sr.setDirection(comp_ra*np.pi/180.,comp_dec*np.pi/180.)
                # use station 0 to compute the beam and get mid channel
                beam = sr.evaluateStation(time,0)
                r = abs(beam[int(len(beam)/2.)])
                beam = ( r[0][0] + r[1][1] ) / 2.
                #print "Beam:", beam,
                comp_flux *= beam

            # calculate the average weighted patch center, and patch flux
            patch_flux  += comp_flux
            patch_ra    += comp_ra * comp_flux
            patch_dec   += comp_dec * comp_flux
            weights_ra  += comp_flux
            weights_dec += comp_flux

        print 'Center RA', patch_ra/weights_ra,
        print 'Center DEC', patch_dec/weights_dec,
        print 'Flux', patch_flux

        patches.append(Patch(patch_name, patch_ra/weights_ra, patch_dec/weights_dec, patch_flux))

    return patches


def create_clusters(patches, Q, show_plot):
    """
    Clusterize all the patches of the skymodel iteratively around the brightest patches
    """

    # sort the patches by brightest first
    idx = np.argsort([patch.flux for patch in patches])[::-1] # -1 to reverse sort
    patches = list(np.array(patches)[idx])

    # initialize clusters with the brightest patches
    clusters = []
    for i, patch in enumerate(patches[0:Q]):
        clusters.append(Cluster('CLUSTER_'+str(i), patch))

    # Iterate until no changes in which patch belongs to which cluster
    count = 1
    patches_seq_old = []
    while True:

        for cluster in clusters:
            cluster.update_cluster_coords()
            # reset patches
            #print cluster.name, ':',
            #print [patch.name for patch in cluster.patches]
            cluster.patches = []
    
        # select patches after the first Q
        for patch in patches[Q:]:
            # finding closest cluster for each patch
            angulardistance = np.inf
            for cluster in clusters:
                adis = abs(angsep2(cluster.centroid_ra, cluster.centroid_dec, patch.ra, patch.dec))
                if adis < angulardistance:
                    angulardistance = adis
                    closest_cluster = cluster
            # add this patch to the closest cluster
            closest_cluster.add_patch(patch)

        patches_seq = []
        for cluster in clusters:
            patches_seq.extend(cluster.patches)

        print 'Iteration', count
        count += 1
        if patches_seq == patches_seq_old: break
        patches_seq_old = patches_seq

    if show_plot:
        # check out http://phrogz.net/css/distinct-colors.html
        color=itertools.cycle(["#d96c6c", "#b23000", "#402310", "#ffd9bf", "#998273", "#b25f00", "#ffaa00", "#736039",
                               "#332f26", "#fff240", "#838c00", "#d5d9a3", "#aaff00", "#334d00", "#657356", "#bfffc8",
                               "#00e65c", "#008c38", "#00ffee", "#004d47", "#59adb3", "#00c2f2", "#267399", "#738c99",
                               "#0066ff", "#79aaf2", "#264599", "#0000ff", "#0000cc", "#000033", "#323040", "#220080",
                               "#b866cc", "#9d7ca6", "#ee00ff", "#59164c", "#ff40a6", "#330d21", "#8c4662", "#ff0022",
                               "#73000f", "#ffbfc8"])
        for cluster in clusters:
            c = color.next()

            # plot patches
            for patch in cluster.patches:
                matplotlib.pyplot.plot(patch.ra,patch.dec,'o',color=c,markersize=10)

            # plot clusters
            matplotlib.pyplot.plot(cluster.centroid_ra,cluster.centroid_dec,'*',color=c,markersize=25)
            matplotlib.pyplot.xlabel('RA [deg]')
            matplotlib.pyplot.ylabel('DEC [deg]')

            matplotlib.pyplot.text(cluster.centroid_ra+0.05, cluster.centroid_dec+0.05, str(cluster.name))

        # Reverse RA axis because of RA definition
        ax = matplotlib.pyplot.gca()
        ax.invert_xaxis()

        matplotlib.pyplot.show()

    return clusters, patches


def write_skymodel(clusters, data, outfilename, to_patch):
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

    for cluster in clusters:

        if to_patch:
            outfile.write(", , " + cluster.name + ", " + ratohms_string(cluster.mid_ra) + ", " + dectodms_string(cluster.mid_dec) + "\n")
        outfile.write("# Cluster centroid: " + ratohms_string(cluster.centroid_ra) + ", " + dectodms_string(cluster.centroid_dec) + "\n")
        outfile.write("# Total cluster flux="+str(round(cluster.total_cluster_flux(),4))+ " Jy\n")


        ra_size = (np.max([patch.ra for patch in cluster.patches])-np.min([patch.ra for patch in cluster.patches]))*3600.
        dec_size = (np.max([patch.dec for patch in cluster.patches])-np.min([patch.dec for patch in cluster.patches]))*3600.*np.cos(np.pi*cluster.mid_dec/180.)

        outfile.write("# RAsize - DECsize [arcsec] = " + str(np.int(ra_size)) + " - " + str(np.int(dec_size)) + "\n")

        # write patches defining cluster first
        for patch in [cluster.init_patch]+cluster.patches:
            idx = np.where(data['Patch'] == patch.name)[0]
            if not to_patch:
                outfile.write(", , " + cluster.name + "_" + patch.name + ", " + ratohms_string(cluster.mid_ra) + ", " + dectodms_string(cluster.mid_dec) + "\n")
            for el in idx:
                if to_patch:
                    linep1 = data['Name'][el]+', '+data['Type'][el]+', '+cluster.name+', '+ratohms_string(data['Ra'][el])+', '+dectodms_string(data['Dec'][el])+', '
                else:
                    linep1 = data['Name'][el]+', '+data['Type'][el]+', '+cluster.name+'_'+patch.name+', '+ratohms_string(data['Ra'][el])+', '+dectodms_string(data['Dec'][el])+', '

                linep2 = str(data['I'][el])+', '+str(data['Q'][el]) +', '+str(data['U'][el]) +', '+str(data['V'][el])+', '
                linep3 = str(data['MajorAxis'][el])+', '+str(data['MinorAxis'][el]) +', '+str(data['Orientation'][el]) +', '
                linep4 = str(data['ReferenceFrequency'][el])+', '+ data['SpectralIndex'][el]
                outfile.write(linep1+linep2+linep3+linep4+'\n')

        outfile.write("\n")
        outfile.write("\n")

    outfile.close()


if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <PyBDSM BBS catalog>')
    parser.add_option('-o', '--outfile', help='Name of the output skymodel file [default catalogue.clusters.skymodel]', metavar='VAL', default='catalog.clusters.skymodel')
    parser.add_option('-n', '--numclusters', dest='Q', help='Number of clusters [default: 10]', metavar='#', default=10, type=int)
    parser.add_option('-s', action='store_true', dest='show_plot', help='Show plot of the cluster disposition', default=False)
    parser.add_option('-p', action='store_true', dest='to_patch', help='Save clusters as patches insted of prefix [default=False]', default=False)
    parser.add_option('-b', '--beamMS', help='Give an MS used to find freq/dir/ant and correct for beam effect when clustering [default=None]', default=None)
    (o, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit()

    data = read_skymodel(args[0])
    patches = compute_patch_center(data, o.beamMS)
    clusters, patches = create_clusters(patches, o.Q, o.show_plot)
    write_skymodel(clusters, data, o.outfile, o.to_patch)

    print "Done."
