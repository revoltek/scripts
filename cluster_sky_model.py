#!/usr/bin/python
#
# Copyright (C) 2013 - Reinout van Weeren
# Copyright (C) 2011 - Francesco de Gasperin
# 
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
#
import matplotlib
matplotlib.use('GTK')
import numpy as np
import pylab as pl
import matplotlib, itertools
from lib_coordinates_mode import *
import os
import sys
pi = np.pi

def load_bbs_skymodel(infilename):
        tmp_input = infilename + '.tmp'
	
	# remove empty lines sed '/^$/d'
	# remove format line grep -v 'format'
	# remove comment lines  grep -v '#'
        os.system("grep -v '00:00:00, +00.00.00' SB200.skymodel | grep -v '#' | grep -v 'format' | sed '/^$/d'>" + tmp_input) # to remove patches headers from skymodel

        types = np.dtype({'names':['Name', 'Type','Patch','Ra', 'Dec', 'I', 'Q', 'U', 'V', 'Maj', 'Min', 'PA', 'RefFreq', 'Spidx'],\
	       'formats':['S100','S100','S100','S100','S100',np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,'S100']})

	    data  = np.loadtxt(tmp_input, comments='format', unpack=True, delimiter=', ', dtype=types)
	    os.system('rm ' + tmp_input)
        return data	


def compute_patch_center(data):
   print 'These are the input patches', numpy.unique(data['Patch'])
   patches = numpy.unique(data['Patch'])
   ra_patches   = numpy.zeros(len(patches))
   dec_patches  = numpy.zeros(len(patches))
   flux_patches = numpy.zeros(len(patches))

   for (patch_id,patch) in enumerate(patches):
    idx = numpy.where(data['Patch'] == patch)
    print 'Patch', patch, 'has', len(idx[0]), 'components'
    #print numpy.shape( data['Patch'][idx])
  
    # set to zero
    ra_patch   = 0.
    dec_patch  = 0.
    ra_weights = 0.
    dec_weights= 0.
    flux_patch = 0.
    
    for component in idx[0]:
    
      # conver RA, DEC to degrees for component
      ra_comp  = data['Ra'][component]
      dec_comp = data['Dec'][component]
      ra_comp  = (ra_comp.split(':'))
      dec_comp = (dec_comp.split('.'))
      flux_comp= numpy.float(data['I'][component])
      ra_comp  = hmstora(numpy.float(ra_comp[0]),numpy.float(ra_comp[1]),numpy.float(ra_comp[2]))
 
      if len(dec_comp) == 4:  # decimal arcsec in Dec
        dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(str(dec_comp[2]+"."+dec_comp[3]))) 
      else:
        dec_comp = dmstodec(numpy.float(dec_comp[0]),numpy.float(dec_comp[1]),numpy.float(dec_comp[2]))
       
      # calculate the average weighted patch center, and patch flux 
      flux_patch = flux_patch + flux_comp
      #print ra_comp, dec_comp, flux_comp
      ra_patch = ra_patch + (flux_comp*ra_comp)
      dec_patch= dec_patch+ (flux_comp*dec_comp)
      ra_weights = ra_weights + flux_comp
      dec_weights= dec_weights + flux_comp
    
    print 'Center RA, Center DEC, flux', ra_patch/ra_weights, dec_patch/dec_weights,  flux_patch
    
    ra_patches[patch_id]  = ra_patch/ra_weights
    dec_patches[patch_id] = dec_patch/dec_weights
    flux_patches[patch_id]= flux_patch
    
   return patches, ra_patches,dec_patches, flux_patches
	
	


def create_initial_clusters(patches,ra_patches, dec_patches, flux_patches, Q, show_plot):
  
  # sort the patches by brightest first
  idx = numpy.argsort(flux_patches)[::-1] #.sort()
  flux_patches = flux_patches[idx]
  ra_patches   = ra_patches[idx]
  dec_patches  = dec_patches[idx]
  patches      = patches[idx]
  #data         = data[idx]
   

  patch_cluster_id = numpy.zeros(len(patches[Q:]), dtype=int)
 
  # find the Q brightest patches
  clusters = numpy.copy(patches[0:Q])
  clusters_ra = numpy.copy(ra_patches[0:Q])
  clusters_dec = numpy.copy(dec_patches[0:Q])
  clusters_flux= numpy.copy(flux_patches[0:Q])

  print clusters 
   
   
  # find the closest cluster to a given patch
  patch_cluster_id = find_closest_cluster(clusters,clusters_ra,clusters_dec,\
                     patches,ra_patches,dec_patches,Q)
  #print patch_cluster_id
 

  clusters_ra_new,clusters_dec_new = update_cluster_centers(clusters,clusters_ra,\
                                     clusters_dec,clusters_flux,ra_patches,      \
  				     dec_patches,flux_patches,patch_cluster_id,Q)
				     
  patch_cluster_id_new = find_closest_cluster(clusters,clusters_ra_new,clusters_dec_new,\
                           patches,ra_patches,dec_patches,Q)
  count = 1
  while (numpy.array_equal(patch_cluster_id,patch_cluster_id_new)) == False:
  
    clusters_ra_new,clusters_dec_new = numpy.copy(update_cluster_centers(clusters,clusters_ra_new,\
                                     clusters_dec_new,clusters_flux,ra_patches,      \
  				     dec_patches,flux_patches,patch_cluster_id,Q))
    patch_cluster_id = numpy.copy(patch_cluster_id_new)			      	     
    patch_cluster_id_new = find_closest_cluster(clusters,clusters_ra_new,clusters_dec_new,\
                           patches,ra_patches,dec_patches,Q)
    print 'Iteration', count
    count = count+1 
  
  print clusters_flux
  clusters_flux = compute_total_cluster_flux(clusters_flux,flux_patches,patch_cluster_id_new,Q) 
  print clusters_flux

  #show_plot = True
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
      idx = numpy.where(patch_cluster_id == cluster_id)
      for el in range(len(idx[0])):
         matplotlib.pyplot.plot(ra_patches[idx[0][el]+Q],dec_patches[idx[0][el]+Q],'o',color=c,markersize=5)
    
      # plot clusters  
      matplotlib.pyplot.plot(clusters_ra_new[cluster_id],clusters_dec_new[cluster_id],'*',color=c,markersize=15)
      matplotlib.pyplot.xlabel('RA [deg]')
      matplotlib.pyplot.ylabel('DEC [deg]')
      
      matplotlib.pyplot.text(clusters_ra_new[cluster_id]+0.05, clusters_dec_new[cluster_id]+0.05, str(cluster_id))
      #matplotlib.pyplot.text(clusters_ra_new[cluster_id]+0.05, clusters_dec_new[cluster_id]+0.05, str(cluster))
      
    # Reverse RA axis because of RA definition
    # REVERSE THE X-AXIS
    ax = matplotlib.pyplot.gca()
    ax.invert_xaxis()
      
    matplotlib.pyplot.show()
 

    ## clusters_ra_new are the centroids, note return new values (they do not match anymore with the defining patch!!)
     
  return clusters,clusters_ra_new,clusters_dec_new,clusters_flux,patches,ra_patches,dec_patches,patch_cluster_id_new




def write_skymodel(clusters,clusters_ra,clusters_dec,clusters_flux,patches,patch_cluster_id_new,\
                   data,Q,outfilename,cluster_ra_min_max, cluster_dec_min_max) :

 outfile =open(outfilename,'w')
 outfile.write("format = Name, Type, Patch, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, Orientation, ReferenceFrequency='1.49803e+08', SpectralIndex='[]'\n")
 outfile.write("\n")
 outfile.write("\n")
 outfile.write("## updated clustered skymodel\n")
 outfile.write("\n")
 outfile.write("\n") 
 
 for cluster_id,cluster in enumerate(clusters):

   clustername = 'CLUSTER_' + str(cluster_id)

   # compute RA, DEC patch in sexagesimal format (from degr)
   mid_ra = numpy.mean(cluster_ra_min_max[cluster_id,:])
   mid_dec= numpy.mean(cluster_dec_min_max[cluster_id,:])
   
   ihr = int(mid_ra/15.)
   xmin =abs(mid_ra*4.0-ihr*60.0)
   imin = int(xmin)
   xsec = (xmin-imin)*60.0
   #print 'RA', ihr, imin, xsec
   ideg = int(mid_dec)
   xmn = abs(mid_dec-ideg)*60.0
   imn = int(xmn)
   xsc = (xmn-imn)*60.0
   #print 'DEC', ideg, imn, xsc

   cl_ra_str = str(ihr) + ':' + str(imin) + ':' +str(round(xsec,2))
   cl_dec_str= str(ideg)+ '.' + str(imn)  + '.' +str(round(xsc,1))

   outfile.write(", , " + clustername + ", "+ cl_ra_str + ", " + cl_dec_str+"\n") 
   outfile.write("# NDPPP patch coordinates: " + str(ihr) + 'h' + str(imin) + 'm' +str(round(xsec,2))+', '+\
                 str(ideg)+ 'd' + str(imn)  + 'm' +str(round(xsc,1))+"\n")
   outfile.write("# Total cluster flux="+str(round(clusters_flux[cluster_id],4))+ " Jy\n")
   
   dec_size = (numpy.max(cluster_dec_min_max[cluster_id,:])-numpy.min(cluster_dec_min_max[cluster_id,:]))*3600.
   ra_size =  (numpy.max(cluster_ra_min_max[cluster_id,:])-numpy.min(cluster_ra_min_max[cluster_id,:]))*3600.*numpy.cos(pi*mid_dec/180.)
   
   
   mmstr = 'RAsize,DECsize [arcsec] = '+ str([numpy.int(ra_size),numpy.int(dec_size)])
   
   outfile.write("# "+ mmstr + "\n")
   
   # write paatches defining cluster first
   idx = numpy.where(data['Patch'] == cluster)
   for el in (idx[0]):
   
     linep1 = data['Name'][el]+', '+data['Type'][el]+', '+clustername+', '+data['Ra'][el]+', '+data['Dec'][el]+', '
     linep2 = str(data['I'][el])+', '+str(data['Q'][el]) +', '+str(data['U'][el]) +', '+str(data['V'][el])+', '
     linep3 = str(data['Maj'][el])+', '+str(data['Min'][el]) +', '+str(data['PA'][el]) +', '
     linep4 = str(data['RefFreq'][el])+', '+ data['Spidx'][el]
     outfile.write(linep1+linep2+linep3+linep4 +'\n')
     
   
   
   
   #write other patches in this cluster
   idx = numpy.where(patch_cluster_id_new == cluster_id)
   patches_list = patches[idx[0]+Q]
   for patch in patches_list:
    
     idx = numpy.where(data['Patch'] == patch)
     for el in (idx[0]):
       linep1 = data['Name'][el]+', '+data['Type'][el]+', '+clustername+', '+data['Ra'][el]+', '+data['Dec'][el]+', '
       linep2 = str(data['I'][el])+', '+str(data['Q'][el]) +', '+str(data['U'][el]) +', '+str(data['V'][el])+', '
       linep3 = str(data['Maj'][el])+', '+str(data['Min'][el]) +', '+str(data['PA'][el]) +', '
       linep4 = str(data['RefFreq'][el])+', '+ data['Spidx'][el]
       outfile.write(linep1+linep2+linep3+linep4 +'\n')
   
   
   
   
   
   outfile.write("\n") 
   outfile.write("\n") 
    
 outfile.close()







def compute_total_cluster_flux(clusters_flux,flux_patches,patch_cluster_id_new,Q):
  clusters_flux_new = numpy.copy(clusters_flux) # copy into new array

  for cluster_id in range(len(clusters_flux)):
    idx = numpy.where(patch_cluster_id_new == cluster_id)
     
    # add the flux from the mathcing patches 
    for el in range(len(idx[0])): 
      clusters_flux_new[cluster_id] = clusters_flux_new[cluster_id] + flux_patches[idx[0][el]+Q]

  return clusters_flux_new




def find_closest_cluster(clusters,clusters_ra,clusters_dec,patches,ra_patches,dec_patches,Q):  
  patch_cluster_id = numpy.zeros(len(patches[Q:]), dtype=int)
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
  clusters_ra_new  = numpy.copy(clusters_ra) * 0.
  clusters_dec_new =  numpy.copy(clusters_dec)* 0.
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
     idx = numpy.where(patch_cluster_id == cluster_id)
     for el in range(len(idx[0])):
       cl_ra      = cl_ra  + ra_patches[idx[0][el]+Q]*flux_patches[idx[0][el]+Q]
       cl_dec     = cl_dec + dec_patches[idx[0][el]+Q]*flux_patches[idx[0][el]+Q]
       weight_ra  = weight_ra + flux_patches[idx[0][el]+Q]
       weight_dec = weight_dec+ flux_patches[idx[0][el]+Q]
  
     cl_ra = cl_ra/weight_ra
     cl_dec= cl_dec/weight_dec
     #print 'new cluster coordinate', cl_ra, cl_dec
     clusters_ra_new[cluster_id] = cl_ra
     clusters_dec_new[cluster_id] = cl_dec

  return clusters_ra_new,clusters_dec_new


def compute_cluster_info(clusters,clusters_ra,clusters_dec,patches,ra_patches,dec_patches,patch_cluster_id_new,Q):
  cluster_ra_min_max  = numpy.zeros((len(clusters),2))
  cluster_dec_min_max = numpy.zeros((len(clusters),2))


  for cluster_id,cluster in enumerate(clusters):

    idx = numpy.copy(numpy.where(patch_cluster_id_new == cluster_id))
    idx = idx[0] + Q

    idx = numpy.sort(numpy.append(idx, cluster_id))
    
    ra_min = numpy.min(ra_patches[idx])
    ra_max = numpy.max(ra_patches[idx])
    dec_min= numpy.min(dec_patches[idx])
    dec_max= numpy.max(dec_patches[idx])
    
    cluster_ra_min_max[cluster_id,:] = [ra_min,ra_max]
    cluster_dec_min_max[cluster_id,:] = [dec_min,dec_max]

  return cluster_ra_min_max, cluster_dec_min_max

		
if __name__ == '__main__':
    import sys
    argc=len(sys.argv)
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <PyBDSM BBS catalog>')
    parser.add_option('--outfile', dest='file_clusters', help='Name of the output skymodel file [default catalogue.clusters.skymodel]', \
                      metavar='VAL', default='catalog.clusters.skymodel')
    parser.add_option('--numclusters', dest='Q', help='Number of clusters [default: 10]', metavar='VAL', default=10)
    parser.add_option('-p', action='store_true', dest='show_plot', help='Show plot of the cluster disposition', default=False)
    (o, args) = parser.parse_args()

    if len(args) < 1: sys.exit("Missing BBS-format sky model.")
    
    data = load_bbs_skymodel(args[0])
    patches,ra_patches,dec_patches, flux_patches =  compute_patch_center(data)

    clusters,clusters_ra,clusters_dec,clusters_flux,patches,ra_patches,dec_patches,patch_cluster_id_new = \
               create_initial_clusters(patches, ra_patches, dec_patches, flux_patches, numpy.int(o.Q), o.show_plot)
   
    cluster_ra_min_max, cluster_dec_min_max = compute_cluster_info(clusters,clusters_ra,clusters_dec,patches,ra_patches, dec_patches,patch_cluster_id_new, numpy.int(o.Q))
   
    write_skymodel(clusters,clusters_ra,clusters_dec,clusters_flux,patches,patch_cluster_id_new,data,numpy.int(o.Q),o.file_clusters, cluster_ra_min_max, cluster_dec_min_max)
    

    print "Done."
