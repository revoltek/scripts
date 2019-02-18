#!/usr/bin/env python

"""
A script to create a mosaic of LOFAR images
(The script was written to generate MSSS mosaics,
but it can be used in other cases as well)

HISTORY
=======
v0.1	G. Heald	Created script
v0.2	R. Breton	Fixed avgpb behavior and weighting implementation
v0.3	G. Heald	Add weights, fix image naming, output sensitivity map
v0.4    S. van Velzen   Add beam & frequency info to header
v0.5    S. van Velzen   Use pyrap to save fits files
v0.6    J. Swinbank     Select a single Stokes parameter from input maps
v0.7	G. Heald	Fix RA behavior, add NCP option, pyfits tweak
v0.8	G. Heald	Fix behavior near RA=0
v0.9    G. Heald        Fix pyfits behavior (updating beam info)

TO DO
=====
Add NCP option (flag added in v0.7 but does not do anything)

"""

version = '0.9 2013-08-14'

import pyrap.images as pim
from pyrap import quanta
import numpy as np
import argparse
#import pylab as plt
import pyfits
import os
import time

def main(args):

 # Generate lists of input images and check that they exist
 images=[]
 avgpbs=[]
 psf_fwhm = [] # resolution
 frequency = [] # frequency of images (should be equal?)
 imagesbase=args.images.split(',')
 for base in imagesbase:
  images.append(base+'.'+args.extension)
  if not os.path.exists(images[-1]):
   print("Error: Image",images[-1],"does not exist")
   return 1
  avgpbs.append(base+'.'+args.avgpbext)
  if not os.path.exists(avgpbs[-1]):
   print("Error: PB image",avgpbs[-1],"does not exist")
   return 1

 # Collect weights and invert if requested
 if args.weights == '':
  weights = np.ones(len(images))
 else:
  weights = np.array(args.weights.split(',')).astype('float')
  if args.invertwt:
   weights = 1./weights
  if len(weights) != len(images):
   print("Error: List of weights is not the same length as list of images.")
   return 1
 print("Combining images")
 formstr = '{0:45s}  {1:45s} {2:s}  {3:s} {4:s} {5:s}'
 print(formstr.format("-----","--------","------------","-------","-------","------"))
 print(formstr.format("Image", "PB image","Norm. weight", "Maj(ac)", "Min(ac)","PA(deg)"))
 print(formstr.format("-----","--------","------------","-------","-------","------"))

 for i in range(len(images)):
  this_pim = pim.image(images[i])
  info_dict = this_pim.info()['imageinfo']['restoringbeam']
  # get beam info
  bpar_ma = quanta.quantity(info_dict['major']).get_value('deg')
  bpar_mi = quanta.quantity(info_dict['minor']).get_value('deg')
  bpar_pa = quanta.quantity(info_dict['positionangle']).get_value('deg')
  psf_fwhm.append([bpar_ma, bpar_mi, bpar_pa])
  frequency.append(this_pim.info()['coordinates']['spectral2']['restfreq'])
  print('{0:45.45s}  {1:45.45s} {2:0.2f}          {3:0.2f}    {4:0.2f}    {5:0.2f}'.format(images[i], avgpbs[i], weights[i]/sum(weights), bpar_ma*60, bpar_mi*60,bpar_pa))

 psf_fwhm = np.array(psf_fwhm)
 frequency = np.array(frequency)
 mean_psf_fwhm = np.mean(psf_fwhm, axis=0)
 mean_frequency = np.mean(frequency)
 print('\nmean Beam: {0:0.3f} maj (arcmin), {1:2.3f} min (arcmin), {2:0.2f} pa (deg)'.format(mean_psf_fwhm[0]*60, mean_psf_fwhm[1]*60, mean_psf_fwhm[2]))
 print('(Frequency (MHz):', mean_frequency*1e-6)

 if np.max(mean_frequency-frequency)/mean_frequency > 1e-6:
  print('\n\nWARNING.\nAre you using  images from different bands?')
  print('Frequencies (Hz):', frequency)

 # Initialize some vectors
 declims = [] # store the limits of the declination axes
 ralims = [] # store the limits of the r.a. axes
 rainc = [] # store the r.a. increments in case they differ
 decinc = [] # store the dec increments in case they differ
 pims = [] # stores the pyrap images of the data
 ppbs = [] # stores the pyrap images of the pb images


# Get image frames for input images
 for im, pb in zip(images, avgpbs):
  image = pim.image(im)
  sptcoords = image.coordinates().get_coordinate('spectral')
  nc = sptcoords.get_axis_size()
  assert(sptcoords.get_image_axis() == 0)

  # Get Stokes axis. Ensure we are working with the Stokes parameter requested.
  stkcoords = image.coordinates().get_coordinate('stokes')
  assert(stkcoords.get_image_axis() == 1)
  if stkcoords.get_axis_size() == 1:
   assert(stkcoords.get_stokes()[0] == args.stokes)
  else:
   stks = stkcoords.get_stokes().index(args.stokes)
   image = image.subimage(blc=(0, stks), trc=(nc-1, stks), dropdegenerate=False)
  ns = 1

  dircoords = image.coordinates().get_coordinate('direction')
  nx = dircoords.get_axis_size(axis=1)
  ny = dircoords.get_axis_size(axis=0)
  inc = dircoords.get_increment()
  ref = dircoords.get_referencepixel()
  val = dircoords.get_referencevalue()
  ra_axis = (list(range(nx))-ref[1])*inc[1]+val[1]
  dec_axis = (list(range(ny))-ref[0])*inc[0]+val[0]
  rainc.append(inc[1])
  decinc.append(inc[0])
  declims.append(min(dec_axis))
  declims.append(max(dec_axis))
  mean_ra = np.mean(ra_axis)
  ralims.append((min(ra_axis)-mean_ra)*np.cos(val[0])+mean_ra)
  ralims.append((max(ra_axis)-mean_ra)*np.cos(val[0])+mean_ra)
  pims.append(image)
  ppbs.append(pim.image(pb))


 # Generate the mosaic coordinate frame
 master_dec = np.arange(min(declims),max(declims),min(decinc))
 if max(ralims)-min(ralims) > 5.*np.pi/3.: # crossed RA=0
  print("Warning: I think the mosaic crosses RA=0, treating the coordinates as such.")
  #ralims[ralims>np.pi] -= 2.*np.pi
  for i in range(len(ralims)):
   if ralims[i]>np.pi: ralims[i] = ralims[i]-2.*np.pi
 master_ra = np.arange(max(ralims),min(ralims),max(rainc))
 if args.verbose:
  print("Found ra,dec pixel increments (arcsec):")
  print(np.array(rainc)*206265.)
  print(np.array(decinc)*206265.)
 ma = pims[-1].coordinates()
 ma['direction'].set_referencepixel([len(master_dec)/2,len(master_ra)/2])
 ma['direction'].set_increment([min(decinc),max(rainc)])
 ma['direction'].set_referencevalue([master_dec[len(master_dec)/2],master_ra[len(master_ra)/2]])
 if args.NCP:
  print('Setting NCP projection is not yet working ....')
  #ma['direction'].set_projection('ZEA')

 # Initialize the arrays for the output image, sensitivity, and weights
 master_im = np.zeros((len(master_dec),len(master_ra)))
 master_weight = np.zeros((len(master_dec),len(master_ra)))
 master_sens = np.zeros((len(master_dec),len(master_ra)))

 # Reproject the images onto the master grid, weight and normalize
 for i in range(len(pims)):
  im = pims[i].regrid([2,3],ma,outshape=(nc,ns,len(master_dec),len(master_ra)))
  pb = ppbs[i].regrid([2,3],ma,outshape=(nc,ns,len(master_dec),len(master_ra)))
  imdata = np.squeeze(im.getdata())
  pbdata = np.squeeze(pb.getdata())
  newim = imdata
  newpb = pbdata
  newwt = (weights[i]*newpb)**2
  master_im += newim*newwt
  master_sens += newpb*newwt
  master_weight += newwt
 inds = master_weight != 0.
 master_im[inds] /= master_weight[inds]
 master_sens[inds] /= master_weight[inds]

 # Show image if requested
 if args.plotimg:
  plt.imshow(master_im,vmin=0.,vmax=0.5)
  plt.show()

 # Write fits files
 arrax = np.zeros( (1,1, len(master_im[:,0]), len(master_im[0,:])) )
 arrax[0,0,:,:] = master_im


 # Open new casa image for mosaic
 new_pim = pim.image('',shape=(1,1, len(master_dec),len(master_ra)), coordsys=ma)
 new_pim.putdata(arrax)
 # Write fits
 new_pim.tofits(args.outfits, overwrite=True)
 # Same for sensitivity
 new_pim_sens = pim.image('',shape=(1,1,len(master_dec),len(master_ra)),coordsys=ma)
 arrax[0,0,:,:] = master_sens
 new_pim_sens.putdata(arrax)  #
 new_pim_sens.tofits(args.sensfits, overwrite=True)

 # need to add new beam info (not sure if this is possible with pyrap)
 hdu = pyfits.open(args.outfits,mode='update')
 header = hdu[0].header
 header.update('BMAJ',mean_psf_fwhm[0])
 header.update('BMIN',mean_psf_fwhm[1])
 header.update('BPA',mean_psf_fwhm[2])
 header.update('BUNIT',pims[-1].info()['unit'])
 header.update('RESTFRQ',mean_frequency)
 header.update('RESTFREQ',mean_frequency)
 newhdu = pyfits.PrimaryHDU(data=hdu[0].data, header=header)
 newhdu.writeto(args.outfits,clobber=True)

 return

print("LOFAR mosaic generator, v"+version+'\n')
parser = argparse.ArgumentParser(description="Mosaic MSSS images.")
parser.add_argument('images',metavar='list_of_images',help='Input image names without extension, as a comma separated list (no spaces). Note that the "0" in the avgpb filenames will be accounted for in the script, do not include it here.')
parser.add_argument('-x','--extension',help='Image extension to combine [default restored.corr]',default='restored.corr')
parser.add_argument('-a','--avgpbext',help='Extension for primary beam images [default avgpb]',default='avgpb')
parser.add_argument('-w','--weights',help='Image weights, as a comma separated list (this must be of the same length and in same order as the IMAGES list) [default equal weights]',default='')
parser.add_argument('-n','--invertwt',help='Invert weights before applying? (Can be useful if image noise values are given as input to WEIGHTS) [default False]',action='store_true',default=False)
parser.add_argument('-v','--verbose',help='Give some verbose output [default False]',action='store_true',default=False)
parser.add_argument('-o','--outfits',help='Output name of mosaic fits file [default mosaic.fits]',default='mosaic.fits')
parser.add_argument('-s','--sensfits',help='Output name of sensitivity fits file [default sensitivity.fits]',default='sensitivity.fits')
parser.add_argument('-p','--plotimg',help='Display image on screen? [default False]',action='store_true',default=False)
parser.add_argument('-S','--stokes',help='Stokes parameter to use?  [default I]',default='I')
parser.add_argument('-N','--NCP',help='Use NCP instead of SIN? This option does not work yet. [default False]',default=False,action='store_true')
args = parser.parse_args()
main(args)
