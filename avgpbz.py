#!/usr/bin/env python

"""
Script to zero the corners of avgpb images
"""

import pyrap.images as pim
import argparse
import numpy as np

def main(args):

 image = args.image
 pb = pim.image(image)
 pbdata = pb.getdata()
 (nx,ny) = pbdata.shape[2:]
 pbrad = args.radius*nx
 cy = ny/2
 cx = nx/2
 pbcounter = 0
 for j in range(nx,cx,-1):
  k = ny
  while np.sqrt((j-cx)**2+(k-cy)**2)>pbrad and k>cy:
   pbcounter += 1
   pbdata[:,:,k-1,j-1]=0.
   pbdata[:,:,k-1,nx-j]=0.
   pbdata[:,:,ny-k,j-1]=0.
   pbdata[:,:,ny-k,nx-j]=0.
   k -= 1
  if k == ny:
   print j,nx
   break
 print pbcounter,'x 4 =',pbcounter*4,'zeros replaced'
 if args.output == '':
  while image[-1]=='/':
   image=image[:-1]
  outim = image+'z'
 else:
  outim = args.output
 print 'Writing',outim
 pout = pim.image(outim,values=pbdata,coordsys=pb.coordinates())
 
print "AVGPB editor"
parser = argparse.ArgumentParser(description="Zero corners of avgpb images")
parser.add_argument('image',help="Image to adjust")
parser.add_argument('-o','--output',help="Output image name [default is to add a 'z' to the end of the input filename",default='')
parser.add_argument('-r','--radius',help='Radius beyond which to zero avgpb values (expressed as fraction of image width) [default 0.5 = half of image width, i.e. zero outside of a circle with diameter NAXIS1]',default=0.5,type=float)
args=parser.parse_args()
main(args)


