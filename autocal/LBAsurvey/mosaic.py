#!/usr/bin/python
# 1. Can mosaic facets or entire fields
# 2. In the first case it blanks outside facet region and no weights
# 3. In the second case, it estimates the weights
# 4. Can apply primary beam correction

import os, sys
import astropy.io.fits as pyfits
from lib_pipeline_img import *
from lib_pipeline import *

class Direction(object):
    def __init__(self, fitsfile, region):
        self.fitsfile = fitsfile
        self.region = region
        
def main(args):

    # make dirs
    check_rm('mosaic')
    os.path.mkdirs('mosaic')

    # collect directions
    directions = []
    for fitsfile in glob.glob('peel/ddcal*'):
        dirname = dirname.replace('peel/','')
        fitsfile = 'peel/'+dirname+'/facetM-facet-MFS-image.fits'
        region = 'region/'+dirname+'-facet.reg'
        directions.append(Direcrion(fitsfile, region))


    os.system

# regrid image at facet resolution using Andre' code

# regrid facets

if __name__ == '__main__':
    import argparse
    print "LOFAR LBA mosaic generator, v"+version+'\n'
    parser = argparse.ArgumentParser(description="Mosaic PiLL images.")
    parser.add_argument('images',metavar='list_of_images',help='Input image names without extension, as a comma separated list (no spaces).')
    

    parser.add_argument('-x','--extension',help='Image extension to combine [default restored.corr]',default='restored.corr')
    parser.add_argument('-a','--avgpbext',help='Extension for primary beam images [default avgpb]',default='avgpb')
    parser.add_argument('-w','--weights',help='Image weights, as a comma separated list (this must be of the same length and in same order as the IMAGES list) [default equal weights]',default='')
    parser.add_argument('-o','--outfits',help='Output name of mosaic fits file [default mosaic.fits]',default='mosaic.fits')
    parser.add_argument('-s','--sensfits',help='Output name of sensitivity fits file [default sensitivity.fits]',default='sensitivity.fits')
    args = parser.parse_args()
    main(args)
