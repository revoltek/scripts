#!/usr/bin/python

import os, sys, glob
from astropy.io import fits

def make_beam(imageroot, outfile='beam.fits'):

    imagefiles = sorted(glob.glob(imageroot+'0*-image.fits'))
    beamfiles = sorted(glob.glob(imageroot+'0*-beam-I.fits'))

    # get weight
    weights = []
    for imagefile in imagefiles:
        with fits.open(imagefile) as imagefits:
            header = imagefits[0].header
            weights.append(header['WSCIMGWG'])

    # combine beams
    with fits.open(beamfiles[0]) as beamfits:
        beam = beamfits[0].data * weights[0]
        header = beamfits[0].header
    for i, beamfile in enumerate(beamfiles[1:]):
        with fits.open(beamfiles) as beamfits:
            beam += beamfits[0].data * weigths[i+1]

    beam /= np.sum(weights)
    fits.writeto(outfile, beam, header)

if __name__=='__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v|-V] image-root-name \n Francesco de Gasperin', version='1.0')
    opt.add_option('-o', '--outfile', help='Filename for the combined beam fits (default=beam.fits)', default='beam.fits')
    (options, args) = opt.parse_args()


    outfile = options.outfile

    make_beam(args, outfile)
