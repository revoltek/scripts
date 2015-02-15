def correctPB(imgname, freq=0, phaseCentre=None):
    """Given an image "img" create the "img.pbcorr" which
    has each pixel corrected for the WSRT primary beam effect.
    phaseCentre: [ra,dec] in deg of the pointing direction
    """
    import numpy as np
    img = ia.open(imgname)
    cs = ia.coordsys()
    if freq == 0: freq = cs.restfrequency()['value'][0]

    print "Frequency is", freq/1.e6, "MHz"

    # if not specified assuming pointing in the centre of the image
    if phaseCentre == None:
        print "Assume pointing in the image centre."
        pixPhaseCentre = ia.topixel( () )['numeric'][0:2]
    else: 
        pixPhaseCentre = ia.topixel( qa.quantity(str(phaseCentre[0])+'deg'), qa.quantity(str(phaseCentre[1])+'deg') )['numeric'][0:2]
        print "Phase centre is at pix: ", pixPhaseCentre

    # function to initialize the beam-array
    assert abs(cs.increment()['numeric'][0]) == abs(cs.increment()['numeric'][1])
    pix2deg = abs(cs.increment()['numeric'][0])*180./np.pi # increment is in rad
    def beam_creator(i,j):
        # distance from phase centre pixel in deg
        d = np.sqrt( (pixPhaseCentre[0] - i)**2 + (pixPhaseCentre[1] - j)**2 ) * pix2deg
        c = 0.68
        nu = freq/1.e9 # freq in GHz
        return np.cos(c*nu*d)**6

    beam = np.fromfunction(beam_creator, ia.shape()[0:2])

    # write new image
    impbcor(imagename=imgname, pbimage=beam, outfile=imgname+'.pbcorr', mode='divide', overwrite=True)
