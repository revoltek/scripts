def correctPB(imgname, freq=0, phaseCentre=None):
    """Given an image "img" create the "img.pbcorr" which
    has each pixel corrected for the GMRT primary beam effect.
    phaseCentre: [ra,dec] in deg of the pointing direction
    """
    import numpy as np
    img = ia.open(imgname)
    cs = ia.coordsys()
    if freq == 0: freq = cs.restfrequency()['value'][0]

    # find the correct freq
    freq = min([153,235,325,610,1400], key=lambda x:abs(x-freq/1.e6))
    print "Frequency is", freq, "MHz"

    # from http://gmrt.ncra.tifr.res.in/gmrt_hpage/Users/doc/manual/UsersManual/node27.html
    parm = {153: [-4.04,76.2,-68.8,22.03],
    235: [-3.366,46.159,-29.963,7.529],
    325: [-3.397,47.192,-30.931,7.803],
    610: [-3.486,47.749,-35.203,10.399],
    1400: [-2.27961,21.4611,-9.7929,1.80153]}[freq]

    # if not specified assuming pointing in the centre of the image
    if phaseCentre == None:
        pixPhaseCentre = ia.topixel( () )['numeric'][0:2]
    else: 
        pixPhaseCentre = ia.topixel( qa.quantity(str(phaseCentre[0])+'deg'), qa.quantity(str(phaseCentre[1])+'deg') )['numeric'][0:2]
        print "Phase centre is at pix: ", pixPhaseCentre

    # function to initialize the beam-array
    assert abs(cs.increment()['numeric'][0]) == abs(cs.increment()['numeric'][1])
    pix2deg = abs(cs.increment()['numeric'][0])*180./np.pi # increment is in rad
    def beam_creator(i,j):
        # get distance from phase centre pixel in deg
        d = np.sqrt( (pixPhaseCentre[0] - i)**2 + (pixPhaseCentre[1] - j)**2 ) * pix2deg
        # from http://www.aips.nrao.edu/cgi-bin/ZXHLP2.PL?PBCOR (converto to arcmin and multiply by freq in GHz)
        d = d * 60 * freq/1.e3
        return 1 + (parm[0]/10**3)*d**2 + (parm[1]/10**7)*d**4 + \
             (parm[2]/10**10)*d**6 + (parm[3]/10**13)*d**8 

    beam = np.fromfunction(beam_creator, ia.shape()[0:2])

    # write new image
    impbcor(imagename=imgname, pbimage=beam, outfile=imgname+'.pbcorr', mode='divide', overwrite=True)
