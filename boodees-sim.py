"""
Simulate 1 VLA pointing of the BooDeeS survey in A+B+C configuration
"""

import sys, os, glob
from math import *
import numpy as np

direction='J2000 14h30m00.0s 34d00m00.0s'

def do_model(direction, size, typ='empty'):
    """
    Generate a model

    typ: empty 
    typ: random
    """
    print("Making model...")
    ia.fromshape("sky.model", [256,256,1,10], overwrite=True)
    cs=ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])
    cell_rad = qa.convert(qa.quantity("0.1arcsec"),"rad")['value']
    cs.setincrement([-cell_rad,cell_rad],'direction')
    cs.setreferencevalue([qa.convert("14h30m",'rad')['value'],qa.convert("34deg",'rad')['value']],type="direction")
    cs.setreferencevalue("1.5GHz",'spectral')
    cs.setincrement('0.5MHz','spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    if typ == 'random':
        cl.done()
        cl.addcomponent(dir=direction, flux=1.0, fluxunit='Jy', freq='1.5GHz', shape="Gaussian",
                    majoraxis="0.1arcmin", minoraxis='0.05arcmin', positionangle='45.0deg')
        ia.modify(cl.torecord(),subtract=False)
    exportfits(imagename='sky.model',fitsimage='sky.model.fits',overwrite=True)
    return 'sky.model'



def do_simulate(project, integration, direction, antennalist, totaltime='3600s', hourangle='transit'):
    """
    run the simulation
    """
    print("Simulating %s..." % project)
    simobserve(
           project      = project, 
           skymodel     = 'sky.model', 
           incenter     = '1.5GHz', 
           inwidth      = '60 MHz' , 
           setpointings = False, 
           ptgfile      = 'pointings.txt',
           integration  = integration,  
           direction    = direction,  
           obsmode      = 'int', 
           antennalist  = antennalist, 
           hourangle    = hourangle, 
           totaltime    = totaltime,  
           thermalnoise = 'tsys-atm', 
           graphics     = 'file', 
           overwrite    = True
           )

    listobs(vis=project+'/'+project+'.'+antennalist[:5]+'.noisy.ms')


def do_clean(weighting, briggs=0.5):
    """
    run the cleaning
    """
    if not os.path.exists('sim'): os.system('mkdir sim')

    print("Cleaning...")
    allvis = glob.glob('boodees*/boodees*.noisy.ms')
    clean( 
           vis          = allvis,
           imagename    = 'sim/%s' % weighting,
           niter        = 0,
           imsize       = [2048],
           cell         = '0.2arcsec',
           weighting    = weighting,
           robust       = briggs
           )

def do_stat(imagename):
    """
    print image statistics
    """
    # extract beam
    head = imhead(imagename=imagename, mode='summary')
    bmaj = head['restoringbeam']['major']['value']
    bmin = head['restoringbeam']['minor']['value']
    print('%s - Beam: %.1f" x %.1f"' % (imagename,bmaj,bmin))
    # extract noise
    stat = imstat(imagename=imagename)
    print('%s - RMS noise: %.1f uJy/b' % (imagename,stat['rms'][0]*1e6))


# create model
do_model(direction, '1deg')

# simulate
#hr = '-4h'
#do_simulate('boodees-A', '2s', direction, 'vla.a.cfg', '7200s', hr)
#do_simulate('boodees-B', '3s', direction, 'vla.b.cfg', '3600s', hr)
#do_simulate('boodees-C', '5s', direction, 'vla.c.cfg', '1800s', hr)

for hr in np.linspace(-6,0,6):
    hr = str(hr)+'h'
    do_simulate('boodees-A_'+hr, '2s', direction, 'vla.a.cfg', '1200s', hr)
    do_simulate('boodees-B_'+hr, '3s', direction, 'vla.b.cfg', '600s', hr)
    do_simulate('boodees-C_'+hr, '5s', direction, 'vla.c.cfg', '300s', hr)

# clean
#do_clean('natural')
do_clean('briggs', 0.)
#do_clean('uniform')

# stats
#do_stat('sim/natural.image')
do_stat('sim/briggs.image')
#do_stat('sim/uniform.image')


# PROBLEM
#for hr in np.linspace(-5,-3,6):
#sim/natural.image - Beam: 2.6" x 1.9"
#sim/natural.image - RMS noise: 7.5 uJy/b
#sim/briggs.image - Beam: 2.0" x 1.4"
#sim/briggs.image - RMS noise: 8.3 uJy/b
#sim/uniform.image - Beam: 1.3" x 1.0"
#sim/uniform.image - RMS noise: 19.9 uJy/b

# BLOCKS of 6 hrs
#for hr in np.linspace(-3,3,6):
#sim/natural.image - Beam: 2.0" x 1.9"
#sim/natural.image - RMS noise: 7.1 uJy/b
#sim/briggs.image - Beam: 1.5" x 1.4"
#sim/briggs.image - RMS noise: 8.3 uJy/b
#sim/briggs0.image - Beam: 1.2" x 1.2"
#sim/briggs0.image - RMS noise: 10.5 uJy/b
#sim/uniform.image - Beam: 1.0" x 1.0"
#sim/uniform.image - RMS noise: 19.2 uJy/b

#for hr in np.linspace(-6,0,6):
#sim/natural.image - Beam: 2.3" x 1.9"
#sim/natural.image - RMS noise: 7.1 uJy/b
#sim/briggs.image - Beam: 1.8" x 1.4"
#sim/briggs.image - RMS noise: 8.4 uJy/b
#sim/briggs0.image - Beam: 1.3" x 1.2"
#sim/briggs0.image - RMS noise: 10.6 uJy/b
#sim/uniform.image - Beam: 1.1" x 1.0"
#sim/uniform.image - RMS noise: 19.6 uJy/b

#for hr in np.linspace(-6,6,6):
#sim/natural.image - Beam: 2.4" x 2.0"
#sim/natural.image - RMS noise: 7.1 uJy/b
#sim/briggs.image - Beam: 1.8" x 1.5"
#sim/briggs.image - RMS noise: 8.4 uJy/b
#sim/uniform.image - Beam: 1.1" x 1.0"
#sim/uniform.image - RMS noise: 20.1 uJy/b
