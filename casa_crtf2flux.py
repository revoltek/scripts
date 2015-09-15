#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
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

# Usage: casapy  --nogui --nologger -c ~/bin/scripts/casa_crtf2flux.py regions.crtf noise.crtf image1 image2 image3 ...
# select flux in image1,image2... above a certain sigma threshold
# images must be casa_img_armonize-ed
# for every region in regions.crtf extract the flux
# if more then 1 img: do integrated spidx of that region

# to set reffreq:
# imhead(imagename,'put','restfreq',{'value': 323098000.0, 'unit': 'Hz'})

# to solve for LELAttribute: coordinates of operands mismatch
# reorder axis in casa with imtrans (typically order='0132' to switch freq and stokes)

# to add from scratch missing axes:
# ia.open(filename)
# im2 = ia.adddegaxes(stokes='I') # if stokes axis missing
# im3 = im2.adddegaxes(spectral=T) # if spectral axis missing
# im3.subimage('test.image')
# then maybe reorder as explained above

# number of sigma below which the image is masked out
# if 0 then the entire region is used 
nsigma = 0

import os, sys, glob
import numpy as np
import linearfit
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
images = sys.argv[7:]
regionfile = sys.argv[5]
noiseregionfile = sys.argv[6]

# get the right freq not the rest freq
def getFreq(image):
    h = imhead(imagename=image, mode='list')
    for i in xrange(10):
        if h['ctype'+str(i+1)] == 'Frequency':
            return float(h['crval'+str(i+1)])

print "Split region file"
lines = sum(1 for line in open(regionfile))
with open(regionfile) as f:
    for i, l in enumerate(f):
        if i == 0: head = l
        else:
            with open('__reg'+str(i)+'.crtf', 'w+') as fo:
                fo.write(head)
                fo.write(l)

print "Calculate RMSs:"
rms = []
freqs = []
lelexpr = []
for i, image in enumerate(images):
    rms.append(imstat(imagename=image, region=noiseregionfile)['rms'][0])
    freqs.append(getFreq(image))
    print image+ "(freq:", freqs[i], ") - rms:", rms[i], 'Jy/b'
    # make mask that selects only pixels above the noise in all or one of the images
    lelexpr.append('"'+image+'" > '+str(nsigma)+'*'+str(rms[i]))

print "Calculate fluxes:"
flux = {}
err = {}
for region in sorted(glob.glob('__reg*crtf')):
    flux[region] = []
    err[region] = []
    print "### Working region:", region
    for i, image in enumerate(images):

        # get beam in pixel
        bmaj = qa.convert(imhead(imagename=image,mode='get',hdkey='beammajor'),'deg')
        bmin = qa.convert(imhead(imagename=image,mode='get',hdkey='beamminor'),'deg')
        ia.open(image)
        cs = ia.coordsys()
        xpix = cs.increment()['numeric'][0]*180/np.pi
        ypix = cs.increment()['numeric'][1]*180/np.pi
        cs.done()
        ia.close()
        pixperbeam = (bmaj['value']*bmin['value']*1.1331)/abs(xpix*ypix)

        # in some cases casa crashes if the first image of the lel expression is not the one in imagename
        lelexpr[0], lelexpr[i] = lelexpr[i], lelexpr[0]

        stat = imstat(imagename=image, region=region, mask='&&'.join(lelexpr)) # cut at 3 sigma in all images
        print "Use this cut:"+str(lelexpr)
        #stat = imstat(imagename=image, region=region, mask=lelexpr[0]) # cut at 3 sigma in THIS image
        #print "Cut at 3 sigma on this image"
        flux[region].append(stat['flux'][0])
        err[region].append(rms[i]*np.sqrt(stat['npts'][0]/pixperbeam)) 
        print image+ " - flux:", flux[region][i], '±', err[region][i], 'Jy'
        lelexpr[i], lelexpr[0] = lelexpr[0], lelexpr[i]

if len(images) > 1:
    print "Calculate spidx"
    for region in sorted(glob.glob('__reg*crtf')):
        (a, b, sa, sb) = linearfit.linear_fit(x=np.log10(freqs), y=np.log10(flux[region]), yerr=0.434*np.array(err[region])/np.array(flux[region]))
        print region+" - spidx:", a, '±', sa
        fig = plt.figure(figsize=(8, 8))
        fig.subplots_adjust(wspace=0)
        ax = fig.add_subplot(110)
        ax.set_xlabel(r'Log Frequency [Hz]')
        ax.set_ylabel(r'Log Flux [Jy]')
        ax.errorbar(np.log10(freqs), np.log10(flux[region]), yerr=0.434*np.array(err[region])/np.array(flux[region]), fmt='ko')
        ax.plot(np.log10(freqs), [b+a*x for x in np.log10(freqs)], 'r-', label=r'$\alpha$={:.2f}$\pm${:.2f}'.format(a,sa))
        ax.legend(loc=1)
        fig.savefig(region+'.pdf', bbox_inches='tight')
        fig.clf()

os.system('rm __reg*crtf')
print "Done."
