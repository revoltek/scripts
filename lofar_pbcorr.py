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

# Use awimager to create an image of the primary beam and then use it to correct another image

import os, sys, random, optparse
import numpy as np
import pyrap.images as pim
import pyrap.tables as pt

import pyfits as pf
import pywcs

def pbcorr(mosaicfits, pbfits, mosaicpbfits='', mosaicpbcutfits='', Pcut = 0.4):
    
    if not os.path.exists(mosaicfits):
        raise Exception, "missing file: {m}".format(m=mosaicfits)
    if not os.path.exists(pbfits):
        raise Exception, "missing file: {m}".format(m=pbfits)

    if mosaicpbfits=='': mosaicpbfits=mosaicfits+'_PBCOR'
    if mosaicpbcutfits=='': mosaicpbcutfits=mosaicfits+'_PBCUT'
    
    if os.path.isfile(mosaicpbfits):
        print "warning: overwriting {m}".format(m=mosaicpbfits)
        os.system("rm -rf %s" %(mosaicpbfits))
    if os.path.isfile(mosaicpbcutfits):
        print "warning: overwriting {m}".format(m=mosaicpbcutfits)
        os.system("rm -rf %s" %(mosaicpbcutfits))
        
    # work on image
    mosaichead = pf.getheader(mosaicfits)
    mosaicdat = pf.getdata(mosaicfits)
    hdulistout = pf.open(mosaicfits)
    wcsout = pywcs.WCS(mosaichead)
    ra0out = mosaichead.get('OBSRA')
    dec0out = mosaichead.get('OBSDEC')
    
    S = mosaicdat.shape
    if len(S) == 4:
#        print "dim 4"
        nc, nf, ny, nx = S
        dim = 4
    elif len(S) == 2:
#        print "dim 2"
        ny, nx = S
        dim = 2
    else:
        raise Exception, "I don't know how to handle an image with this shape: "+str(S)
    
    mosaiccut = mosaicdat.copy()
    pbcorout = np.zeros_like(mosaicdat)
    
    # work on beam fits
    pbhead = pf.getheader(pbfits)
    pbcordat = pf.getdata(pbfits)
    pbhdulist = pf.open(pbfits)
    pbwcs = pywcs.WCS(pbhead)
    
    pbS = pbcordat.shape
    if len(pbS) == 4:
        print "dim 4"
        pbnc, pbnf, pbny, pbnx = pbS
        ndim = 4
    elif len(pbS) == 2:
        print "dim 2"
        pbny, pbnx = pbS
        ndim = 2
    else:
        raise Exception, "I don't know how to handle an image with this shape: "+str(pbS)
    
    # rescale? 
    if pbnx != nx or pbny != ny:
        #tt[[0,0,1,1,2,2],:][:,[0,0,1,1,2,2]]
        print nx,ny
        print pbnx, pbny
    
        # resample the pb image on the input image grid
        pbcordat_resampled = np.nan*np.ones_like(mosaicdat)
        
        if ndim == 4:
            print "resampling pb image... this may take a while"
            # loop by column 
            prog = 0
            for xi in range(nx):
                if 100*xi/nx > prog:
                    prog = 100*xi/nx
                    if prog%10 == 0:
                        #print prog ,
                        sys.stdout.write(str(prog))
                        sys.stdout.flush()
                    else:
                        #print '.',
                        sys.stdout.write('.')
                        sys.stdout.flush()
                c = np.zeros(ny)
                f = np.zeros(ny)
                x = xi*np.ones(ny)
                y = np.arange(ny)
                pixcrd = np.array([x, y, c, f]).transpose()
                ra, dec, c, f = wcsout.wcs_pix2sky(pixcrd, 0).transpose()
                worldcrd = np.array([ra, dec, c*0, f*0]).transpose()
                pbx, pby, pbc, pbf = pbwcs.wcs_sky2pix(worldcrd, 0).transpose()
                pbx = np.array(pbx, dtype=int)
                pby = np.array(pby, dtype=int)
                
                for i in range(ny):
                    # outside the pb coverage
                    if (pbx[i] < 0) or (pby[i] < 0) or (pbx[i] >= pbnx) or (pby[i] >= pbny):
                        #pbcordat_resampled[0,0,x[i],y[i]] = np.nan
                        pass
                    else:
                        #pbcordat_resampled[0,0,x[i],y[i]] = pbcordat[0,0,pbx[i],pby[i]]
                        #pbcordat_resampled[0,0,x[i],y[i]] = pbcordat[0,0,pbx[i],pby[i]]
                        pbcordat_resampled[0,0,y[i],x[i]] = pbcordat[0,0,pby[i],pbx[i]]
                        # odd pixels are transposed? dont fully understand, seems to work but may be problematic later
                        # not actually - was a bug
                        # still a bug... pby and pbx should be transposed
            
        elif ndim == 2:
            print "resampling pb image... this may take a while"
            # loop by column 
            prog = 0
            for xi in range(nx):
                if 100*xi/nx > prog:
                    prog = 100*xi/nx
                    if prog%10 == 0:
                        #print prog ,
                        sys.stdout.write(str(prog))
                        sys.stdout.flush()
                    else:
                        #print '.',
                        sys.stdout.write('.')
                        sys.stdout.flush()
                x = xi*np.ones(ny)
                y = np.arange(ny)
                pixcrd = np.array([x, y]).transpose()
                ra, dec = wcsout.wcs_pix2sky(pixcrd, 0).transpose()
                worldcrd = np.array([ra, dec]).transpose()
                pbx, pby = pbwcs.wcs_sky2pix(worldcrd, 0).transpose()
                pbx = np.array(pbx, dtype=int)
                pby = np.array(pby, dtype=int)
                
                for i in range(ny):
                    # outside the pb coverage
                    if (pbx[i] < 0) or (pby[i] < 0 ) or (pbx[i] >= pbnx ) or (pby[i] >= pbny):
                        #pbcordat_resampled[0,0,x[i],y[i]] = np.nan
                        pass
                    else:
                        pbcordat_resampled[y[i],x[i]] = pbcordat[pbx[i],pby[i]]
            
#        pbcordat_resampled[pbcordat_resampled**0.5<Pcut] = np.nan  # set nan beyond 
        mosaiccut[pbcordat_resampled**0.5<Pcut] = np.nan  # set nan beyond 
#        mosaiccut[np.isnan(pbcordat_resampled)] = np.nan  # set nan beyond 
        mosaiccor = mosaiccut/(pbcordat_resampled**0.5)   #assuming awimager output is avgpb
        
        # save rescaled pb for inspection
        pb_rescaled_fits=mosaicfits.replace('.fits','')+'.pb.fits'
        if os.path.isfile(pb_rescaled_fits):
            print "warning: overwriting {m}".format(m=pb_rescaled_fits)
            os.system("rm -rf %s" %(pb_rescaled_fits))
        pf.writeto(pb_rescaled_fits, pbcordat_resampled**0.5, header=mosaichead)

    else:
        mosaiccut[pbcordat**0.5<Pcut] = np.nan
        mosaiccor = mosaiccut /(pbcordat**0.5)   #assuming awimager output is avgpb
        
    pf.writeto(mosaicpbcutfits, mosaiccut, header=mosaichead)
    pf.writeto(mosaicpbfits, mosaiccor, header=mosaichead)


def main(opts, args):
    
    if opts.output == None:
        opts.output = args[1].rstrip('/')+'_PBcor'
    if opts.pbcut > 1 or opts.pbcut < 0:
        print "ERROR: PBcut must be between 0 and 1."
        sys.exit(1)

    ms = args[0].rstrip('/')
    img_in = args[1].rstrip('/')
    v = opts.verbose
    fov = opts.fov
    img_out = opts.output.rstrip('/')
    data_col = opts.datacol
    pbcut = opts.pbcut
    npix = opts.npix

    img_tmp = 'awimg_'+str(random.randint(0,1e9))
    cellsize = np.ceil(fov*3600./npix)
    print 'running: awimager ms='+ms+' image='+img_tmp+' niter=0 data='+data_col+' weight=briggs robust=0 npix='+str(npix)+' cellsize='+str(cellsize)+'arcsec padding=1. stokes=I operation=mfclark FindNWplanes=True PBCut=1e-2'
    os.system('awimager ms='+ms+' image='+img_tmp+' niter=0 data='+data_col+' weight=briggs robust=0 npix='+str(npix)+' cellsize='+str(cellsize)+'arcsec padding=1. stokes=I operation=mfclark FindNWplanes=True PBCut=1e-2')
    os.system('avgpbz.py '+img_tmp+'0.avgpb')
    os.system('image2fits in='+img_tmp+'0.avgpbz out='+img_tmp+'.fits')

    pbcorr(mosaicfits = img_in, pbfits = img_tmp+'.fits', Pcut = pbcut)

    os.system('rm -rf '+img_tmp+'*')

######################################
usage = "%prog [options] msfile image"
opt = optparse.OptionParser(usage)
opt.add_option('-v','--verbose',help='Verbose mode', default=False, action='store_true')
opt.add_option('-f','--fov',help='Linear size of the field of view (degrees) [5]', default=5, type='float')
opt.add_option('-o','--output',help='Output image name [input image + "_PBcor"]', default=None, type='string')
opt.add_option('-d','--datacol',help='Data column to image [CORRECTED_DATA]', default='CORRECTED_DATA', type='string')
opt.add_option('-p','--pbcut',help='Cut at PB value, must be a number from 0 to 1 [0.4]', default='0.4', type='float')
opt.add_option('-n','--npix',help='Number of pixel in the pb image made by awimager [512]', default='512', type='int')
opts,args = opt.parse_args()
if len(args) != 2:
    print 'Need required argument msname, image'
else:
    main(opts,args)
