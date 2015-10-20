#!/usr/bin/env python

# create a mask using bdsm of an image

def make_mask(image_name, threshpix=5, threshisl=3, atrous_do=False, mask_name=None, rmsbox=(55,12), mask_combine=None):

    import sys, os
    import numpy as np
    import pyfits, pyrap
    import lofar.bdsm as bdsm

    # wavelets are required to fit gaussians
    if atrous_do: stop_at = None
    else: stop_at = 'isl'

    # DO THE SOURCE DETECTION
    img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox, \
        thresh_pix=int(threshpix), thresh_isl=int(threshisl), atrous_do=atrous_do, atrous_jmax=3, \
        adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(80,20), \
        stop_at=stop_at, blank_limit=1e-5)

    # WRITE THE MASK FITS
    if mask_name == None: mask_name = image_name+'.newmask'
    if os.path.exists(mask_name): os.system('rm -r ' + mask_name)
    print 'Making mask:', mask_name
    img.export_image(img_type='island_mask', img_format='casa', outfile=mask_name)
    del img

    # do an pixel-by-pixel "OR" operation with a given mask
    if mask_combine != None:
        print "Combining with "+mask_combine
        img = pyrap.images.image(mask_name)
        pixels_mask = img.getdata()
        imgcomb = pyrap.images.image(mask_combine)
        assert imgcomb.shape() == img.shape()
        pixels_mask[np.where(imgcomb.getdata() == 1.)] = 1.
        img.putdata(pixels_mask)
        img.unlock()
        imgcomb.unlock()
        del img
        del imgcomb

    return mask_name

if __name__=='__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v|-V] imagename \n Francesco de Gasperin', version='1.0')
    opt.add_option('-p', '--threshpix', help='Threshold pixel (default=5)', type='int', default=5)
    opt.add_option('-i', '--threshisl', help='Threshold island (default=3)', type='int', default=3)
    opt.add_option('-t', '--atrous_do', help='BDSM extended source detection (default=False)', action='store_true', default=False)
    opt.add_option('-m', '--newmask', help='Mask name (default=imagename with mask in place of image)', default=None)
    opt.add_option('-c', '--combinemask', help='Mask name of a mask to add to the found one (default=None)', default=None)
    opt.add_option('-r', '--rmsbox', help='rms box size (default=55,12)', default='55,12')
    (options, args) = opt.parse_args()
    
    rmsbox = (int(options.rmsbox.split(',')[0]),int(options.rmsbox.split(',')[1]))
    make_mask(args[0].rstrip('/'), options.threshpix, options.threshisl, options.atrous_do, options.newmask, rmsbox, options.combinemask)
