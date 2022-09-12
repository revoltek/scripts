import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox
import matplotlib.hatch
from matplotlib.patches import Polygon
from matplotlib.colorbar import ColorbarBase
from astropy import units as u
import matplotlib.image as mpimg
import PIL



class ArrowHatch(matplotlib.hatch.Shapes):
    """
    Arrow hatch. Use with hatch="arr{angle}{size}{density}"
                 angle: integer number between 0 and 360
                 size: some integer between 2 and 20
                 density: some integer >= 1
                 https://stackoverflow.com/questions/48334315/how-to-fill-the-bars-of-a-pyplot-barchart-with-arrows
    """
    filled = True
    size = 1

    def __init__(self, hatch, density):
        v1 = [[.355,0], [.098, .1], [.151,.018], [-.355,.018]]
        v2 = np.copy(v1)[::-1]
        v2[:,1] *= -1
        v = np.concatenate((v1,v2))
        self.path = Polygon(v, closed=True, fill=False).get_path()
        self.num_lines = 0
        if len(hatch) >= 5:
            if hatch[:3] == "arr":
                h = hatch[3:].strip("{}").split("}{")
                angle = np.deg2rad(float(h[0]))
                self.size = float(h[1])/10.
                d = int(h[2])
                self.num_rows = 2*(int(density)//6*d)
                self.num_vertices = (self.num_lines + 1) * 2

                R = np.array([[np.cos(angle), -np.sin(angle)],
                              [np.sin(angle), np.cos(angle)]])
                self.shape_vertices = np.dot(R,self.path.vertices.T).T
                self.shape_codes = self.path.codes
        matplotlib.hatch.Shapes.__init__(self, hatch, density)

def setSize(ax, wcs, ra, dec, size_ra, size_dec):
    """
    Properly set bottom left and top right pixel assuming a center and a size in deg
    """
    # bottom
    dec_b = dec - size_dec/2.
    # top
    dec_t = dec + size_dec/2.
    # bottom left
    ra_l = ra-size_ra/np.cos(dec_b*np.pi/180)/2.
    # top right
    ra_r = ra+size_ra/np.cos(dec_t*np.pi/180)/2.

    x,y = wcs.wcs_world2pix([ra_l,ra_r]*u.deg, [dec_b,dec_t]*u.deg, 1, ra_dec_order=True)
    ax.set_xlim(x[1], x[0])
    ax.set_ylim(y[0], y[1])
    return x.astype(int), y.astype(int)


def addScalebar(ax, wcs, z, kpc,  fontsize, color='black'):
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    print("-- Redshift: %f" % z)
    degperpixel = np.abs(wcs.all_pix2world(0,0,0)[1] - wcs.all_pix2world(0,1,0)[1]) # delta deg for 1 pixel
    degperkpc = cosmo.arcsec_per_kpc_proper(z).value/3600.
    pixelperkpc = degperkpc/degperpixel
    fontprops = fm.FontProperties(size=fontsize)
    scalebar = AnchoredSizeBar(ax.transData, kpc*pixelperkpc, '%i kpc' % kpc, 'lower right', fontproperties=fontprops, pad=0.5, color=color, frameon=False, sep=5, label_top=True, size_vertical=2)
    ax.add_artist(scalebar)


def addBeam(ax, hdr, edgecolor='black'):
    """
    hdr: fits header of the file
    """
    from radio_beam import Beam

    bmaj = hdr['BMAJ']
    bmin = hdr['BMIN']
    bpa = hdr['BPA']
    beam = Beam(bmaj*u.deg,bmin*u.deg,bpa*u.deg)

    assert np.abs(hdr['CDELT1']) == np.abs(hdr['CDELT2'])
    pixscale = np.abs(hdr['CDELT1'])
    posx = ax.get_xlim()[0]+bmaj/pixscale
    posy = ax.get_ylim()[0]+bmaj/pixscale
    r = beam.ellipse_to_plot(posx, posy, pixscale *u.deg)
    r.set_edgecolor(edgecolor)
    r.set_facecolor('white')
    ax.add_patch(r)

def addRegion(regionfile, ax, header, alpha=1.0):
    import pyregion
    reg = pyregion.open(regionfile)
    reg = reg.as_imagecoord(header)
    patch_list, artist_list = reg.get_mpl_patches_texts()
    [p.set_alpha(alpha) for p in patch_list]
    for p in patch_list:
        ax.add_patch(p)
    for a in artist_list:
        ax.add_artist(a)

def addCbar(fig, plottype, im, header, int_max, fontsize, cbanchor=[0.127, 0.89, 0.772, 0.02]):
    cbaxes = fig.add_axes(cbanchor)
    cbar = ColorbarBase(cbaxes, cmap=im.get_cmap(), norm=im.norm, orientation='horizontal',  alpha=1.0) #
    # cbar = fig.colorbar(im, cax=cbaxes, orientation='horizontal', pad=0.35, alpha=1.0)
    if plottype == 'stokes':
        cbaxes.xaxis.set_label_text(r'Flux density (mJy beam$^{-1}$)', fontsize=fontsize)
        log_start, log_stop = -2, np.floor(np.log10(int_max)).astype(int)
        cbar.set_ticks([0.001,0.005,0.01,0.05,0.1,0.5,1.,5.,10.])
        # cbar.set_ticks(10**(np.linspace(log_start, log_stop+1, log_stop - log_start +2 )))  # horizontal colorbar
    elif plottype in ['si', 'si+err']:
        cbaxes.xaxis.set_label_text(r'$\alpha_{'+f'{(header["FREQLO"]*1e-6):.0f}'+'MHz}^{'f'{(header["FREQHI"]*1e-6):.0f}'+'MHz}$', fontsize=fontsize)
    elif plottype == 'sierr':
        cbaxes.xaxis.set_label_text(r'$\Delta\alpha_{'+f'{(header["FREQLO"]*1e-6):.0f}'+'MHz}^{'f'{(header["FREQHI"]*1e-6):.0f}'+'MHz}$', fontsize=fontsize)
    else:
        cbaxes.xaxis.set_label_text('Curvature', fontsize=fontsize)
    cbaxes.xaxis.tick_top()
    cbaxes.xaxis.set_label_position('top')
    cbar.ax.tick_params(labelsize=fontsize-3)

def get_sdss_cutout(ra, dec, size=[6, 6], scale=0.5):
    """
    Get an SDSS cutout
    Parameters
    ----------
    ra: RA in deg
    dec: DEC in deg
    size: [int, int], size in arcmin. default = 6
    scale: scale, arcsec per pixel. default = 0.5

    Returns
    -------

    """
    outfile = 'test'
    print(f"wget 'http://skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg?ra={ra}&dec={dec}&scale={scale}&height={size[0]}&width={size[1]}&opt=' -O {outfile}.jpeg")
    os.system(f"wget 'http://skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg?ra={ra}&dec={dec}&scale={scale}&height={size[0]*60}&width={size[1]*60}&opt=' -O {outfile}.jpeg")
    # jpeg_image = PIL.Image.open(f"{outfile}.jpeg")
    # data = np.array(jpeg_image.getdata())
    data = mpimg.imread(f"{outfile}.jpeg")
    print(data)
    return data
