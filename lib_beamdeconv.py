#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
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


import numpy as np

# Imports for the custom kernels - https://github.com/radio-astro-tools/radio_beam/blob/master/radio_beam/beam.py
from astropy.modeling.models import Ellipse2D, Gaussian2D
from astropy.convolution import Kernel2D
from astropy.convolution.kernels import _round_up_to_odd_integer
class EllipticalGaussian2DKernel(Kernel2D):
    """
    2D Elliptical Gaussian filter kernel.
    The Gaussian filter is a filter with great smoothing properties. It is
    isotropic and does not produce artifacts.
    Parameters
    ----------
    stddev_maj : float
        Standard deviation of the Gaussian kernel in direction 1
    stddev_min : float
        Standard deviation of the Gaussian kernel in direction 1
    position_angle : float
        Position angle of the elliptical gaussian
    x_size : odd int, optional
        Size in x direction of the kernel array. Default = support_scaling *
        stddev.
    y_size : odd int, optional
        Size in y direction of the kernel array. Default = support_scaling *
        stddev.
    support_scaling : int
        The amount to scale the stddev to determine the size of the kernel
    mode : str, optional
        One of the following discretization modes:
            * 'center' (default)
                Discretize model by taking the value
                at the center of the bin.
            * 'linear_interp'
                Discretize model by performing a bilinear interpolation
                between the values at the corners of the bin.
            * 'oversample'
                Discretize model by taking the average
                on an oversampled grid.
            * 'integrate'
                Discretize model by integrating the
                model over the bin.
    factor : number, optional
        Factor of oversampling. Default factor = 10.
    See Also
    --------
    Box2DKernel, Tophat2DKernel, MexicanHat2DKernel, Ring2DKernel,
    TrapezoidDisk2DKernel, AiryDisk2DKernel, Gaussian2DKernel,
    EllipticalTophat2DKernel
    Examples
    --------
    Kernel response:
     .. plot::
        :include-source:
        import matplotlib.pyplot as plt
        from radio_beam import EllipticalGaussian2DKernel
        gaussian_2D_kernel = EllipticalGaussian2DKernel(10, 5, np.pi/4)
        plt.imshow(gaussian_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """
    _separable = True
    _is_bool = False

    def __init__(self, stddev_maj, stddev_min, position_angle,
                 support_scaling=8, **kwargs):
        self._model = Gaussian2D(1. / (2 * np.pi * stddev_maj * stddev_min), 0,
                                 0, x_stddev=stddev_maj, y_stddev=stddev_min,
                                 theta=position_angle)

        try:
            from astropy.modeling.utils import ellipse_extent
        except ImportError:
            raise NotImplementedError("EllipticalGaussian2DKernel requires"
                                      " astropy 1.1b1 or greater.")

        max_extent = \
            np.max(ellipse_extent(stddev_maj, stddev_min, position_angle))
        self._default_size = \
            _round_up_to_odd_integer(support_scaling * 2 * max_extent)
        super(EllipticalGaussian2DKernel, self).__init__(**kwargs)
        self._truncation = np.abs(1. - 1 / self._array.sum())


def deconvolve_ell(t_bmaj,t_bmin,t_bpa,bmaj,bmin,bpa):
    """
    wrapper for deconvolve() to work in elliptical coordinates
    """
    a,b,c = elliptic2quadratic(bmaj,bmin,bpa)
    t_a,t_b,t_c = elliptic2quadratic(t_bmaj,t_bmin,t_bpa)
    c_a,c_b,c_c = deconvolve(t_a,t_b,t_c,a,b,c)
    return quadratic2elliptic(c_a,c_b,c_c)


def quadratic2elliptic(A,B,C,D=0,E=0,F=-np.log(2)):
    """Invert:
    (A0 cos^2 phi + c0 sin^2 phi)k = A
    (A0-C0)sin 2phi k = B
    (a0 sin^2 phi + c0 cos^2 phi) k = C
    returns bmaj,bmin,bpa[,xc,y if D,E != 0]"""
    if (np.isinf(A) and np.isinf(B) and np.isinf(C)):#delta function
        return 0., 0., 0.
    assert (B**2 - 4*A*C) != 0, "It is parabolic, not elliptic or hyperbolic"
    if A!=C:#not a circle so should be able to solve second equation
        #A-C = A0 cos^2 phi + c0 sin^2 phi - a0 sin^2 phi - c0 cos^2 phi
        #A-C = A0 (cos^2 phi - sin^2 phi)- c0 (-sin^2 phi  + c0 cos^2 phi) = (A0 - C0) cos 2phi
        #(cos^2 phi - sin^2 phi) = cos 2 phi        
        phi = np.arctan2(B,(A-C))/2.#choose your own modulus
    else:#circle
        phi = np.pi/4.#arbitrary, nonzero cos and sin for the rest
        phi = 0.
    #rotate x,y to phi to x',y' =  xcos - ysin, xsin + ycos
    #then expand A(x' - xc)^2 + B(x'-xc)(y'-yc) + C(y' - yc)^2 + D(x' - xc) + E(y' - yc) + F = 0
    #then now the coord sys is unrotated relative to the quadratic paramters
    c = np.cos(phi)
    s = np.sin(phi)
    c2 = c*c
    s2 = s*s
    A1 = A*c2 + B*c*s + C*s2
    B1 = 2.*(C-A)*s*c+B*(c2-s2)#should be zero since there's no cross term in the unroated
    #print "Rotated cross term: {0} =? 0".format(B1)
    C1 = A*s2 - B*c*s + C*c2
    D1 = D*c + E*s
    E1 = -D*s + E*c
    assert (A1 != 0) and (C1 != 0), "degenerate between ellipse and hyperbola"
    #complete square for x's and y's
    #A1(x-xc)^2 + C1(y-yc)^2 + F = A1xx - 2A1xxc + A1xcxc + C1yy - 2C1yyc + C1ycyc = A1() + D1() + C1() + E1() + A1xcxc + C1ycyc + F =0 
    xc1 = D1/(-2.*A1)
    yc1 = E1/(-2.*C1)#solutions (in rotated frame)
    #now unrotate them
    xc = xc1*c - yc1*s#might be xc1*c + yc1*s
    yc = xc1*s + yc1*c#might be -xc1*s + yc1*c
    #bring the remainder of completed squares to rhs with F
    rhs = -F + D1**2/(4.*A1)  + E1**2/(4.*C1)
    #(x-xc)^2/(bmaj/2)^2 + (y-yc)^2/(bmin/2)^2 = 1
    #A1(x-xc)^2 + C1*y-yc)^2 = rhs
    
    A0 = A1/rhs
    C0 = C1/rhs
    bmaj = np.sign(A0)*2.*np.sqrt(1./np.abs(A0))
    bmin = np.sign(C0)*2.*np.sqrt(1./np.abs(C0))
    assert bmaj*bmin > 0, "Hyperbolic solution ;) inversion success but not physical."
        #return None,None,None
    if bmin > bmaj:
        temp = bmin
        bmin = bmaj
        bmaj = temp
    bpa = phi# - np.pi/2.#starts at y
    if E==0 and D==0:
        return bmaj,bmin,bpa*180./np.pi
    return bmaj,bmin,bpa*180./np.pi,xc,yc
        

def elliptic2quadratic(bmaj,bmin,pa,xc=0,yc=0,k=np.log(2)):
    '''a*x**2 + b*x*y + c*y**2 + d*x + e*y + f = 0
    pa in deg
    return A,B,C[,D,E,F if xc,yc!=0]'''

    #unrotated solution
    a0 = k/(bmaj/2.)**2
    c0 = k/(bmin/2.)**2
    theta = (pa + 90.)*np.pi/180.
    #Rotated Solution
    cos2 = np.cos(theta)**2
    sin2 = np.sin(theta)**2
    A = (a0*cos2 + c0*sin2)
    C = (c0*cos2 + a0*sin2)
    B = (a0 - c0 )*np.sin(2.*theta)
    #Now move center
    D = -2.*A*xc - B*yc
    E = -2.*C*yc - B*xc
    F = A*xc**2 + B*xc*yc + C*yc**2 - 1./k
    if xc==0 and yc==0:
        return A,B,C
    return A,B,C,D,E,F


def deconvolve(A1,B1,C1,A2,B2,C2):
    '''Solves analytically G(A1,B1,C1) = convolution(G(A2,B2,C2), G(Ak,Bk,Ck))
    Returns Ak,Bk,Ck
    A,B,C are quadratic parametrization.
    If you have bmaj,bmin,bpa, then get A,B,C = ecliptic2quadratic(0,0,bmaj,bmin,bpa)
    
    Returns (np.inf,np.inf,np.inf) if solution is delta function'''
    D = B1**2 - 2*B1*B2 + B2**2 - 4*A1*C1 + 4* A2* C1 + 4* A1* C2 - 4* A2* C2
    if (np.abs(D) < 10*(1-2./3.-1./3.)):
        return np.inf,np.inf,np.inf#delta function
    if (D<0.):
        #print "Inverse Gaussian, discriminant D:",D
        #ie. hyperbolic solution, still valid but elliptic representation is impossible instead you get hyperbolic parameters: negative bmaj/bmin
        pass
    Ak = (-A2* B1**2 + A1* B2**2 + 4* A1* A2* C1 - 4* A1* A2* C2)/D
    Bk = (-B1**2 *B2 + B1* B2**2 + 4* A1* B2* C1 - 4* A2* B1* C2)/D
    Ck = (B2**2 *C1 - B1**2 *C2 + 4* A1* C1* C2 - 4* A2* C1* C2)/D
    assert (Bk*Bk - 4*Ak*Ck) != 0, "Indifinite deconvolution det = 0"
    return Ak,Bk,Ck


def convolve(A1,B1,C1,A2,B2,C2):
    '''
        Convolves two gaussians with quadratic parametrization:
        A,B,C are quadratic parametrization.
        If you have bmaj,bmin,bpa, then get A,B,C = elliptic2quadratic(0,0,bmaj,bmin,bpa)
        Where g = factor*Exp(-A*X**2 - B*X*Y - C*Y**2)
    '''
    D1 = 4.*A1*C1 - B1**2
    D2 = 4.*A2*C2 - B2**2
    D3 = -2.*B1 * B2 + 4.*A2*C1 + 4.*A1*C2 + D1+D2
    D4 = C2*D1+C1*D2
    #Non-solvable cases
    if (D1*D2*D3*D4 == 0):
        print ("Can't convolve...")
        return (None,None,None)
    if (D3 < 0):#always imaginary
        print(("D3 < 0, Imaginary solution",D3))
        return (None,None,None)
    factor = 2.*np.pi*np.sqrt(D1 + 0j)*np.sqrt(D2 + 0j)/np.sqrt(D3/D4 + 0j)/np.sqrt(D4/(D1*D2) + 0j)
    if np.abs(np.imag(factor)) > 10.*(7./3 - 4./3 - 1.):
        print ("Imaginary result somehow...")
        return (None,None,None)
    factor = np.real(factor)
    A = (A2*D1 + A1 * D2)/D3
    B = (B2*D1+B1*D2)/D3
    C = D4/D3
    k = np.log(factor*2.)
    return A,B,C#,factor

def findCommonBeam(beams, debugplots=False, confidence=0.005):
    '''Given a list `beams` where each element of beams is a list of elliptic beam parameters (bmaj_i,bmin_i, bpa_i)
    with bpa in degrees
    return the beam parameters of the common beam of minimal area.
    
    Common beam means that all beams can be convolved to the common beam.
    
    `confidence` parameter is basically how confident you want solution. So 0.01 is knowing solution to 1%.
    Specifically it's how long to sample so that there are enough statistics to properly sample likelihood with required accuracy.
    default is 0.005. Computation time scale inversely with it.'''
    def beamArea(bmaj,bmin,bpa=None):
        return bmaj*bmin*np.pi/4./np.log(2.)
    def isCommonBeam(beamCandQuad,beamsQuad):
        for beamQuad in beamsQuad:
            try:
                Ak,Bk,Ck = deconvolve(beamCandQuad[0],beamCandQuad[1],beamCandQuad[2],beamQuad[0],beamQuad[1],beamQuad[2])
                bmaj,bmin,bpa = quadratic2elliptic(Ak,Bk,Ck)
            except:
                return False
        return True 
    def samplePrior(beamLast, beamsQuad):
        iter = 0
        while True:
            std = 1.5
            beam = [beamLast[0]*np.exp(np.log(std)*np.random.uniform(low=-1,high=1.)),
                     beamLast[1]*np.exp(np.log(std)*np.random.uniform(low=-1,high=1.)),
                     beamLast[2] + np.random.uniform(low=-5,high=5)]
            #beam[0] = np.abs(beam[0])
            #beam[1] = np.abs(beam[1])
            if beam[1] > beam[0]:
                temp = beam[1]
                beam[1] = beam[0]
                beam[0] = temp
            while beam[2] > 90.:
                beam[2] -= 180.
            while beam[2] < -90.:
                beam[2] += 180.
            A,B,C = elliptic2quadratic(*beam)
            if isCommonBeam((A,B,C),beamsQuad):
                return beam
            iter += 1
        
    # TODO: is it OK to minimize only area? Or we should put a penalty on maj/min to avoid very long beams?
    def misfit(beam,areaLargest):
        area = beamArea(*beam)
        L2 = (area - areaLargest)**2/2.
        return L2
    #Get beam areas
    N = len(beams)
    areas = []
    beamsQuad = []
    i = 0
    while i < N:
        areas.append(beamArea(*beams[i]))
        beamsQuad.append(elliptic2quadratic(*beams[i]))
        i += 1
    beam0 = beams[np.argmax(areas)]
    areaLargest = np.max(areas)
    beam0Quad = elliptic2quadratic(*beam0)
    if isCommonBeam(beam0Quad,beamsQuad):
        return beam0
    else:
        bmajMax = np.max(beams,axis=0)[0]
        beam0 = [bmajMax,bmajMax,0.]
    #MC search, 1/binning = confidence
    binning = int(1./confidence)
    Nmax = 1e6
    beamsMH = np.zeros([binning*binning,3],dtype=np.double)
    beamsMul = np.zeros(binning*binning,dtype=np.double)
    beamsMH[0,:] = beam0
    beamsMul[0] = 1
    accepted = 1
    Si = misfit(beam0,areaLargest)
    Li = np.exp(-Si)
    maxL = Li
    maxLBeam = beam0
    iter = 0
    while accepted < binning**2 and iter < Nmax:
        beam_j = samplePrior(beamsMH[accepted-1], beamsQuad)
        Sj = misfit(beam_j,areaLargest)
        Lj = np.exp(-Sj)
        #print("Sj = {}".format(Sj))
        if Sj < Si or np.log(np.random.uniform()) < Si - Sj:
            Si = Sj
            beamsMH[accepted,:] = beam_j
            beamsMul[accepted] += 1
            #print("Accepted")
            accepted += 1
        else:
            beamsMul[accepted-1] += 1
        if Lj > maxL:
            maxL = Lj
            maxLBeam = beam_j
        iter += 1
    if accepted == binning**2:
        pass
        #print("Converged in {} steps with an acceptance rate of {}".format(iter,float(accepted)/iter))
    else:
        beamsMH = beamsMH[:iter,:]
        beamsMul = beamsMul[:iter]
    if debugplots:
        import pylab as plt
#         plt.hist(beamsMH[:,0],bins=binning)
#         plt.show()
#         plt.hist(beamsMH[:,1],bins=binning)
#         plt.show()
#         plt.hist(beamsMH[:,2],bins=binning)
#         plt.show()
        from matplotlib.patches import Ellipse
        ax = plt.subplot(1,1,1)
        ax.add_artist(Ellipse(xy=(0,0), width=maxLBeam[0], height=maxLBeam[1], angle=maxLBeam[2], facecolor="none",edgecolor='red',alpha=1,label='common beam'))
        for beam in beams:
            ax.add_artist(Ellipse(xy=(0,0), width=beam[0], height=beam[1], angle=beam[2], facecolor="none",edgecolor='black',ls='--',alpha=1))
        ax.set_xlim(-0.5,0.5)
        ax.set_ylim(-0.5,0.5)
        plt.legend(frameon=False)
        plt.show()
        
#     meanBeam = np.sum(beamsMH.T*beamsMul,axis=1)/np.sum(beamsMul)
#     stdBeam = np.sqrt(np.sum(beamsMH.T**2*beamsMul,axis=1)/np.sum(beamsMul) - meanBeam**2)
#     print ("(Gaussian) beam is {} +- {}".format(meanBeam,stdBeam))
#     logmeanBmajBmin = np.sum(np.log(beamsMH[:,:2]).T*beamsMul,axis=1)/np.sum(beamsMul)
#     logstdBmajBmin = np.sqrt(np.sum(np.log(beamsMH[:,:2]).T**2*beamsMul,axis=1)/np.sum(beamsMul) - logmeanBmajBmin**2)
#     logstdBmajBminu = np.exp(logmeanBmajBmin + logstdBmajBmin) - np.exp(logmeanBmajBmin)
#     logstdBmajBminl = np.exp(logmeanBmajBmin) - np.exp(logmeanBmajBmin - logstdBmajBmin)
#     logmeanBmajBmin = np.exp(logmeanBmajBmin)
#     print("(Lognormal) bmaj/bmin is {} + {} - {}".format(logmeanBmajBmin,logstdBmajBminu,logstdBmajBminl))
#     print ("Max Likelihood beam is {}".format(maxLBeam))
    return maxLBeam
    
def fftGaussian(A,B,C,X,Y):
    D = 4*A*C-B**2
    return 2*np.pi/np.sqrt(D)*np.exp(-4*np.pi/D*(-C*X**2 +B*X*Y -A*Y**2))

def gaussian(A,B,C,X,Y):
    return np.exp(-A*X**2 - B*X*Y - C*Y**2)

def psfTGSS1(dec):
    '''input declination in degrees
    return bmaj(arcsec), bmin(arcsec), bpa(degrees)'''
    if dec > 19.0836824:
        return 25.,25.,0.
    else:
        return 25.,25./np.cos(np.pi*(dec-19.0836824)/180.),0.

def test_elliptic2quadratic():
    
    for i in range(100):
        bpa = np.random.uniform()*180.-90.#deg
        bmaj = np.random.uniform()
        bmin = np.random.uniform()*bmaj
    
        A,B,C = elliptic2quadratic(bmaj,bmin,bpa)   
        bmaj2,bmin2,bpa2 = quadratic2elliptic(A,B,C)
        assert np.isclose(bmaj,bmaj2) and np.isclose(bmin,bmin2) and np.isclose(bpa,bpa2), "Failed to pass {},{},{} != {},{},{}".format(bmaj,bmin,bpa,bmaj2,bmin2,bpa2)
    return True

def test_convolvedeconvolve(N=100):
    for i in range(N):
        bpa = np.random.uniform()*180.-90.#deg
        bmaj = np.random.uniform()
        bmin = np.random.uniform()*bmaj
    
        A1,B1,C1 = elliptic2quadratic(bmaj,bmin,bpa) 
        
        bpa2 = np.random.uniform()*180.-90.#deg
        bmaj2 = np.random.uniform()
        bmin2 = np.random.uniform()*bmaj2
    
        A2,B2,C2 = elliptic2quadratic(bmaj2,bmin2,bpa2)  
        
        Ac,Bc,Cc = convolve(A1,B1,C1,A2,B2,C2)
        
        Ak,Bk,Ck = deconvolve(Ac,Bc,Cc,A1,B1,C1)
        
        bmaj2_,bmin2_,bpa2_ = quadratic2elliptic(Ak,Bk,Ck)
        
        assert np.isclose(bmaj2_,bmaj2) and np.isclose(bmin2_,bmin2) and np.isclose(bpa2_,bpa2), "Failed to pass {},{},{} != {},{},{}".format(bmaj2_,bmin2_,bpa2_,bmaj2,bmin2,bpa2)
    return True

def test_deltaFunctionDeconvolve():
    bpa = np.random.uniform()*180.-90.#deg
    bmaj = np.random.uniform()
    bmin = np.random.uniform()*bmaj

    A1,B1,C1 = elliptic2quadratic(bmaj,bmin,bpa) 
    #deconv same beam
    Ak,Bk,Ck = deconvolve(A1,B1,C1,A1,B1,C1)
    bmaj_d, bmin_d, bpa_d = quadratic2elliptic(Ak,Bk,Ck)
    assert bmaj_d==0 and bmin_d==0 and bpa_d==0,"Supposed to be the delta"
    return True

def test_timing():
    from time import clock
    i = 0
    t1 = clock()
    for i in range(10000):
        bpa = np.random.uniform()*180.-90.#deg
        bmaj = np.random.uniform()
        bmin = np.random.uniform()*bmaj
    
        A1,B1,C1 = elliptic2quadratic(bmaj,bmin,bpa) 
        
        bpa2 = np.random.uniform()*180.-90.#deg
        bmaj2 = np.random.uniform()
        bmin2 = np.random.uniform()*bmaj2
    
        A2,B2,C2 = elliptic2quadratic(bmaj2,bmin2,bpa2)  
        
        Ac,Bc,Cc = convolve(A1,B1,C1,A2,B2,C2)
        
        Ak,Bk,Ck = deconvolve(Ac,Bc,Cc,A1,B1,C1)
        
        bmaj2_,bmin2_,bpa2_ = quadratic2elliptic(Ak,Bk,Ck)
    print(("Time avg. ~ {} seconds".format((clock()-t1)/10000)))
        
def test_findCommonBeam():
    np.random.seed(1234)
    for i in range(10):
        beams = []
        for i in range(3):
            bpa = np.random.uniform()*180.-90.#deg
            bmaj = np.random.uniform()
            bmin = np.random.uniform()*bmaj
            beams.append((bmaj,bmin,bpa))
        commonBeam = findCommonBeam(beams,debugplots=True)
        print(("Common beam amongst {} is {}".format(beams,commonBeam)))
