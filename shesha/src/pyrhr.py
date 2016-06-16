# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 09:54:08 2016

@author: J.Meslem - E.Gendron - D.Gratadour
"""

import numpy as np
import matplotlib.pylab as plt
plt.ion()
#from iterkolmo import *
import iterkolmo
import pylab

"""
DOCUMENT
Generate Koolmogorov Phase Screen;
Input :
    - n : Pixels number

Output : Random phase screen p

Note :
    - Better if n is a power of 2
    - Fried parameter r0 is constant and equal to the value of one pixel in meter (For E-ELT : D=39m, r0 = 0.32m)
"""


def generate_kolmo(n):

    Zx,Zy,Xx,Xy,istencil = iterkolmo.create_stencil(n)        # create the indicial matrix
    pxx = iterkolmo.Cxx(n,Zx,Zy,Xx,Xy)                        # covariance matrixes definition
    pzz = iterkolmo.Czz(n,Zx,Zy,istencil)
    pxz = iterkolmo.Cxz(n,Zx,Zy,Xx,Xy,istencil)
    A,B = iterkolmo.AB(pxz,pxx,pzz,n)                        # Matrix creation for Kolmogorov law application
    p = np.zeros((n,n))                             # declaration of the phase screen

    for i in range(2*n):
        p = iterkolmo.extrude(A,B,p,istencil)               #Kolmogorov phase screen

    while np.max(np.abs(p)) > 100:
        p/=10

    return p


def rebin(a, shape):
    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    return a.reshape(sh).mean(-1).mean(1)


def pyr_analysis(n, mod, N, Dtel, obs, nrebin, l, _pix, Amp, mv, Pangle, p, pup=False, noise=False, disp=True):
    """
    DOCUMENT
    Simulate perfect pyramid wavefront sensor analysis without noise;
    Input :
        - n : Pixels number
        - mod : Modulation amplitude in lamda/D unit
        - N : Modulation points number
        - Dtel : Telescope diameter (m)
        - obs : Central obstruction of the telescope
        - nrebin : Binning factor between Pyramid High resolution image and Low resolution (ex: 16)
        - l : Wavelength (m)
        - pix : Angle seen by one pixel (arcsec)
        - Amp : Field amplitude
        - mv : Magnitude in V Band
        - p : Phase screen

    Output :
        - PSF : PSFs on the pyramid after modulation
        - PUPIM_HR : Pupil image seen by the detector behind the pyramid (High resolution)
        - PUPIMB_LR = Pupil image seen by the detector behind the pyramid (Low resolution)
        - dwx, dwy : Wavefront perturbation along axe x and y respectively mesured by pyramid WFS
        - gradx, grady : Real wavefront perturbation of the phase screen

    Note :
        - n must be the same for Kolmogorov
        - Fried parameter r0 is constant and equal to 0.32m with E-ELT

    SCRIPT EXAMPLE:
        from tools import *
        import astropy.io.fits as pyfits
        pupCOMPASS = pyfits.getdata("puppyrCOMPASS.fits")

        imhrCOMPASS = pyfits.getdata("imageCompassPYrHR4.fits")
        phasehrCOMPASS = pyfits.getdata("PhaseCompassPYrHR4.fits")
        phaselrCOMPASS = pyfits.getdata("imageCompassPYrLR4.fits")

        ncompass = phasehrCOMPASS.shape[0]
        n = 1024 # wao.wfs.get_pyrimghr(0).shape
        nrebin = 16 # wao.config.p_wfs0._nrebin
        pup_sep = 16 # wao.config.p_wfs0.nxsub
        mod = 3. # wao.config.p_wfs0.pyr_ampl
        N = 16 # wao.config.p_wfs0.pyr_npts
        Dtel = 8. # wao.config.p_tel.diam
        obs = 0.12 # wao.config.p_tel.cobs
        l = 500e-9  # wao.config.p_wfs0.Lambda*1e-6
        _pix = 0.003222875064238906 # wao.config.p_wfs0._qpixsize
        Amp = 1.
        mv = 10.
        Pangle = pup_sep * nrebin
        p = np.zeros((n,n))
        pup = p.copy()
        p[(n-ncompass)/2:(n+ncompass)/2, (n-ncompass)/2:(n+ncompass)/2] = phasehrCOMPASS*2*np.pi/l*1e-6
        pup[(n-ncompass)/2:(n+ncompass)/2, (n-ncompass)/2:(n+ncompass)/2] = pupCOMPASS
        PSF, PUPIM_HR, PUPIM_LR, dwx, dwy, IA, IB, IC, ID, ppup, pyr = pyr_analysis(n,mod,N,Dtel,obs, nrebin, l,_pix,Amp, mv, Pangle, p, pup=pup)

    """
    if(disp):
        fig1, axes1 = plt.subplots(2, 2) # fig and 2 x 2 nparray of axes
        fig2, axes2 = plt.subplots(1, 2) # fig and 2 x 2 nparray of axes
        fig3, axes3 = plt.subplots(2, 2) # fig and 2 x 2 nparray of axes

    pix = (np.pi*_pix)/(3600*180) #Conversion en rad
    z=l/(n*pix) #pix (m)
    D=Dtel/z #Diametre miroir (pix)
    D = np.round(D)  # on veut un entier coute que coute
    if (D%2 == 1):
        D+=1
    r = D/2.0
    larg=l/Dtel/pix #Lamda/D en pix

    E = np.zeros((n,n), dtype=np.complex64) #Field initialisation

    # on choisit les tqblequx x et y de cette fqcon, pour que, par
    # convention, on adresse les tableaux selon tab[x,y]
    y = np.tile(np.arange(n) - n/2, (n,1))
    x = y.T
    d2 = x*x + y*y
    if(pup is False):
        pup = (d2 < r*r) & (d2 > (obs*r)**2) #Pupille et obstruction centrale

    ppup = pup[int(n/2-D/2)+1 : int(n/2+D/2)+1 , int(n/2-D/2)+1 : int(n/2+D/2)+1]

    #Pangle = int(3*D/4)+1 #Pyramid angle
    #Pangle = int((D/2)*3)+1 #Pyramid angle
    K = (2*np.pi*Pangle)/n
    pyr = K*(np.abs(x)+np.abs(y))  #Pyramide
    pyr = np.fft.fftshift(pyr)
    theta = (np.arange(N))*2*np.pi/N #Angle de modulation

    a=mod*np.cos(theta)
    b=mod*np.sin(theta)
    PUPIM = 0
    PSF = 0
    PUPYR = 0
    # Coeff multiplicatif des fonctions x et y pour etre en unites lam/D
    magic = larg * (2*np.pi/n)

    x = x * magic
    y = y * magic
    if(disp):
        axes1[0, 0].matshow(p) # Phase
        axes1[0, 0].set_title("Phase")
        axes1[0, 1].matshow(p*pup) # Phase in pupil
        axes1[0, 1].set_title("Phase in HR pupil")
        # PSF on top of Pyramid
        PSF0 = np.fft.fftshift(np.abs( np.fft.fft2(pup*Amp*np.exp(1j*(p))) )**2) #PSF not modulated
        axes1[1, 1].matshow(PSF0, cmap="gist_earth")
        axes1[1, 1].plot((PSF0.shape[0]/2,PSF0.shape[0]/2), (0, PSF0.shape[0]), color="blue")
        axes1[1, 1].plot((0, PSF0.shape[0]), (PSF0.shape[0]/2,PSF0.shape[0]/2), color="blue")
        axes1[1, 1].set_title("PSF on top of Pyramid")

    # Modulation LOOP !!!
    for i in range(N):
        phi=a[i]*x + b[i]*y #Modulation Tip/Tilt
        E = pup*Amp*np.exp(1j*(p+phi))
        EPSF = np.fft.fft2(E)
        PSF += np.abs( EPSF )**2 #PSF
        PUPYR = np.abs(np.fft.ifft2((EPSF*np.exp(1j*pyr))))**2
        PUPIM += PUPYR
        if(disp):
            axes1[1, 0].scatter(a[i], b[i], color="red", marker='x', s=20) # Modulation points
            diffractionCircle = pylab.Circle((a[i],b[i]), radius=0.5, alpha=0.5)
            axes1[1, 0].add_patch(diffractionCircle)
            axes1[1, 0].set_title("Modulation points (in lambda/D)")
        print "Modulation iter=", i
        print "min=", np.min(PSF), " max=", np.max(PSF)

    PSF = np.fft.fftshift(PSF)
#    pli(PUPIM,win=1) # pli

    Nph = 10**(-0.4*mv+3) #photons/cm^2/s/Angstrom
    Nph *= 3000*1/200.*(np.pi*(Dtel*100/2.)**2-np.pi*(0.28*Dtel*100/2.)**2) #photons

    S=np.sum(PUPIM)
    PUPIM /= S
    PUPIM *= Nph
    #n = nrebin
    #PUPIM0=np.copy(PUPIM)
    if(noise):
        PUPIM = np.random.poisson(PUPIM) #Bruit de photons
        PUPIM = PUPIM*1.0 + np.random.normal(0,0.5,(n,n)) #Bruit lecture gaussien e- rms/pix

    IA = PUPIM[int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1 , int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1]
    IB = PUPIM[int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1 , int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1]
    IC = PUPIM[int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1 , int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1]
    ID = PUPIM[int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1 , int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1]

    Itot = IA + IB + IC + ID
    Sx = ppup*(IB+ID -(IA+IC))/(Itot)
    Sy = ppup*(IC+ID -(IA+IB))/(Itot)
    dwx = mod*np.sin(0.5*np.pi*Sx)
    dwy = mod*np.sin(0.5*np.pi*Sy) # en unite lam/D
    PUPIM_LR = rebin(PUPIM, (PUPIM.shape[0]/nrebin, PUPIM.shape[1]/nrebin))

    dwx_lr = rebin(dwx, (dwx.shape[0]/nrebin, dwx.shape[1]/nrebin))
    dwy_lr = rebin(dwy, (dwy.shape[0]/nrebin, dwy.shape[1]/nrebin))

    if(disp):
        axes2[0].matshow(PUPIM, cmap="gist_earth")
        axes2[0].set_title("PYRAMID image High Resolution")
        axes2[1].matshow(PUPIM_LR, cmap="gist_earth")
        axes2[1].set_title("PYRAMID image Low Resolution")

        axes3[0, 0].matshow(dwx, cmap="gist_earth")
        axes3[0, 1].matshow(dwy, cmap="gist_earth")
        axes3[0, 0].set_title("GRAD (x) High Resolution")
        axes3[0, 1].set_title("GRAD (y) High Resolution")
        axes3[1, 0].matshow(dwx_lr, cmap="gist_earth")
        axes3[1, 1].matshow(dwy_lr, cmap="gist_earth")
        axes3[1, 0].set_title("GRAD (x) Low Resolution")
        axes3[1, 1].set_title("GRAD (y) Low Resolution")

        fig1.tight_layout(pad=0.1, h_pad=0.5, w_pad=0.5,rect=None)
        fig2.tight_layout(pad=0.1, h_pad=0.5, w_pad=0.5,rect=None)
        fig3.tight_layout(pad=0.1, h_pad=0.5, w_pad=0.5,rect=None)

    """
    gradx=np.zeros((n,n)) #Gradient theorique du front d'onde
    grady=np.zeros((n,n))

    g = np.gradient(p)
    gradx=g[0][int(n/2-D/2)+1 : int(n/2+D/2)+1 , int(n/2-D/2)+1 : int(n/2+D/2)+1]*ppup
    gradx *= l/2/np.pi  #  metres
    gradx /= z  # radians d'angle  (note: z=Dtel/D)
    gradx /= l/Dtel   #   unites lam/D

    grady=g[1][int(n/2-D/2)+1 : int(n/2+D/2)+1 , int(n/2-D/2)+1 : int(n/2+D/2)+1]*ppup
    grady *= l/2/np.pi  #  metres
    grady /= z  # radians d'angle  (note: z=Dtel/D)
    grady /= l/Dtel   #   unites lam/D
    """
    return PSF, PUPIM , PUPIM_LR, dwx, dwy, IA, IB, IC, ID, ppup, pyr


def compPyrCOMPASS(wao):
        #imhrCOMPASS = wao.wfs.get_hrimg_pyr(0)
        phasehrCOMPASS = wao.wfs.get_phase(0)
        pupCOMPASS = wao.config.p_geom.get_mpupil()
        #phaselrCOMPASS = wao.wfs.get_pyrimg(0)
        ncompass = phasehrCOMPASS.shape[0]
        n = wao.wfs.get_pyrimghr(0).shape[0]
        nrebin = wao.config.p_wfs0._nrebin
        pup_sep = wao.config.p_wfs0.nxsub
        mod = wao.config.p_wfs0.pyr_ampl
        N = wao.config.p_wfs0.pyr_npts
        Dtel = wao.config.p_tel.diam
        obs = wao.config.p_tel.cobs
        l = 500e-9  # wao.config.p_wfs0.Lambda*1e-6
        _pix = wao.config.p_wfs0._qpixsize
        Amp = 1.
        mv = 10.
        Pangle = pup_sep * nrebin
        p = np.zeros((n,n))
        pup = p.copy()
        p[(n-ncompass)/2:(n+ncompass)/2, (n-ncompass)/2:(n+ncompass)/2] = phasehrCOMPASS*2*np.pi/l*1e-6
        pup[(n-ncompass)/2:(n+ncompass)/2, (n-ncompass)/2:(n+ncompass)/2] = pupCOMPASS

        PSF, PUPIM_HR, PUPIM_LR, dwx, dwy, IA, IB, IC, ID, ppup, pyr = pyr_analysis(n,mod,N,Dtel,obs, nrebin, l,_pix,Amp, mv, Pangle, p, pup=pup)

"""
from iterkolmo import *
n=1024
A,B,istencil = AB(n)
phase = np.zeros((n,n))

for i in range(2*n):
    phase=extrude(phase,2,A,B,istencil)
pangle = 16*16
mod = 4.
N = 32
Dtel = 8.
l = 600e-9
pix = 0.01/4.
Amp = 1.
mv = 10.
p = 0.
p = phase
PSF, PUPIM, PUPIM_LR, dwx, dwy, IA, IB, IC, ID, ppup, pyr = pyr_analysis(n,20,N,Dtel,l,pix,Amp, mv, pangle, p)
"""
