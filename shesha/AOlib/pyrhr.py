# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 09:54:08 2016

@author: J.Meslem - E.Gendron - D.Gratadour - FV
"""

import numpy as np
import matplotlib.pylab as plt
plt.ion()
#from iterkolmo import *
import shesha.util.iterkolmo
import pylab
import sys
import os
sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/src/")
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams['image.cmap'] = 'viridis'
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

    Zx, Zy, Xx, Xy, istencil = iterkolmo.create_stencil(n)  # create the indicial matrix
    pxx = iterkolmo.Cxx(n, Zx, Zy, Xx, Xy)  # covariance matrixes definition
    pzz = iterkolmo.Czz(n, Zx, Zy, istencil)
    pxz = iterkolmo.Cxz(n, Zx, Zy, Xx, Xy, istencil)
    A, B = iterkolmo.AB(pxz, pxx, pzz,
                        n)  # Matrix creation for Kolmogorov law application
    p = np.zeros((n, n))  # declaration of the phase screen

    for i in range(2 * n):
        p = iterkolmo.extrude(A, B, p, istencil)  #Kolmogorov phase screen

    while np.max(np.abs(p)) > 100:
        p /= 10

    return p


def rebin(a, shape):
    sh = int(shape[0]), int(a.shape[0]) // int(shape[0]), int(
            shape[1]), int(a.shape[1]) // int(shape[1])
    return a.reshape(sh).mean(-1).mean(1)


def pyr_analysis(n, mod, N, Dtel, obs, nrebin, l, _pix, Amp, mv, Pangle, p, pupsub,
                 pup=False, noise=False, disp=True):
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
        n = 1024 # wao.supervisor.wfs.get_pyrimghr(0).shape
        nrebin = 16 # wao.supervisor.config.p_wfs0._nrebin
        pup_sep = 16 # wao.supervisor.config.p_wfs0.nxsub
        mod = 3. # wao.supervisor.config.p_wfs0.pyr_ampl
        N = 16 # wao.supervisor.config.p_wfs0.pyr_npts
        Dtel = 8. # wao.supervisor.config.p_tel.diam
        obs = 0.12 # wao.supervisor.config.p_tel.cobs
        l = 500e-9  # wao.supervisor.config.p_wfs0.Lambda*1e-6
        _pix = 0.003222875064238906 # wao.supervisor.config.p_wfs0._qpixsize
        Amp = 1.
        mv = 10.
        Pangle = pup_sep * nrebin
        p = np.zeros((n,n))
        pup = p.copy()
        p[(n-ncompass)/2:(n+ncompass)/2, (n-ncompass)/2:(n+ncompass)/2] = phasehrCOMPASS*2*np.pi/l*1e-6
        pup[(n-ncompass)/2:(n+ncompass)/2, (n-ncompass)/2:(n+ncompass)/2] = pupCOMPASS
        PSF, PUPIM_HR, PUPIM_LR, dwx, dwy, IA, IB, IC, ID, ppup, pyr = pyr_analysis(n,mod,N,Dtel,obs, nrebin, l,_pix,Amp, mv, Pangle, p, pup=pup)

    """
    if (disp):
        plt.figure(11)
        fig1, axes1 = plt.subplots(2, 2)  # fig and 2 x 2 nparray of axes
        plt.figure(12)
        fig2, axes2 = plt.subplots(1, 2)  # fig and 2 x 2 nparray of axes
        plt.figure(13)
        fig3, axes3 = plt.subplots(2, 2)  # fig and 2 x 2 nparray of axes

    pix = (np.pi * _pix) / (3600 * 180)  #Conversion en rad
    z = l / (n * pix)  #pix (m)
    D = Dtel / z  #Diametre miroir (pix)
    D = np.round(D)  # on veut un entier coute que coute
    if (D % 2 == 1):
        D += 1
    r = D / 2.0
    larg = l / Dtel / pix  #Lamda/D en pix

    E = np.zeros((n, n), dtype=np.complex64)  #Field initialisation

    # on choisit les tqblequx x et y de cette fqcon, pour que, par
    # convention, on adresse les tableaux selon tab[x,y]
    y = np.tile(np.arange(n) - n / 2, (n, 1))
    x = y.T
    d2 = x * x + y * y
    if (pup is False):
        pup = (d2 < r * r) & (d2 > (obs * r)**2)  #Pupille et obstruction centrale

    ppup = pup[int(n / 2 - D / 2) + 1:int(n / 2 + D / 2) + 1,
               int(n / 2 - D / 2) + 1:int(n / 2 + D / 2) + 1]

    #Pangle = int(3*D/4)+1 #Pyramid angle
    #Pangle = int((D/2)*3)+1 #Pyramid angle
    K = (2 * np.pi * Pangle) / n
    pyr = K * (np.abs(x) + np.abs(y))  #Pyramide
    pyr = np.fft.fftshift(pyr)
    theta = (np.arange(N)) * 2 * np.pi / N  #Angle de modulation

    a = mod * np.cos(theta)
    b = mod * np.sin(theta)
    PUPIM = 0
    PSF = 0
    PUPYR = 0
    # Coeff multiplicatif des fonctions x et y pour etre en unites lam/D
    magic = larg * (2 * np.pi / n)
    print("lambda/D=%f pixels" % larg)
    x = x * magic
    y = y * magic
    if (disp):
        axes1[0, 0].matshow(p)  # Phase
        axes1[0, 0].set_title("Phase")
        axes1[0, 1].matshow(p * pup)  # Phase in pupil
        axes1[0, 1].set_title("Phase in HR pupil")
        # PSF on top of Pyramid
        PSF0 = np.fft.fftshift(np.abs(np.fft.fft2(pup * Amp * np.exp(1j * (p))))
                               **2)  #PSF not modulated
        axes1[1, 1].matshow(PSF0, cmap="gist_earth")
        axes1[1, 1].plot((PSF0.shape[0] / 2, PSF0.shape[0] / 2), (0, PSF0.shape[0]),
                         color="blue")
        axes1[1, 1].plot((0, PSF0.shape[0]), (PSF0.shape[0] / 2, PSF0.shape[0] / 2),
                         color="blue")
        axes1[1, 1].set_title("PSF on top of Pyramid")

    # Modulation LOOP !!!
    for i in range(N):
        phi = a[i] * x + b[i] * y  #Modulation Tip/Tilt
        E = pup * Amp * np.exp(1j * (p + phi))
        EPSF = np.fft.fft2(E)
        PSF += np.abs(EPSF)**2  #PSF
        PUPYR = np.abs(np.fft.ifft2((EPSF * np.exp(1j * pyr))))**2
        PUPIM += PUPYR

        if (disp):
            axes1[1, 1].scatter(a[i] * larg + PSF0.shape[0] / 2,
                                b[i] * larg + PSF0.shape[1] / 2, color="red", marker='x',
                                s=20)  # Modulation points
            axes1[1, 0].scatter(a[i], b[i], color="red", marker='x',
                                s=20)  # Modulation points
            diffractionCircle = pylab.Circle((a[i], b[i]), radius=0.5, alpha=0.5)
            axes1[1, 0].add_patch(diffractionCircle)
            axes1[1, 0].set_title("Modulation points (in lambda/D)")
        print("Modulation iter=", i, "/", N)
        print("min=", np.min(PSF), " max=", np.max(PSF))

    PSF = np.fft.fftshift(PSF)
    #    pli(PUPIM,win=1) # pli

    Nph = 10**(-0.4 * mv + 3)  #photons/cm^2/s/Angstrom
    Nph *= 3000 * 1 / 200. * (np.pi * (Dtel * 100 / 2.)**2 - np.pi *
                              (0.28 * Dtel * 100 / 2.)**2)  #photons

    S = np.sum(PUPIM)
    PUPIM /= S
    PUPIM *= Nph
    #n = nrebin
    PUPIM0 = np.copy(PUPIM)

    if (noise):
        PUPIM = np.random.poisson(PUPIM)  #Bruit de photons
        PUPIM = PUPIM * 1.0 + np.random.normal(
                0, 0.5, (n, n))  #Bruit lecture gaussien e- rms/pix
    """
    IA1 = PUPIM[int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1 , int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1]
    IB1 = PUPIM[int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1 , int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1]
    IC1 = PUPIM[int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1 , int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1]
    ID1 = PUPIM[int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1 , int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1]
    """
    # nouvelle version (corrected by Rico & Fab)
    x1 = int(n / 2 - D / 2 - Pangle)
    y1 = x1 + int(D)
    x2 = n - y1
    y2 = n - x1
    IA = PUPIM[x1:y1, x1:y1]
    IB = PUPIM[x2:y2, x1:y1]
    IC = PUPIM[x1:y1, x2:y2]
    ID = PUPIM[x2:y2, x2:y2]
    """
    IA = PUPIM_compass[int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1 , int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1]
    IB = PUPIM_compass[int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1 , int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1]
    IC = PUPIM_compass[int(n/2-D/2-Pangle)+1 : int(n/2+D/2-Pangle)+1 , int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1]
    ID = PUPIM_compass[int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1 , int(n/2-D/2+Pangle)+1 : int(n/2+D/2+Pangle)+1]
    """

    PUPIM_LR = rebin(PUPIM0, (PUPIM0.shape[0] / nrebin, PUPIM0.shape[1] / nrebin))

    Itot = IA + IB + IC + ID
    Sx = ppup * ((IA + IC) - (IB + ID)) / (Itot)
    Sy = ppup * ((IA + IB) - (IC + ID)) / (Itot)
    dwx = mod * np.sin(0.5 * np.pi * Sx)
    dwy = mod * np.sin(0.5 * np.pi * Sy)  # en unite lam/D

    # Low resolution:
    #ppup_LR = rebin(ppup, (ppup.shape[0]/nrebin, ppup.shape[1]/nrebin))
    xx1 = int(n / 2 - D / 2 - Pangle) // int(nrebin)
    yy1 = int(xx1 + int(D) / nrebin)
    xx2 = int(n / nrebin - yy1)
    yy2 = int(n / nrebin - xx1)
    IA_LR = PUPIM_LR[xx1:yy1, xx1:yy1]
    IB_LR = PUPIM_LR[xx2:yy2, xx1:yy1]
    IC_LR = PUPIM_LR[xx1:yy1, xx2:yy2]
    ID_LR = PUPIM_LR[xx2:yy2, xx2:yy2]
    Itot_LR = IA_LR + IB_LR + IC_LR + ID_LR
    Sx_LR = pupsub * ((IA_LR + IC_LR) - (IB_LR + ID_LR)) / (Itot_LR)
    Sy_LR = pupsub * ((IA_LR + IB_LR) - (IC_LR + ID_LR)) / (Itot_LR)
    dwx_LR = mod * np.sin(0.5 * np.pi * Sx_LR)
    dwy_LR = mod * np.sin(0.5 * np.pi * Sy_LR)  # en unite lam/D

    ind = np.where(pupsub.flatten())[0]
    slopesx = dwx_LR.flatten()[ind]
    slopesy = dwy_LR.flatten()[ind]

    slopes = np.zeros(slopesx.shape[0] * 2)
    slopes[:len(ind)] = slopesx
    slopes[len(ind):] = slopesy

    if (disp):
        axes2[0].matshow(PUPIM0, cmap="gist_earth")
        axes2[0].set_title("PYRAMID image High Resolution")
        axes2[1].matshow(PUPIM_LR, cmap="gist_earth")
        axes2[1].set_title("PYRAMID image Low Resolution")
        #divider1 = make_axes_locatable(axes2[0])
        #cax1 = divider1.append_axes("right", size="20%", pad=0.05)
        #cbar1 = plt.colorbar(PUPIM0, cax=cax1)
        #divider2 = make_axes_locatable(axes2[1])
        #cax2 = divider2.append_axes("right", size="20%", pad=0.05)
        #cbar2 = plt.colorbar(PUPIM_LR, cax=cax2)

        axes3[0, 0].matshow(dwx, cmap="gist_earth")
        #divider3 = make_axes_locatable(axes3[0, 0])
        #cax3 = divider3.append_axes("right", size="20%", pad=0.05)
        #cbar3 = plt.colorbar(dwx, cax=cax3)
        axes3[0, 1].matshow(dwy, cmap="gist_earth")
        #divider4 = make_axes_locatable(axes3[0, 1])
        #cax4 = divider4.append_axes("right", size="20%", pad=0.05)
        #cbar4 = plt.colorbar(dwy, cax=cax4)
        axes3[0, 0].set_title("GRAD (x) High Resolution")
        axes3[0, 1].set_title("GRAD (y) High Resolution")

        axes3[1, 0].matshow(dwx_LR, cmap="gist_earth")
        #divider5 = make_axes_locatable(axes3[1, 0])
        #cax5 = divider5.append_axes("right", size="20%", pad=0.05)
        #cbar5 = plt.colorbar(dwx_LR, cax=cax5)
        axes3[1, 1].matshow(dwy_LR, cmap="gist_earth")
        #divider6 = make_axes_locatable(axes3[1, 1])
        #cax6 = divider6.append_axes("right", size="20%", pad=0.05)
        #cbar6 = plt.colorbar(dwy_LR, cax=cax5)
        axes3[1, 0].set_title("GRAD (x) Low Resolution")
        axes3[1, 1].set_title("GRAD (y) Low Resolution")

        fig1.tight_layout(pad=0.1, h_pad=0.5, w_pad=0.5, rect=None)
        fig2.tight_layout(pad=0.1, h_pad=0.5, w_pad=0.5, rect=None)
        fig3.tight_layout(pad=0.1, h_pad=0.5, w_pad=0.5, rect=None)
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
    return PUPIM0, PUPIM, dwx, dwy, IA, IB, IC, ID, slopes, pyr


def compPyrCOMPASS(wao):
    #imhrCOMPASS = wao.supervisor.wfs.get_hrimg_pyr(0)
    phasehrCOMPASS = wao.supervisor._sim.wfs.get_phase(0)
    pupCOMPASS = wao.supervisor.getMpupil()
    #phaselrCOMPASS = wao.supervisor.wfs.get_pyrimg(0)
    ncompass = phasehrCOMPASS.shape[0]
    n = wao.supervisor._sim.wfs.get_pyrimghr(0).shape[0]
    nrebin = wao.supervisor._sim.config.p_wfs0._nrebin
    pup_sep = wao.supervisor._sim.config.p_wfs0.pyr_pup_sep
    mod = wao.supervisor._sim.config.p_wfs0.pyr_ampl
    N = wao.supervisor._sim.config.p_wfs0.pyr_npts
    Dtel = wao.supervisor._sim.config.p_tel.diam
    obs = wao.supervisor._sim.config.p_tel.cobs
    l = wao.supervisor._sim.config.p_wfs0.Lambda * 1e-6
    _pix = wao.supervisor._sim.config.p_wfs0._qpixsize
    Amp = 1.
    mv = 10.
    Pangle = pup_sep * nrebin
    p = np.zeros((n, n))
    pup = p.copy()
    p[(n - ncompass) // 2:(n + ncompass) // 2, (n - ncompass) // 2:(n + ncompass) //
      2] = phasehrCOMPASS * 2 * np.pi / l * 1e-6
    pup[(n - ncompass) // 2:(n + ncompass) // 2, (n - ncompass) // 2:(n + ncompass) //
        2] = pupCOMPASS
    pupsub = wao.supervisor._sim.config.p_wfs0._isvalid[1:-1, 1:-1]
    PUPIM0, PUPIM, dwx, dwy, IA, IB, IC, ID, slopes, pyrmask = pyr_analysis(
            n, mod, N, Dtel, obs, nrebin, l, _pix, Amp, mv, Pangle, p, pupsub, pup=pup)
    return PUPIM0, PUPIM, dwx, dwy, IA, IB, IC, ID, slopes, pyrmask


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
