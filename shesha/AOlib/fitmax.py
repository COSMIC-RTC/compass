#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 16:13:34 2019

@author: micado
"""

import numpy as np


def fitmax(im, xmax, ymax, full=False, size=3):
    """
    <im> : the image
    <xmax>, <ymax> : coordinates of the maximum in the image
    <full=> : True / False(default). When False, returns only (x0, y0).
              When True, returns (x0, y0, fwhmx, fwhmy, valmax, t)

    Usage:
        posmax = fitmax(image,xmax,ymax)

    Fits the square of 3x3 pixels of the image centered on pixel
    of coordinates (xmax,ymax), by a paraboloid of equation
    z = A(x-x0)**2 + B(y-y0)**2 + C(x-x0).(y-y0) + D
    and returns the tuple (x0, y0).

    Tue Oct  6 08:09:20 CEST 2009
    Correction pgm : il n'y avait pas de terme croise C.(x-x0).(y-y0) dans le
    modèle, ce qui interdit le fit par un paraboloide oblique... c'est con.

    (%i45) linsolve([-C*y0-2*A*x0=p(4),-C*x0-2*B*y0=p(5)],[x0,y0]);
                         p(5) C - 2 p(4) B         p(4) C - 2 p(5) A
    (%o45)       [x0 = - -----------------, y0 = - -----------------]
                             2                         2
                            C  - 4 A B                C  - 4 A B

    """
    z = im[xmax - 1:xmax + 2, ymax - 1:ymax + 2].flatten()
    v = np.array(
            [[18, -36, 18, 18, -36, 18, 18, -36, 18],
             [18, 18, 18, -36, -36, -36, 18, 18, 18], [27, 0, -27, 0, 0, 0, -27, 0, 27],
             [-18, 0, 18, -18, 0, 18, -18, 0, 18], [-18, -18, -18, 0, 0, 0, 18, 18, 18],
             [-12, 24, -12, 24, 60, 24, -12, 24, -12]]) / 108.00

    # fit coeffs of p[0]*x**2+p[1]*y**2+p[2]*x*y+p[3]*x+p[4]*y+p[5]
    p = np.dot(v, z)
    A = p[0]
    B = p[1]
    C = p[2]
    denom = C**2. - 4 * A * B
    if (denom == 0):
        x0, y0 = 0.0, 0.0
    else:
        y0 = (2. * B * p[3] - p[4] * C) / denom
        x0 = (2. * A * p[4] - p[3] * C) / denom

    D = p[5] - (x0 * y0 * C + A * x0**2. + B * y0**2.)

    x0 += xmax
    y0 += ymax
    if (full is True):
        valmax = D  # top of the polynom
        if ((B - A) == 0):
            t = np.pi() / 4
            if (C == 0):
                t = 0.
        else:
            t = np.arctan(C / (B - A)) / 2
        AA = B * np.sin(t)**2 - C * np.cos(t) * np.sin(t) + A * np.cos(t)**2
        BB = A * np.sin(t)**2 + C * np.cos(t) * np.sin(t) + B * np.cos(t)**2
        fwhmx = 1.66 * np.sqrt(-D / 2. / AA)
        fwhmy = 1.66 * np.sqrt(-D / 2. / BB)

        return (x0, y0, fwhmx, fwhmy, valmax, t)  # t = angle
    else:
        return (x0, y0)


def fitmax_cross1d(im, xmax, ymax):
    x0, xvtop = fit3ptsParab(im[xmax - 1:xmax + 2, ymax])
    y0, yvtop = fit3ptsParab(im[xmax, ymax - 1:ymax + 2])
    '''
    u = np.array([-1,0,1])
    y0, yvtop = fit3ptsParab( im[(xmax+u, ymax+u)] )
    y0, yvtop = fit3ptsParab( im[(xmax+u, ymax+u)] )
    '''
    x0 = x0 + xmax
    y0 = y0 + ymax
    vtop = (xvtop + yvtop) / 2.0
    return x0, y0, vtop


def fit3ptsParab(abc):
    c = abc[1]
    a2 = abc[0] + abc[2] - 2 * c
    b = (abc[2] - abc[0]) / 2.0
    if a2 == 0:
        xmax = 0.0
    else:
        xmax = -b / a2
    vtop = ((a2 / 2 * xmax) + b) * xmax + c
    return xmax, vtop
