#!/usr/local/bin/python3.6
# encoding: utf-8
'''
Created on 1 aout 2017

@author: fferreira
'''
import numpy as np


def rebin(a, shape):
    """
    TODO: docstring

    """

    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    return a.reshape(sh).mean(-1).mean(1)


def fft_goodsize(s):
    """find best size for a fft from size s

    :parameters:
         s: (int) size
    """
    return 2**(int(np.log2(s)) + 1)


def bin2d(data_in, binfact):
    """
    Returns the input 2D array "array", binned with the binning factor "binfact".
    The input array X and/or Y dimensions needs not to be a multiple of
    "binfact"; The final/edge pixels are in effect replicated if needed.
    This routine prepares the parameters and calls the C routine _bin2d.
    The input array can be of type int, float or double.
    Last modified: Dec 15, 2003.
    Author: F.Rigaut
    SEE ALSO: _bin2d

    :parmeters:
        data_in: (np.ndarray) : data to binned

        binfact: (int) : binning factor

    """
    if (binfact < 1):
        raise ValueError("binfact has to be >= 1")

    nx = data_in.shape[0]
    ny = data_in.shape[1]
    fx = int(np.ceil(nx / float(binfact)))
    fy = int(np.ceil(ny / float(binfact)))

    data_out = np.zeros((fx, fy), dtype=data_in.dtype)

    for i1 in range(fx):
        for j1 in range(fy):
            for i2 in range(binfact):
                for j2 in range(binfact):
                    i = i1 * binfact + i2
                    j = j1 * binfact + j2
                    if (i >= nx):
                        i = nx - 1
                    if (j >= ny):
                        j = ny - 1
                    data_out[i1, j1] += data_in[i, j]

    return data_out


def MESH(Range, Dim):
    """
    TODO: docstring

    """

    last = (0.5 * Range - 0.25 / Dim)
    step = (2 * last) // (Dim - 1)

    return np.tile(np.arange(Dim) * step - last, (Dim, 1))


def pad_array(A, N):
    """
    TODO: docstring

    """

    S = A.shape
    D1 = (N - S[0]) // 2
    D2 = (N - S[1]) // 2
    padded = np.zeros((N, N))
    padded[D1:D1 + S[0], D2:D2 + S[1]] = A
    return padded


def dist(dim, xc=-1, yc=-1):
    """
    TODO: docstring

    """

    if (xc < 0):
        xc = int(dim / 2.)
    else:
        xc -= 1.
    if (yc < 0):
        yc = int(dim / 2.)
    else:
        yc -= 1.

    dx = np.tile(np.arange(dim) - xc, (dim, 1))
    dy = np.tile(np.arange(dim) - yc, (dim, 1)).T

    d = np.sqrt(dx**2 + dy**2)
    return d


def makegaussian(size, fwhm, xc=-1, yc=-1, norm=0):
    """
    Returns a centered gaussian of specified size and fwhm.
    norm returns normalized 2d gaussian

    :param size: (int) :
    :param fwhm: (float) :
    :param xc: (float) : (optional) center position on x axis
    :param yc: (float) : (optional) center position on y axis
    :param norm: (int) : (optional) normalization

    """
    tmp = np.exp(-(dist(size, xc, yc) / (fwhm / 1.66))**2.)
    if (norm > 0):
        tmp = tmp / (fwhm**2. * 1.140075)
    return tmp
