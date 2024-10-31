#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
#
# COMPASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COMPASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COMPASS. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2024 COSMIC Team


import numpy as np

import shlex
from subprocess import Popen, PIPE  # , call
from sys import stdout
from time import sleep

# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import time
from math import factorial as fac

def clr(*figs):
    """
    THE Fab function

    clears the current figure (no arg) or specified window

    """
    if figs:
        for fig in figs:
            # fig = fig[i]
            plt.figure(num=fig)
            plt.clf()
    else:
        plt.clf()


def system(cmd, output=False):
    """
    Execute the external command
    system("ls")


    out = system("ls", out=True)
    out = system("ls -l", out=True)




    and get its stdout exitcode and stderr.
    """
    args = shlex.split(cmd)
    proc = Popen(args, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    exitcode = proc.returncode
    #
    if "\n" in out:
        out = out.split("\n")[:-1]

    for i in range(len(out)):
        print((out[i]))

    if output:
        # print("here")
        return out, exitcode, err


def pli(
    data,
    color="gist_earth",
    cmin=9998,
    cmax=9998,
    win=1,
    origin=None,
    aspect="equal",
):
    """
    plots the transpose of the data

    color maps keywords can be found in
    http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps

    """
    options = ""
    if cmin != 9998:
        exec('options += ",vmin=cmin"')

    if cmax != 9998:
        exec('options += ",vmax=cmax"')

    if color == b"yorick":
        color = "gist_earth"
    if origin is None:
        origin = ""
    if aspect != "auto":
        aspect = "'" + aspect + "'"
    else:
        aspect = "'auto'"

    exec("plt.matshow(data, aspect=" + aspect + ", fignum=win, cmap=color" + options + origin + ")")


def binning(w, footprint):
    # the averaging block
    # prelocate memory
    binned = np.zeros(w.shape[0] * w.shape[1]).reshape(w.shape[0], w.shape[1])
    # print(w)
    for i in range(w.shape[0]):
        for j in range(w.shape[1]):
            binned[i, j] = w[i, j].sum() / (footprint * footprint + 0.0)

    return binned


def minmax(tab):
    tabVect = np.reshape(tab, tab.size)
    return [np.min(tabVect), np.max(tabVect)]


def plg(
    data,
    x="",
    win=1,
    xlog=0,
    ylog=0,
    color="black",
):
    """


    color = "green"
    color = "0.71" [0-1] gray scale
    color = '#eeefff'
    See also:

    http://matplotlib.org/api/colors_api.html

    """
    fig = plt.figure(win)
    ax = fig.add_subplot(1, 1, 1)
    try:
        data.ndim
        if data.ndim > 1:
            print(("Warning %dD dimensions. Cannot plot data. Use pli instead. " % data.ndim))
    except BaseException:
        return
    if x == "":
        ax.plot(data, color=color)
    else:
        ax.plot(x, data, color=color)

    if xlog == 1:
        ax.set_xscale("log")
    else:
        ax.set_xscale("linear")

    if ylog == 1:
        ax.set_yscale("log")
    else:
        ax.set_yscale("linear")
    fig.show()
    return fig, ax


def zcen(data):
    data = np.array(data)
    if len(data.shape) > 1:
        print("oups zcen with dims > 1 not coded yet...")
        return 0
    tmp = tmp2 = []
    for i in range(len(data) - 1):
        tmp = (float(data[i]) + float(data[i + 1])) / 2.0
        tmp2 = np.append(tmp2, tmp)
    return tmp2


def getValidSubapArray(nssp, rext, rint, return2d=False):
    # The Grata case, tip-tilt sensor only.
    if nssp == 1:
        return [1]
    # to avoid some bug that eliminates useful central subapertures when
    # obs=0.286
    if (nssp == 7) and (rint > 0.285 and rint < 0.29):
        rint = 0.285
        print("cas particulier")
    x = zcen(np.linspace(-1, 1, num=nssp + 1))
    xx = []
    for i in range(nssp):
        xx = np.hstack((xx, x))
    x = np.reshape(xx, (nssp, nssp))
    y = np.transpose(x)
    r = np.sqrt(x * x + y * y)
    valid2dext = (r < rext) * 1
    valid2dint = (r >= rint) * 1
    valid2d = valid2dint * valid2dext

    if return2d:
        return valid2d
    else:
        valid = np.reshape(valid2d, [nssp * nssp])

    return valid.tolist()


"""
def plsh(slopesvector,  nssp=14,  rmax=0.98, obs=0, win=1, invertxy=False):

    tmp = getValidSubapArray( nssp, rmax, obs);
    X,Y = meshgrid(np.linspace(-1, 1, nssp), np.linspace(-1, 1, nssp))
    vx = np.zeros([nssp*nssp])
    vy = np.zeros([nssp*nssp])
    hart = where(tmp)[0]
    vx.flat[hart] = slopesvector.flat[0:len(slopesvector)/2]
    vy.flat[hart] = slopesvector.flat[len(slopesvector)/2+1:]
    vx = vx.reshape([nssp, nssp])
    vy = vy.reshape([nssp, nssp])

    figure(num=win)
    if(invertxy):
        Q = quiver(X,Y, vy, vx)
    else:
        Q = quiver(X,Y, vx, vy)
    #qk = quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
    l,r,b,t = axis()
    dx, dy = r-l, t-b
    axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy]) # MUST DO OTHERWISE THE AUTOSCALE CAN MISS SOME ARROWS
    #title('Minimal arguments, no kwargs')
"""


def plpyr(slopesvector, validArray):
    """
    wao.config.p_wfss[0]._isvalid
    """
    nslopes = slopesvector.shape[0] / 2
    x, y = np.where(validArray.T)
    plt.quiver(x, y, slopesvector[0:nslopes], slopesvector[nslopes:])


def plsh(
    slopesvector,
    nssp,
    validint,
    sparta=False,
    invertxy=False,
    returnquiver=False,
):
    """
    <slopesvector> is the input vector of slopes
    <nssp> is the number of subapertures in the diameter of the pupil
    <validint> is the normalized diameter of central obscuration (between 0 and 1.00)
    <sparta> when==1, slopes are ordered xyxyxyxy...
             when==0, slopes are xxxxxxyyyyyyy
    <xy> when==1, swap x and y. Does nothing special when xy==0.

    The routine plots a field vector of subaperture gradients defined in
    vector <slopesvector>.
    The routine automatically adjusts/finds what are the valid subapertures
    for plotting, depending on the number of elements in <slopesvector>. Only the
    devalidated subapertures inside the central obscuration cannot be
    known, that’s why <validint> has to be passed in the argument list.

    """
    nsub = slopesvector.shape[0] // 2
    x = np.linspace(-1, 1, nssp)
    x, y = np.meshgrid(x, x)
    r = np.sqrt(x * x + y * y)
    # defines outer and inner radiuses that will decide of validity of subapertures
    # inner radius <validint> is passed as an argument.
    # outer one will be computed so that it will match the number of
    # subapertures in slopesvector
    rorder = np.sort(r.reshape(nssp * nssp))
    # number of subapertures not valid due to central obscuration
    ncentral = nssp * nssp - np.sum(r >= validint, dtype=np.int32)
    # determine value of external radius so that the test (validint < r < validext)
    # leads to the correct number of subapertures
    validext = rorder[ncentral + nsub]
    # get the indexes of valid subapertures in the nsspxnssp map
    valid = (r < validext) & (r >= validint)
    ivalid = np.where(valid)
    # feeding data <slopesvector> into <vv>
    vx = np.zeros([nssp, nssp])
    vy = np.zeros([nssp, nssp])
    if sparta is False:
        # Canary, compass, etc..  slopes ordered xxxxxxxyyyyyyy
        vy[ivalid] = slopesvector[0:nsub]
        vx[ivalid] = slopesvector[nsub:]
    else:
        # SPARTA case, slopes ordered xyxyxyxyxyxy...
        vx[ivalid] = slopesvector[0::2]
        vy[ivalid] = slopesvector[1::2]
    if invertxy is True:
        # swaps X and Y
        tmp = vx
        vx = vy
        vy = tmp
    if returnquiver:
        return x, y, vx, vy
    else:
        plt.quiver(x, y, vx, vy, pivot="mid")


def pl3d(im):
    """
    ir = pyfits.get_data("/home/fvidal/data/Run2015/June2015_27_onsky/ir/ir_2015-06-28_06h27m40s_script44_gain.fits")

    JAMAIS TESTEE !!!!!!!!!!!!!!

    """
    X = np.arange(-5, 5, 0.25)
    Y = np.arange(-5, 5, 0.25)
    X, Y = np.meshgrid(X, Y)
    Z = im
    plt.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet)
    plt.show()


def FFThz(signal, fe, freq=0):
    """PSD = FFThz( signal, fe )   OU  f = FFThz( 1024, fe, freq=1 )
    On the first form, returns the power spectral density of signal.
    If signal has units 'u', the PSD has units 'u^2/Hz'.
    The frequency axis can be get by using the keyword freq=1."""
    if freq == 1:
        n = signal.size
        d = np.linspace(0, fe, n + 1)[0 : n / 2 + 1]
        return d[1:]
    else:
        n = signal.size
        d = np.abs(np.fft.fft(signal))[0 : n / 2 + 1]
        d = d**2 / (fe * n / 2)
        d[n / 2] /= 2
        return d[1:]


def computePSD(zerall, fe, izerNum, wfsNum):
    if np.isscalar(wfsNum):
        wfsNum = [wfsNum]

    for ii in wfsNum:
        PSD = FFThz(zerall[ii][izerNum, :], fe)

    PSD /= len(wfsNum)
    if len(wfsNum) > 1:
        ff = FFThz(zerall[wfsNum][izerNum, :], fe, freq=1)
    else:
        ff = FFThz(zerall[wfsNum[0]][izerNum, :], fe, freq=1)

    return PSD, ff


def countExample(seconds):
    for i in range(1, int(seconds)):
        stdout.write("\r%d" % i)
        stdout.flush()
        sleep(1)
    stdout.write("\n")


def plotSubapRectangles(pup, isvalid, istart, jstart):
    fig = plt.matshow(pup)
    pdiam = istart[1] - istart[0]
    for i in istart:
        for j in jstart:
            if isvalid[i // pdiam, j // pdiam]:
                color = "green"
            else:
                color = "red"
            fig.axes.add_patch(
                plt.Rectangle((i - 0.5, j - 0.5), pdiam, pdiam, fill=False, color=color)
            )

def gaussian(x, a, x0, sigma, b):

    """
    This function dives a gaussian profile.

    :param x: 1d array of distances to central pixel [μm]
    :param a: signal peak value
    :param x0: position of the center [μm]
    :param sigma: variance
    :return:
    """

    gaus = a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + b
    return gaus


def min_array(array):
    """
    Compute min of array excluding nan values
    """

    return np.min(array[np.where(np.isnan(array) == False)])


def max_array(array):
    """
    Compute max of array excluding nan values
    """

    return np.max(array[np.where(np.isnan(array) == False)])


def wait_until(predicate, timeout, period):
    mustend = time.time() + timeout

    while time.time() < mustend:
        if predicate:
            return True
        else:
            time.sleep(period)

    return print("Reach Time Out = {} s".format(timeout))


def zernike(n_pix, m, n):
    """
    Noll Zernike term expansion

    :param n_pix: <int> size of zernike support
    :param m: <int> angular meridional frequency
    :param n: <int> radial order ; n>0, n >= m
    :return:
    """

    x = np.linspace(-1, 1, n_pix)
    y = np.linspace(-1, 1, n_pix)
    xx, yy = np.meshgrid(x, y)

    ro, theta = cart2polar(xx, yy)
    w = np.where(ro > 1)
    w2 = np.where(ro <= 1)
    # ro[w] = 0

    if m > n or n < 0:
        print('Error Zernike index should verify: n >= m and n > 0')
        return 0

    elif (m == 0) and (n == 0):
        z = 2 * np.cos(0 * theta)

    elif m >= 0:
        if (n - m) % 2 == 0:
            R = sum(
                (((-1) ** kk * fac(n - kk)) / (fac(kk) * fac((n + m) / 2 - kk) * fac((n - m) / 2 - kk))) * ro ** (
                            n - 2 * kk) for kk in range(int((n - m) / 2) + 1))
        else:
            R = 0
        z = R * np.cos(m * theta)

    else:
        mb = np.abs(m)
        if (n - mb) % 2 == 0:
            R = sum(
                (((-1) ** kk * fac(n - kk)) / (fac(kk) * fac((n + mb) / 2 - kk) * fac((n - mb) / 2 - kk))) * ro ** (
                            n - 2 * kk) for kk in range(int((n - mb) / 2) + 1))
        else:
            R = 0
        z = R * np.sin(mb * theta)

    # z[w] = 0

    # print("mean map = {}", np.mean(z[w2]))

    return z/2


def list_of_zernike(n_pix):
    """

    :param n_pix:
    :return:
    """

    m = np.array([0, 1, -1, 0, 2, -2, 1, -1, 3, -3, 0, 2, -2, 4, -4, 1, -1, 3, -3, 5, -5])
    n = np.array([0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5])
    name = ['piston', 'tip', 'tilt', 'defocus', 'astigmatism', 'astigmatism', 'coma', 'coma', 'trefoil', 'trefoil',
            'spherical', 'astigmatism 2', 'astigmatism 2', 'quadrafoil', 'quadrafoil', 'coma 2', 'coma2', 'trefoil 2',
            'trefoil 2', 'pentafoil', 'pentafoil']

    z = []

    for i in np.arange(len(m)):
        z.append(zernike(n_pix, m[i], n[i]))

    return z


def m_n_zernike(n_z):
    m = np.array([0, 1, -1, 0, 2, -2, 1, -1, 3, -3, 0, 2, -2, 4, -4, 1, -1, 3, -3, 5, -5, 0, 2, -2, 4, -4, 6, -6, 1, -1,
                  3, -3, 5, -5, 7, -7, 0, 2, -2, 4, -4, 6, -6, 8, -8])
    n = np.array([0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7,
                  7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8])

    return m[n_z], n[n_z]


def cart2polar(xx, yy):
    """
    Convert cartesian coordinates into polar coordinates
    :param xx: 2D-square array of x linear coordinates
    :param yy: 2D-square array of y linear coordinates - should have same dimensions as xx
    :return: r, xx/yy like array, normalized radius
            theta, xx/yy like array, angle in radian
    """

    phi = np.arctan2(-yy, -xx)
    theta = phi - np.min(phi)
    r = np.sqrt(xx ** 2 + yy ** 2)

    return r, theta


def rad2arcsec(x):
    return x * 180 / np.pi * 3600


def embed_f(matrix, factor):
    """
    Embed a Matrix in another matrix factor times larger
    :param matrix: matrix to embedded
    :param factor: factor for new matrix
    :return:
    """

    n_matrix = len(matrix)
    n_embed = int(factor * len(matrix))
    matrix_embed = np.zeros((n_embed, n_embed), dtype=complex)
    s = int(n_embed / 2 - n_matrix / 2)
    e = int(n_embed / 2 + n_matrix / 2)
    matrix_embed[s:e, s:e] = matrix

    return matrix_embed


def embed_len(matrix, n_embed):
    """
    Embed a matrix in another one with length n_embed
    :param matrix: matrix to be embedded
    :param n_embed: length of the created matrix
    :return:
    """

    n_matrix = len(matrix)
    matrix_embed = np.zeros((n_embed, n_embed), dtype=complex)
    s = int(n_embed / 2 - n_matrix / 2)
    e = int(n_embed / 2 + n_matrix / 2)
    matrix_embed[s:e, s:e] = matrix

    return matrix_embed


def find_nearest(array, value):
    """
    Find index of the nearest value in an array.
    :param array:
    :param value:
    :return:
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def eq_mod(data, mod_value):
    """
    Equivalent of modulo, except it returns the closest value to zero between data%mod_value and data%-mod_value
    :param data: array of data to be processed
    :mod_value: value for modulo operation
    :return: processed data 
    Example :
        mod(0.68, 0.7) = 0.68
        eq_mod(0.68, 0.7) = 0.02
    """
    set1 = data%mod_value
    set2 = data%-mod_value

    new_data = np.zeros(data.shape)
    w1 = np.where(np.abs(set1) == np.minimum(np.abs(set1), np.abs(set2)))
    w2 = np.where(np.abs(set2) == np.minimum(np.abs(set1), np.abs(set2)))

    new_data[w1] = set1[w1]
    new_data[w2] = set2[w2]

    return new_data


def jump_counter(data, mod_value):
    """
    Number of times when data jump of more than the mod_value.
    Usefull for hopping differential piston.
    """

    set1 = data - np.floor_divide(data, mod_value) * mod_value
    set2 = data - np.floor_divide(data, - mod_value) * (- mod_value)

    # w1 = np.where(np.abs(set1) == np.minimum(np.abs(set1), np.abs(set2)))
    w2 = np.where(np.abs(set2) == np.minimum(np.abs(set1), np.abs(set2)))

    new_set = np.floor_divide(data, mod_value)
    new_set[w2] = - np.floor_divide(data, -mod_value)[w2]

    return np.sum(np.abs(np.diff(new_set, axis=0)))


def compute_pv_strehl(pv, lwfs=0.7):
    """
    pv data des 6 pistons calculés à tous les seeings, toutes les réalisations de phase, toutes les trames en boucle fermée
    lwfs longueur d'onde de l'analyseur
    """
    srpv = np.zeros(pv.shape[0])
    npetal = pv.shape[1]
    for k in range(pv.shape[0]):    # toutes les trames
        r = pv[k]
        srpv[k] = np.sum(np.cos((r[None,:] - r[:,None])*2*np.pi/lwfs)) / (npetal**2)
    stdpv = np.sqrt(-np.log(srpv))
    return srpv, stdpv


def psf_pf(dzl, n, N=1024):
    """
    dzl = D.z / lambda
    n taille du support voulu pour la tache d'Airy
    """
    P = int(np.round(N*dzl))
    x = np.linspace(-1, 1, N)
    #intermède coordonnées polaires pour calcul rayon
    xx, yy = np.meshgrid(x, x)
    r = np.sqrt(xx ** 2 + yy ** 2)
    # création pupille
    pup = np.zeros((N, N))
    pup[r<P/N] = 1
    # tache d'airy
    ft_pup = np.fft.fftshift(np.fft.fft2(pup)) / np.sum(pup)
    ft_pup_crop = ft_pup[(N-n)//2: (N+n)//2, (N-n)//2: (N+n)//2]
    return pup, np.abs(ft_pup_crop)**2


def compute_psf(pup, phi, lambda2rad=1):
    field = pup * np.exp(1j * phi * lambda2rad)
    tfield = np.fft.fftshift(np.fft.fft2(field))
    return (np.abs(tfield) / np.sum(pup))**2


def create_pupil(n_pix, d_tel, n_seg=0, d_obs=0, d_spider=0, d_spider2 = 0, form=None, d_in=None):
    """
    all distances expressed in pixels with respect to n_pix
    """
    x = np.linspace(-1, 1, n_pix)
    xx, yy = np.meshgrid(x, x)
    r = np.sqrt(xx**2 + yy**2)
    from scipy.ndimage import rotate
    my_pup = (r < d_tel/n_pix) * (r >= d_obs/n_pix) * 1
    if d_spider2 != 0:
        my_pup[int(n_pix//2 - d_spider//2):int(n_pix//2 + d_spider2//2), n_pix//2:] = 0
    else:
        my_pup[int(n_pix//2 - d_spider//2):int(n_pix//2 + d_spider//2), n_pix//2:] = 0
    my_pup_rot = my_pup
    for i in range(n_seg):
        my_pup_rot = rotate(my_pup_rot, 360 / n_seg, reshape=False)
        my_pup *= my_pup_rot

    if form is not None:
        my_pup[np.where(r<d_in/n_pix)] = ipup[np.where(r<d_in/n_pix)]

    return my_pup
