# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 15:42:43 2016
Function for DM Python kl
@author: translate by sdurand
Compass Yorick translation
"""

# import
import numpy as np
from scipy import interpolate

#function
#__________________________________________________________________________


def make_radii(cobs, nr):
    d = (1. - cobs * cobs) / nr
    rad2 = cobs**2 + d / 16. + d * (np.arange(nr, dtype=np.float32))
    radp = np.sqrt(rad2)
    return radp


#__________________________________________________________________________


def make_kernels(cobs, nr, radp, funct, outscl):
    #make kernels
    #DOCUMENT res=make_kernels(cobs,nr,rad,funct=,outscl=)
    #This routine generates the kernel used to find the KL modes.
    #The  kernel constructed here should be simply a discretization
    #of the continuous kernel. It needs rescaling before it is treated
    #as a matrix for finding  the eigen-values. The outer scale
    #should be in units of the diameter of the telescope.
    nth = 5 * nr
    kers = np.zeros((nth, nr, nr), dtype=np.float32)
    cth = np.cos((np.arange(nth, dtype=np.float32)) * (2. * np.pi / nth))
    dth = 2. * np.pi / nth
    fnorm = -1. / (2 * np.pi * (1. - cobs**2)) * 0.5
    #the 0.5 is to give  the r**2 kernel, not the r kernel
    for i in range(nr):
        for j in range(i + 1):
            te = 0.5 * np.sqrt(
                    radp[i]**2 + radp[j]**2 - (2 * radp[i] * radp[j]) * cth)
            #te in units of the diameter, not the radius
            if (funct == b"kolmo"):
                #DOCUMENT var=kolstf(dvec)
                #This routine returns the kolmogorov phase variance at spatial
                #dimension (inverse of the spatial frequency) dvec
                te = 6.88 * te**(5. / 3.)

            elif (funct == b"karman"):
                #DOCUMENT var=kolstf(dvec)
                #This routine returns the Von Karman phase variance at spatial
                #dimension (inverse of the spatial frequency) dvec. Same as kolstf
                #but with a correcting factor to account for the outter scale.
                #The latter should be in units of telescope diameter
                te = 6.88 * te**(5. / 3.) * (
                        1 - 1.485 * (te / outscl)**(1. / 3.) + 5.383 *
                        (te / outscl)**(2) - 6.281 * (te / outscl)**(7. / 3.))

            else:

                raise TypeError("kl funct error")

            f = np.fft.fft(te, axis=-1)
            kelt = fnorm * dth * np.float32(f.real)
            kers[:, i, j] = kers[:, j, i] = kelt
    return kers


#__________________________________________________________________________
#__________________________________________________________________________


def piston_orth(nr):
    s = np.zeros((nr, nr), dtype=np.float32)
    for j in range(nr - 1):
        rnm = 1. / np.sqrt(np.float32((j + 1) * (j + 2)))
        s[0:j + 1, j] = rnm
        s[j + 1, j] = -1 * (j + 1) * rnm

    rnm = 1. / np.sqrt(nr)
    s[:, nr - 1] = rnm
    return s


#__________________________________________________________________________
#__________________________________________________________________________


def make_azimuth(nord, npp):
    #DOCUMENT piston_orth(nr)

    azbas = np.zeros((np.int32(1 + nord), npp), dtype=np.float32)
    th = np.arange(npp, dtype=np.float32) * (2. * np.pi / npp)

    azbas[0, :] = 1.0
    for i in np.arange(1, nord, 2):
        azbas[np.int32(i), :] = np.cos((np.int32(i) / 2 + 1) * th)
    for i in np.arange(2, nord, 2):
        azbas[np.int32(i), :] = np.sin((np.int32(i) / 2) * th)

    return azbas


#__________________________________________________________________________
#__________________________________________________________________________


def radii(nr, npp, cobs):
    #r = radii(nr,npp,cobs)
    #DOCUMENT res=radii(NumberOfR,NumberOfPhi,Dim)
    #This routine generates an nr x npp array with npp copies of the
    #radial coordinate array. Radial coordinate span the range from
    #r=cobs to r=1 with successive annuli having equal areas (ie, the
    #area between cobs and 1 is divided into nr equal rings, and the
    #points are positioned at the half-area mark on each ring). There
    #are no points on the border.

    r2 = cobs**2 + (np.arange(nr, dtype=np.float) + 0.) / nr * (1.0 - cobs**2)
    rs = np.sqrt(r2)
    r = np.transpose(np.tile(rs, (npp, 1)))

    return r


#__________________________________________________________________________
#__________________________________________________________________________


def polang(r):
    #p = polang(r)
    #DOCUMENT res=polang(RadialCoordArray)
    #This routine generates an array with the same dimensions as r,
    #but containing the azimuthal values for a polar coordinate system.

    s = r.shape
    nr = s[0]
    np1 = s[1]
    phi1 = np.arange(np1, dtype=np.float) / float(np1) * 2. * np.pi
    p1, p2 = np.meshgrid(np.ones(nr), phi1)
    p = np.transpose(p2)

    return p


#__________________________________________________________________________
#__________________________________________________________________________


def setpincs(ax, ay, px, py, cobs):
    #func setpincs(ax,ay,px,py,cobs,&pincx,&pincy,&pincw)
    #DOCUMENT setpincs(ax,ay,px,py,cobs,&pincx,&pincy,&pincw)
    #This routine determines a set of squares for interpolating
    #from cartesian to polar coordinates, using only those points
    #with cobs < r < 1
    #SEE ALSO : pcgeom

    s = ax.shape
    nc = s[0]
    #s = px.shape# not used
    #nr = s[0]# not used
    #npp = s[1]# not used
    dcar = (ax[nc - 1, 0] - ax[0, 0]) / (nc - 1)
    ofcar = ax[0, 0]
    rlx = (px - ofcar) / dcar
    rly = (py - ofcar) / dcar
    lx = np.int32(rlx)
    ly = np.int32(rly)
    shx = rlx - lx
    shy = rly - ly

    pincx = np.zeros((4, lx.shape[0], lx.shape[1]))
    pincx[[1, 2], :, :] = lx + 1
    pincx[[0, 3], :, :] = lx

    pincy = np.zeros((4, ly.shape[0], ly.shape[1]))
    pincy[[0, 1], :, :] = ly
    pincy[[2, 3], :, :] = ly + 1

    pincw = np.zeros((4, shx.shape[0], shx.shape[1]))
    pincw[0, :, :] = (1 - shx) * (1 - shy)
    pincw[1, :, :] = shx * (1 - shy)
    pincw[2, :, :] = shx * shy
    pincw[3, :, :] = (1 - shx) * shy

    axy = ax**2 + ay**2
    axyinap = np.clip(axy, cobs**2. + 1.e-3, 0.999)
    #sizeaxyinap=axyinap.shape[1]# not used

    #pincw = pincw*axyinap[pincx+(pincy-1)*sizeaxyinap] --->

    for z in range(pincw.shape[0]):
        for i in range(pincw.shape[1]):
            for j in range(pincw.shape[2]):
                pincw[z, i, j] = pincw[z, i, j] * axyinap[
                        np.int32(pincx[z, i, j]),
                        np.int32(pincy[z, i, j])]

    pincw = pincw * np.tile(1.0 / np.sum(pincw, axis=0), (4, 1, 1))

    return pincx, pincy, pincw


#__________________________________________________________________________
#__________________________________________________________________________


def pcgeom(nr, npp, cobs, ncp, ncmar):
    #pcgeom,bas,ncp,ncmar;
    #DOCUMENT pcgeom,&geom,ncp,ncmar
    #This routine builds a geom_struct. px and py are the x and y
    #coordinates of points in the polar arrays.  cr and cp are the
    #r and phi coordinates of points in the cartesian grids. ncmar
    #allows the possibility that there is a margin of ncmar points
    #in the cartesian arrays outside the region of interest
    nused = ncp - 2 * ncmar
    ff = 0.5 * nused
    hw = np.float(ncp - 1) / 2.

    r = radii(nr, npp, cobs)
    p = polang(r)

    px0 = r * np.cos(p)
    py0 = r * np.sin(p)
    px = ff * px0 + hw
    py = ff * py0 + hw
    ax = np.reshape(
            np.arange(int(ncp)**2, dtype=np.float) + 1, (int(ncp), int(ncp)),
            order='F')
    ax = np.float32(ax - 1) % ncp - 0.5 * (ncp - 1)
    ax = ax / (0.5 * nused)
    ay = np.transpose(ax)

    pincx, pincy, pincw = setpincs(ax, ay, px0, py0, cobs)

    dpi = 2 * np.pi
    cr2 = (ax**2 + ay**2)
    ap = np.clip(cr2, cobs**2 + 1.e-3, 0.999)
    #cr = (cr2 - cobs**2) / (1 - cobs**2) * nr - 0.5;
    cr = (cr2 - cobs**2) / (1 - cobs**2) * nr
    cp = (np.arctan2(ay, ax) + dpi) % dpi
    cp = (npp / dpi) * cp

    cr = np.clip(cr, 1.e-3, nr - 1.001)
    #fudge -----, but one of the less bad ones
    cp = np.clip(cp, 1.e-3, npp - 1.001)
    #fudge -----  this is the line which
    #gives that step in the cartesian grid
    #at phi = 0.
    return ncp, ncmar, px, py, cr, cp, pincx, pincy, pincw, ap


#__________________________________________________________________________
#__________________________________________________________________________


def set_pctr(dim, nr, npp, nkl, cobs, nord, ncmar=None, ncp=None):
    #set_pctr,klbasis, ncp= dim
    #DOCUMENT geom=set_pctr(bas, ncp =, ncmar=)
    #This routine calls pcgeom to build a geom_struct with the
    #right initializations. bas is a gkl_basis_struct built with
    #the gkl_bas routine.
    ncp = dim
    if (ncmar == None):
        ncmar = 2
    if (ncp == None):
        ncp = 128
    ncp, ncmar, px, py, cr, cp, pincx, pincy, pincw, ap = pcgeom(
            nr, npp, cobs, ncp, ncmar)
    return ncp, ncmar, px, py, cr, cp, pincx, pincy, pincw, ap


#__________________________________________________________________________

#__________________________________________________________________________


def gkl_fcom(kers, cobs, nf):
    #gkl_fcom,kers,cobs,nkl,evals,nord,npo,ordd,rabas;
    #DOCUMENT gkl_fcom(kers,cobs,nf,&evals,&nord,&npo,&ordd,&rabas)
    #This routine does the work : finding the eigenvalues and
    #corresponding eigenvectors. Sort them and select the right
    #one. It returns the KL modes : in polar coordinates : rabas
    #as well as the associated variance : evals. It also returns
    #a bunch of indices used to recover the modes in cartesian
    #coordinates (nord, npo and ordd).
    nkl = nf
    st = kers.shape
    print(st)
    nr = st[1]
    nt = st[0]
    nxt = 0
    fktom = (1. - cobs**2) / nr
    #fevtos = np.sqrt(2*nr) #not used

    evs = np.zeros((nr, nt), dtype=np.float32)
    #ff isnt used - the normalisation for
    #the eigenvectors is straightforward:
    #integral of surface**2 divided by area = 1,
    #and the cos**2 term gives a factor
    #half, so multiply zero order by
    #sqrt(n) and the rest by sqrt (2n)

    #zero order is a special case...
    #need to deflate to eliminate infinite eigenvalue - actually want
    #evals/evecs of zom - b where b is big and negative
    zom = kers[0, :, :]
    s = piston_orth(nr)

    ts = np.transpose(s)
    #b1 = ((ts(,+)*zom(+,))(,+)*s(+,))(1:nr-1, 1:nr-1)
    btemp = (ts.dot(zom).dot(s))[0:nr - 1, 0:nr - 1]

    #newev = SVdec(fktom*b1,v0,vt)
    v0, newev, vt = np.linalg.svd(fktom * btemp, full_matrices=True)

    v1 = np.zeros((nr, nr), dtype=np.float32)
    v1[0:nr - 1, 0:nr - 1] = v0
    v1[nr - 1, nr - 1] = 1

    vs = s.dot(v1)
    newev = np.concatenate((newev, [0]))
    #print(np.size(newev))
    evs[:, nxt] = np.float32(newev)
    kers[nxt, :, :] = np.sqrt(nr) * vs

    nxt = 1
    while True:
        vs, newev, vt = np.linalg.svd(
                fktom * kers[nxt, :, :], full_matrices=True)
        #newev = SVdec(fktom*kers(,,nxt),vs,vt)
        evs[:, nxt] = np.float32(newev)
        kers[nxt, :, :] = np.sqrt(2. * nr) * vs
        mxn = max(np.float32(newev))
        egtmxn = np.floor(evs[:, 0:nxt + 1] > mxn)
        nxt = nxt + 1
        if ((2 * np.sum(egtmxn) - np.sum(egtmxn[:, 0])) >= nkl):
            break

    nus = nxt - 1
    kers = kers[0:nus + 1, :, :]

    #evs = reform (evs [:, 1:nus], nr*(nus))

    evs = np.reshape(evs[:, 0:nus + 1], nr * (nus + 1), order='F')
    a = np.argsort(-1. * evs)[0:nkl]

    #every eigenvalue occurs twice except
    #those for the zeroth order mode. This
    #could be done without the loops, but
    #it isn't the stricking point anyway...

    no = 0
    ni = 0
    #oind = array(long,nf+1)
    oind = np.zeros(nkl + 1, dtype=np.int32)

    while True:
        if (a[ni] < nr):
            oind[no] = a[ni]
            no = no + 1
        else:
            oind[no] = a[ni]
            oind[no + 1] = a[ni]
            no = no + 2

        ni = ni + 1
        if (no >= (nkl)):
            break

    oind = oind[0:nkl]
    tord = (oind) / nr + 1

    odd = np.arange(nkl, dtype=np.int32) % 2
    pio = (oind) % nr + 1

    evals = evs[oind]
    ordd = 2 * (tord - 1) - np.floor((tord > 1) & (odd)) + 1

    nord = max(ordd)

    rabas = np.zeros((nr, nkl), dtype=np.float32)
    sizenpo = np.int32(nord)
    npo = np.zeros(sizenpo, dtype=np.int32)

    for i in range(nkl):
        npo[np.int32(ordd[i]) - 1] = npo[np.int32(ordd[i]) - 1] + 1
        rabas[:, i] = kers[tord[i] - 1, :, pio[i] - 1]

    return evals, nord, npo, ordd, rabas


#__________________________________________________________________________


def gkl_sfi(_klbas, i):
    #DOCUMENT
    #This routine returns the i'th function from the generalised KL
    #basis bas. bas must be generated first with gkl_bas.
    nr = _klbas.nr
    npp = _klbas.npp
    ordp = _klbas.ordd
    rabasp = _klbas.rabas
    azbasp = _klbas.azbas
    nfunc = _klbas.nfunc

    if (i > nfunc - 1):
        raise TypeError("kl funct order it's so big")

    else:

        ordi = np.int32(ordp[i])
        rabasi = rabasp[:, i]

        azbasp = np.transpose(azbasp)
        azbasi = azbasp[ordi - 1, :]

        sf1 = np.zeros((nr, npp), dtype=np.float64)
        for j in range(npp):
            sf1[:, j] = rabasi

        sf2 = np.zeros((npp, nr), dtype=np.float64)
        for j in range(nr):
            sf2[:, j] = azbasi

        sf = sf1 * np.transpose(sf2)

        return sf


def pol2car(pol, _klbas, mask=1):
    # DOCUMENT cart=pol2car(cpgeom, pol, mask=)
    # This routine is used for polar to cartesian conversion.
    # pol is built with gkl_bas and cpgeom with pcgeom.
    # However, points not in the aperture are actually treated
    # as though they were at the first or last radial polar value
    # -- a small fudge, but not serious  ?*******
    #cd = interpolate.interp2d(cr, cp,pol)
    #    ncp = _klbas.ncp
    cr = _klbas.cr
    cp = _klbas.cp
    nr = _klbas.nr
    npp = _klbas.npp

    r = np.arange(nr, dtype=np.float64)
    phi = np.arange(npp, dtype=np.float64)
    tab_phi, tab_r = np.meshgrid(phi, r)

    cdc = interpolate.griddata((tab_r.flatten(), tab_phi.flatten()),
                               pol.flatten(), (cr, cp),
                               method='cubic')
    cdl = interpolate.griddata((tab_r.flatten(), tab_phi.flatten()),
                               pol.flatten(), (cr, cp),
                               method='linear')

    if (mask == 1):
        ap = _klbas.ap
        cdc = cdc * (ap)
        cdl = cdl * (ap)
        #cdxy = cdxy*(ap)

    return cdc, cdl  #,cdxy


def kl_view(_klbas, mask=1):

    nfunc = _klbas.nfunc
    ncp = _klbas.ncp

    tab_kl_c = np.zeros((nfunc, ncp, ncp), dtype=np.float64)
    tab_kl_l = np.zeros((nfunc, ncp, ncp), dtype=np.float64)

    for i in range(nfunc):

        tab_kl_c[i, :, :], tab_kl_l[i, :, :] = pol2car(
                gkl_sfi(_klbas, i), _klbas, mask)

    return tab_kl_c, tab_kl_l  #,tab_klxy


def kl_view_onepic(_klbas, imax, mask=1):

    tab_kl_c, tab_kl_l = kl_view(_klbas, mask)

    imax = (np.floor(np.sqrt(imax)))**2

    print("number of kl view :", imax)

    part_c = tab_kl_c[0:imax, :, :]
    part_l = tab_kl_l[0:imax, :, :]

    part_r_c = np.zeros((
            np.sqrt(part_c.shape[0]) * part_c.shape[1],
            np.sqrt(part_c.shape[0]) * part_c.shape[2]))
    part_r_l = np.zeros((
            np.sqrt(part_l.shape[0]) * part_l.shape[1],
            np.sqrt(part_l.shape[0]) * part_l.shape[2]))

    for i in range(np.int32(np.sqrt(imax))):
        for j in range(np.int32(np.sqrt(imax))):

            part_r_c[i * part_c.shape[1]:(i + 1) * part_c.shape[1],
                     j * part_c.shape[1]:(j + 1) * part_c.shape[1]] = part_c[
                             i * np.sqrt(imax) + j, :, :]
            part_r_l[i * part_l.shape[1]:(i + 1) * part_l.shape[1],
                     j * part_l.shape[1]:(j + 1) * part_l.shape[1]] = part_l[
                             i * np.sqrt(imax) + j, :, :]

    return part_r_c, part_r_l
