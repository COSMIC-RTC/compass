from cython.operator cimport dereference as deref, preincrement as inc

import make_pupil as mkP
from blaze.expr.expressions import shape

import numpy as np
cimport numpy as np
np.import_array()

import hdf5_utils as h5
import resDataBase as db
import pandas as pd
from scipy import interpolate
import copy as copy
from scipy.sparse import csr_matrix
import shesha_kl as klfunc

#max_extent signature
cdef _dm_init(Dms dms, Param_dm p_dms, list p_wfs, Param_geom p_geom, Param_tel p_tel,int *max_extent):
    """ inits a Dms object on the gpu

    :parameters:
        dms: (Dms) :

        p_dms: (Param_dms) : dm settings

        p_wfs: (list) : list of wfs settings

        p_geom: (Param_geom) : geom settings

        p_tel: (Param_tel) : telescope settings

        max_extend: (int*) :
    """

    cdef float patchDiam
    cdef int extent
    cdef long dim
    cdef long ninflu, influsize, ninflupos, n_npts
    cdef long _nr, _np

    cdef float tmp

    if(p_dms.pupoffset is not None):
        p_dms.puppixoffset = p_dms.pupoffset / p_tel.diam * p_geom.pupdiam


    cdef float irc, coupling, ir
    cdef long pitch, smallsize

    if( p_dms.type_dm=="pzt"):
        if p_dms.file_influ_hdf5 == None:
            # calcul pitch ______________________

            # find out the support dimension for the given mirror.
            norms = [np.linalg.norm([w.xpos, w.ypos]) for w in p_wfs]
            patchDiam = p_geom.pupdiam + 2 * np.max(norms) * \
                4.848e-6 * np.abs(p_dms.alt) / p_tel.diam * p_geom.pupdiam
            # Patchdiam
            p_dms._pitch = long(patchDiam / (p_dms.nact - 1))

            extent = p_dms._pitch * (p_dms.nact + 5)  # + 2.5 pitch each side
            p_dms._n1 = np.floor(p_geom.cent - extent / 2)
            p_dms._n2 = np.ceil(p_geom.cent + extent / 2)
            if(p_dms._n1 < 1):
                p_dms._n1 = 1
            if(p_dms._n2 > p_geom.ssize):
                p_dms._n2 = p_geom.ssize

            #calcul influsize ___________________

            coupling=p_dms.coupling

            irc=1.16136+2.97422*coupling+(-13.2381)*coupling**2+20.4395*coupling**3
            pitch=p_dms._pitch
            ir=irc*pitch

            smallsize=np.ceil(2*ir+10)
            if(smallsize%2!=0):smallsize+=1
            p_dms._influsize=smallsize
            #______________________________________

            #calcul defaut influsize
            make_pzt_dm(p_dms,p_geom,p_tel,irc)
        else :
            read_influ_hdf5 (p_dms,p_tel,p_geom)

        # max_extent
        max_extent[0] = max(max_extent[0], p_dms._n2 - p_dms._n1 + 1)

        dim = max(p_dms._n2-p_dms._n1+1, p_geom._mpupil.shape[0])
        ninflu=long(p_dms._ntotact)
        influsize=long(p_dms._influsize)
        ninflupos=long(p_dms._influpos.size)
        n_npts = long(p_dms._ninflu.size)

        dms.add_dm(p_dms.type_dm, p_dms.alt, dim, ninflu, influsize,
                   ninflupos, n_npts, p_dms.push4imat)
        dms.load_pzt(p_dms.alt, p_dms._influ, p_dms._influpos.astype(np.int32),
                     p_dms._ninflu, p_dms._influstart, p_dms._i1, p_dms._j1, None)

    elif(p_dms.type_dm == "tt"):


        if(p_dms.alt == 0):
            # TODO check next line (different max used ?)
            # (max among all the p_dms, previous p_dms._n modified)
            # adapted: #max_extent
            # max_extent=max(p_dms._n2-p_dms._n1+1)
            extent = int(max_extent[0] * 1.05)
            if (extent % 2 != 0):
                extent += 1
            p_dms._n1 = np.floor(p_geom.cent - extent / 2)
            p_dms._n2 = np.ceil(p_geom.cent + extent / 2)
            if(p_dms._n1 < 1):
                p_dms._n1 = 1
            if(p_dms._n2 > p_geom.ssize):
                p_dms._n2 = p_geom.ssize
        else :
            extent = p_geom.pupdiam + 16
            p_dms._n1 = np.floor(p_geom.cent - extent / 2)
            p_dms._n2 = np.ceil(p_geom.cent + extent / 2)
            if(p_dms._n1 < 1):
                p_dms._n1 = 1
            if(p_dms._n2 > p_geom.ssize):
                p_dms._n2 = p_geom.ssize
            # max_extent
            max_extent[0] = max(max_extent[0], p_dms._n2 - p_dms._n1 + 1)




        dim = long(p_dms._n2 - p_dms._n1 + 1)
        make_tiptilt_dm(p_dms, p_wfs, p_geom, p_tel)
        dms.add_dm(p_dms.type_dm, p_dms.alt, dim, 2, dim,
                   1, 1, p_dms.push4imat)
        dms.load_tt(p_dms.alt, p_dms._influ)

    elif(p_dms.type_dm == "kl"):

        extent = p_geom.pupdiam + 16
        p_dms._n1 = np.floor(p_geom.cent - extent / 2)
        p_dms._n2 = np.ceil(p_geom.cent + extent / 2)
        if(p_dms._n1 < 1):
            p_dms._n1 = 1
        if(p_dms._n2 > p_geom.ssize):
            p_dms._n2 = p_geom.ssize
        # max_extent
        max_extent[0] = max(max_extent[0], p_dms._n2 - p_dms._n1 + 1)


        dim = long(p_dms._n2 - p_dms._n1 + 1)
        
        make_kl_dm(p_dms,p_wfs,p_geom,p_tel)
            
        ninflu = p_dms.nkl
        influsize = long(p_dms._klbas.ncp)
        _nr = long(p_dms._klbas.nr)
        _npp = long(p_dms._klbas.npp)
        ord_L = copy.copy(p_dms._klbas.ordd)
        rabas_L = copy.copy(p_dms._klbas.rabas.flatten('F'))
        azbas_L = copy.copy(p_dms._klbas.azbas.flatten('F'))
        cr_L = copy.copy(p_dms._klbas.cr.flatten('F'))
        cp_L = copy.copy(p_dms._klbas.cp.flatten('F'))
        
        dms.add_dm(p_dms.type_dm, p_dms.alt, dim, ninflu, influsize,
                   _nr, _npp, p_dms.push4imat)                           
        dms.load_kl(p_dms.alt, np.float32(rabas_L), np.float32(azbas_L),
                    np.int32(ord_L), np.float32(cr_L), np.float32(cp_L))


    else :

        raise StandardError("This type of DM doesn't exist ")
        # Verif
        # res1 = pol2car(*y_dm(n)._klbas,gkl_sfi(*y_dm(n)._klbas, 1));
        # res2 = yoga_getkl(g_dm,0.,1);


def dm_init(p_dms, list p_wfs, Sensors sensors, Param_geom p_geom, Param_tel p_tel):
    """Create and initialize a Dms object on the gpu

    :parameters:
        p_dms: (list of Param_dms) : dms settings

        p_wfs: (Param_wfs) : wfs settings

        p_geom: (Param_geom) : geom settings

        p_tel: (Param_tel) : telescope settings
    """
    cdef int max_extent = 0
    if(len(p_dms) != 0):
        dms = Dms(len(p_dms))
        for i in range(len(p_dms)):
             #max_extent
            _dm_init(dms, p_dms[i], p_wfs, p_geom, p_tel, & max_extent)

    if(p_wfs is not None):
        if(sensors is not None):
            for i in range(len(p_wfs)):
                if(not p_wfs[i].openloop):
                    for j in range(p_wfs[i].dms_seen.size):
                        k = p_wfs[i].dms_seen[j]
                        dims = p_dms[k]._n2 - p_dms[k]._n1 + 1
                        dim = p_geom._mpupil.shape[0]
                        if(dim < dims):
                            dim = dims
                        xoff = p_wfs[i].xpos * 4.848e-6 * p_dms[k].alt / p_tel.diam * p_geom.pupdiam
                        yoff = p_wfs[i].ypos * 4.848e-6 * p_dms[k].alt / p_tel.diam * p_geom.pupdiam
                        xoff = xoff + (dim - p_geom._n) / 2
                        yoff = yoff + (dim - p_geom._n) / 2
                        sensors.sensors.d_wfs[i].d_gs.add_layer(p_dms[k].type_dm, p_dms[k].alt, xoff, yoff)
    return dms


cpdef createSquarePattern(float pitch, int nxact ):
    """
    Creates a list of M=nxact^2 actuator positions spread over an square grid.
    Coordinates are centred around (0,0).

    :parameters:
        pitch: (float) : distance in pixels between 2 adjacent actus
        nxact: (int) : number of actu across the pupil diameter
    :return:
        xy: (np.ndarray(dims=2,dtype=np.float32)) : xy[M,2] list of coodinates
    """

    xy = np.tile( np.arange(nxact) - (nxact-1.)/2. , (nxact,1)).astype(np.float32)
    xy = np.array([xy.flatten(), xy.T.flatten()]) * pitch
    xy = np.float32(xy)
    return xy


cpdef createHexaPattern(float pitch, float supportSize):
    """
    Creates a list of M actuator positions spread over an hexagonal grid.
    The number M is the number of points of this grid, it cannot be
    known before the procedure is called.
    Coordinates are centred around (0,0).
    The support that limits the grid is a square [-n/2,n/2].

    :parameters:
        pitch: (float) : distance in pixels between 2 adjacent actus
        n: (float) : size in pixels of the support over which the coordinate list
             should be returned.
    :return:
        xy: (np.ndarray(dims=2,dtype=np.float32)) : xy[M,2] list of coodinates
    """
    V3 = np.sqrt(3)
    nx = int(np.ceil((supportSize/2.0)/pitch) + 1)
    x = pitch * (np.arange(2*nx+1, dtype=np.float32)-nx)
    Nx = x.shape[0]
    ny = int(np.ceil((supportSize/2.0)/pitch/V3) + 1)
    y = (V3 * pitch) * (np.arange(2*ny+1, dtype=np.float32)-ny)
    Ny = y.shape[0]
    x = np.tile(x,(Ny,1)).flatten()
    y = np.tile(y,(Nx,1)).T.flatten()
    x = np.append(x, x + pitch/2.)
    y = np.append(y, y + pitch*V3/2.)
    xy = np.float32(np.array([x,y]))
    return xy


def n_actuator_select(Param_dm p_dm,Param_tel p_tel, xc,yc):
    """
    Fonction for select actuator in fonction of Margin_in, margin_out or ntotact.
    default margin_out=1.44pitch, default for margin_in taking all the actuators.


    :parameters:
        p_dm: (Param_dm) : dm settings
        xc: actuators x positions (origine in center of mirror)
        yc: actuators y positions (origine in center of mirror)

    :return:
        liste_fin: actuator indice selection for xpos/ypos


    """
    # the following determine if an actuator is to be considered or not
    # relative to the pitchmargin parameter.
    dis=np.sqrt(xc**2+yc**2)
    #cdef np.ndarray dis2=np.sqrt(cub[0,:]**2+cub[1,:]**2)
    #cdef float rad_in, rad_out

    #test Margin_in

    if(p_dm.margin_in==0):
        # 1 if valid actuator, 0 if not:

        rad_in=0.0

    else:
        pitchMargin_in=p_dm.margin_in
        rad_in=(((p_dm.nact-1.)/2.)*p_tel.cobs-pitchMargin_in)*p_dm._pitch

    if(p_dm._ntotact==0):
        if(p_dm.margin_out==0):
            pitchMargin_out=1.44
        else:
            pitchMargin_out=p_dm.margin_out
        rad_out=((p_dm.nact-1.)/2.+pitchMargin_out)*p_dm._pitch

        liste_fin = np.where((dis < rad_out) * (dis > rad_in))[0]

    else:
        liste_i = sorted(range(len(dis)), key=lambda k: dis[k])
        liste2 = dis[liste_i] > rad_in


        if(sum(liste2)<p_dm._ntotact):
            print 'ntoact very hight'
            liste_fin = liste_i[(np.size(liste2)-sum(liste2)):]
        else:
            liste_fin = liste_i[(np.size(liste2)-sum(liste2)):p_dm._ntotact+(np.size(liste2)-sum(liste2))]


    return liste_fin


cpdef make_pzt_dm(Param_dm p_dm,Param_geom geom,Param_tel p_tel,irc):

    """Compute the actuators positions and the influence functions for a pzt DM

    :parameters:
        p_dm: (Param_dm) : dm settings

        geom: (Param_geom) : geom settings

        p_tel: (Param_tel) : tel settings

        irc: factor for influence size
    :return:
        influ: (np.ndarray(dims=3,dtype=np.float64)) : cube of the IF for each actuator

    """
    cdef int i
    #best parameters, as determined by a multi-dimensional fit
    #(see coupling3.i)
    cdef float p1,p2, coupling
    coupling=p_dm.coupling

    p1=4.49469+(7.25509+(-32.1948+17.9493*coupling)*coupling)*coupling
    p2=2.49456+(-0.65952+(8.78886-6.23701*coupling)*coupling)*coupling

    cdef float tmp_c=1.0/np.abs(irc)
    cdef float ccc = (coupling - 1.+ tmp_c**p1)/(np.log(tmp_c)*tmp_c**p2)

    # prepare to compute IF on partial (local) support of size <smallsize>
    cdef long pitch=p_dm._pitch
    cdef long smallsize

    smallsize=p_dm._influsize

    # compute location (x,y and i,j) of each actuator:
    cdef long nxact=p_dm.nact
    cdef np.ndarray cub

    if p_dm.type_pattern == None:
        p_dm.type_pattern = <bytes>'square'

    if p_dm.type_pattern == 'hexa':
        print "Patter type : Hexa"
        cub = createHexaPattern( pitch, geom.pupdiam * 1.1)
    elif p_dm.type_pattern == 'square':
        print "Patter default type : Square"
        cub = createSquarePattern( pitch, nxact )
    else :
        raise StandardError("This pattern does not exist for pzt dm")

    inbigcirc = n_actuator_select(p_dm,p_tel,cub[0,:],cub[1,:])

    #print 'inbigcirc',inbigcirc.shape

    # converting to array coordinates:
    cub += geom.cent

    # filtering actuators outside of a disk radius = rad (see above)
    cdef np.ndarray cubval = cub[:,inbigcirc]
    ntotact = cubval.shape[1]
    print "ntact =",ntotact
    xpos    = cubval[0,:]
    ypos    = cubval[1,:]
    i1t      = (cubval[0,:]-smallsize/2+0.5-p_dm._n1).astype(np.int32)
    j1t      = (cubval[1,:]-smallsize/2+0.5-p_dm._n1).astype(np.int32)

    # Allocate array of influence functions
    cdef np.ndarray[ndim=3,dtype=np.float32_t] influ=np.zeros((smallsize,smallsize,ntotact),dtype=np.float32)
    # Computation of influence function for each actuator
    cdef int i1, j1
    cdef np.ndarray[ndim=2,dtype=np.float32_t] x, y, tmp
    for i in range(ntotact):
        i1 = i1t[i]
        x  = np.tile(np.arange(i1,i1+smallsize,dtype=np.float32),(smallsize,1)) # pixel coords in ref frame "dm support"
        x += p_dm._n1                                          # pixel coords in ref frame "pupil support"
        x -= xpos[i]                                     # pixel coords in local ref frame
        x  = np.abs(x)/(irc*pitch)                             # normalized coordiantes in local ref frame

        j1 = j1t[i]
        y  = np.tile(np.arange(j1,j1+smallsize,dtype=np.float32),(smallsize,1)) # idem as X, in Y
        y += p_dm._n1
        y -= ypos[i]
        y  = np.abs(y.T)/(irc*pitch)

        #clip
        x[x<1e-8]=1e-8
        x[x>2]=2.
        y[y<1e-8]=1e-8
        y[y>2]=2.
        tmp = (1.-x**p1+ccc*np.log(x)*x**p2)*(1.-y**p1+ccc*np.log(y)*y**p2)
        tmp = tmp*(x <= 1.0)*(y <= 1.0)
        influ[:, :, i] = tmp

    if(p_dm._puppixoffset is not None):
        xpos +=p_dm._puppixoffset[0]
        ypos +=p_dm._puppixoffset[1]

    influ=influ*float(p_dm.unitpervolt/np.max(influ))

    p_dm._influ = influ

    print  'number of actuator after filter = ',np.size(inbigcirc)
    p_dm._ntotact = np.size(inbigcirc)
    p_dm._xpos = xpos
    p_dm._ypos = ypos

    # i1, j1 calc :

    p_dm._i1 = (p_dm._xpos - p_dm._influsize/2. +0.5 - p_dm._n1).astype(np.int32)
    p_dm._j1 = (p_dm._ypos - p_dm._influsize/2. +0.5 - p_dm._n1).astype(np.int32)

    comp_dmgeom(p_dm,geom)

    cdef long dim=max(geom._mpupil.shape[0],p_dm._n2-p_dm._n1+1)
    cdef float off=(dim-p_dm._influsize)/2

    cdef np.ndarray[ndim=2, dtype=np.float32_t] kernconv=np.zeros((dim,dim),dtype=np.float32)
    kernconv [off:off+p_dm._influsize,off:off+p_dm._influsize] = p_dm._influ[:,:,0]
    kernconv=np.roll(kernconv,kernconv.shape[0]/2,axis=0)
    kernconv=np.roll(kernconv,kernconv.shape[1]/2,axis=1)
    p_dm._influkernel= kernconv


cpdef read_influ_hdf5 (Param_dm p_dm,Param_tel p_tel, Param_geom geom):
    """Read HDF for influence pzt fonction and form

    :parameters:
        p_dm: (Param_dm) : dm settings

        geom: (Param_geom) : geom settings

        p_tel: (Param_tel) : tel settings

    """
    # read h5 file for influence fonction
    h5_tp = pd.read_hdf(p_dm.file_influ_hdf5,'resAll')
    print "Read Ifluence fonction in h5 : ",p_dm.file_influ_hdf5

    # cube_name
    influ_h5 = h5_tp[p_dm.cube_name][0]


    # x_name
    xpos_h5 = h5_tp[p_dm.x_name][0]


    # y_name
    ypos_h5 = h5_tp[p_dm.y_name][0]


    # center_name
    center_h5 = h5_tp[p_dm.center_name][0]


    # influ_res
    res_h5 = h5_tp[p_dm.influ_res][0]
    res_h5_m = (res_h5[0]+res_h5[1])/2.


    # a introduire
    diam_h5 = [2.54,2.54] # metre

    #soustraction du centre introduit
    xpos_h5_0 = xpos_h5-center_h5[0]
    ypos_h5_0 = ypos_h5-center_h5[1]

    # interpolation du centre (ajout du nouveau centre)
    center = geom.cent


    # calcul de la resolution de la pupille
    res_compass = p_tel.diam/geom.pupdiam


    # interpolation des coordonnées en pixel avec ajout du centre
    xpos = (xpos_h5_0*(43.68/diam_h5[0]))/res_compass + center
    ypos = (ypos_h5_0*(43.68/diam_h5[1]))/res_compass + center

    # interpolation des fonction d'influence

    influ_size_h5 = influ_h5.shape[0]
    ninflu = influ_h5.shape[2]


    x = np.arange(influ_size_h5)*res_h5_m*(p_tel.diam/diam_h5[0])
    y = np.arange(influ_size_h5)*res_h5_m*(p_tel.diam/diam_h5[1])
    xmax = max(x)
    ymax = max(y)
    xnew = np.arange(0,xmax,res_compass)
    xnew = xnew+(xmax-max(xnew))/2.
    ynew = np.arange(0,ymax,res_compass)
    ynew = ynew+(ymax-max(ynew))/2.
    influ_size = xnew.shape[0]

    #creation du ouveaux cube d'influance
    influ_new = np.zeros((influ_size,influ_size,ninflu))


    for i in range(ninflu):

        influ = influ_h5[:,:,i]
        f = interpolate.interp2d(x,y,influ,kind='cubic')
        influ_new[:,:,i] = f(xnew, ynew)


    p_dm._xpos = np.float32(xpos)
    p_dm._ypos = np.float32(ypos)

    # number of actuator
    print "Actuator number in H5 data : ",ninflu
    p_dm._ntotact = np.int(ninflu)

    # def influente fonction normalize by unitpervolt
    p_dm._influ = np.float32(influ_new)*p_dm.unitpervolt

    # def influence size
    print "influence size in pupil : ",np.int(influ_size),"pixel"
    p_dm._influsize = np.int(influ_size)


    # Def dm limite (n1 and n2)
    extent = (max(xpos) - min(xpos))+(influ_size*2)
    p_dm._n1 = np.floor(geom.cent - extent / 2)
    p_dm._n2 = np.ceil(geom.cent + extent / 2)
    if(p_dm._n1 < 1):
        p_dm._n1 = 1
    if(p_dm._n2 > geom.ssize):
        p_dm._n2 = geom.ssize

    # refaire la definition du pitch pour n_actuator
    #inbigcirc = n_actuator_select(p_dm,p_tel,xpos-center[0],ypos-center[1])
    #print 'nb = ',np.size(inbigcirc)
    #p_dm._ntotact = np.size(inbigcirc)


    # i1, j1 calc :

    p_dm._i1 = (p_dm._xpos - p_dm._influsize/2. +0.5 - p_dm._n1).astype(np.int32)
    p_dm._j1 = (p_dm._ypos - p_dm._influsize/2. +0.5 - p_dm._n1).astype(np.int32)


    comp_dmgeom(p_dm,geom)

    cdef long dim=max(geom._mpupil.shape[0],p_dm._n2-p_dm._n1+1)
    cdef float off=(dim-p_dm._influsize)/2

    cdef np.ndarray[ndim=2, dtype=np.float32_t] kernconv=np.zeros((dim,dim),dtype=np.float32)
    kernconv [off:off+p_dm._influsize,off:off+p_dm._influsize] = p_dm._influ[:,:,0]
    kernconv=np.roll(kernconv,kernconv.shape[0]/2,axis=0)
    kernconv=np.roll(kernconv,kernconv.shape[1]/2,axis=1)
    p_dm._influkernel= kernconv


cpdef make_tiptilt_dm(Param_dm p_dm,list p_wfs, Param_geom p_geom, Param_tel p_tel):
    """Compute the influence functions for a tip-tilt DM

    :parameters:
        p_dm: (Param_dm) : dm settings

        p_wfs: (list) : list of wfs settings

        p_geom: (Param_geom) : geom settings

        p_tel: (Param_tel) : telescope settings
    :return:
        influ: (np.ndarray(dims=3,dtype=np.float64)) : cube of the IF

    """
    cdef int dim = max(p_dm._n2 - p_dm._n1 + 1, p_geom._mpupil.shape[0])
    norms = [np.linalg.norm([w.xpos, w.ypos]) for w in p_wfs]
    cdef int patchDiam = p_geom.pupdiam + \
        2 * np.max(norms) * 4.848e-6 * \
        abs(p_dm.alt / p_tel.diam * p_geom.pupdiam)

    cdef nzer = 2
    cdef np.ndarray[ndim = 3, dtype = np.float32_t] influ = make_zernike(nzer + 1, dim,
                                                                         patchDiam, p_geom.cent - p_dm._n1 + 1, p_geom.cent - p_dm._n1 + 1, 1)[:, :, 1:]

    # normalization factor: one unit of tilt gives 1 arcsec:
    cdef float current = influ[dim / 2 - 1, dim / 2 - 1, 0] - influ[dim / 2 - 2, dim / 2 - 2, 0]
    cdef float fact = p_dm.unitpervolt * p_tel.diam / p_geom.pupdiam * 4.848 / current

    influ = influ * fact
    p_dm._ntotact = influ.shape[2]

    cdef int i
    p_dm._influ = influ

    return influ

cpdef make_klbas(Param_dm p_dm, int nkl,float cobs, long dim,funct,float outscl=3):
    #DOCUMENT make_klbas(nfunc,cobs,nr=,np=,funct=,outscl=)

    print(funct)
    nr = np.long(5.0*np.sqrt(nkl))
    npp = np.long(5*nr)      
    radp = klfunc.make_radii(cobs,nr)

    kers = klfunc.make_kernels(cobs,nr,radp,funct,outscl)

    evals,nord,npo,ordd,rabas = klfunc.gkl_fcom(kers,cobs,nkl)

    azbas = klfunc.make_azimuth(nord,npp)

    ncp,ncmar,px,py,cr,cp,pincx,pincy,pincw,ap = klfunc.set_pctr(dim,nr,npp,nkl,cobs,nord)

    #make klbas
    p_dm._klbas.nr = nr # number of radial points
    p_dm._klbas.npp = npp # number of elements
    p_dm._klbas.nfunc = nkl# number of KL
    p_dm._klbas.nord = nord
    p_dm._klbas.radp = radp
    p_dm._klbas.evals = evals # veriance of kl number
    p_dm._klbas.npo = npo
    p_dm._klbas.ordd = ordd  #the radial orders of the basis
    p_dm._klbas.rabas = rabas # the radial array of the basis
    p_dm._klbas.azbas = np.transpose(azbas) #the azimuthal array of the basis
    p_dm._klbas.kers = kers
    

    #pcgeom
    p_dm._klbas.ncp = ncp # dim of grid
    p_dm._klbas.ncmar = ncmar # marge
    p_dm._klbas.px = px # x coord in polar array
    p_dm._klbas.py = py # y coord in polar array
    p_dm._klbas.cr = cr # radial coord in cartesien grid
    p_dm._klbas.cp = cp # phi coord in cartesien grid
    p_dm._klbas.pincx = pincx
    p_dm._klbas.pincy = pincy
    p_dm._klbas.pincw = pincw
    p_dm._klbas.ap = ap
    

    
        

cpdef make_kl_dm(Param_dm p_dm, list p_wfs,Param_geom p_geom, Param_tel p_tel):
    """Compute the influence function for a Karhunen-Loeve DM

    :parameters:
        p_dm: (Param_dm) : dm settings

        p_wfs: (list) : wfs settings

        p_geom: (Param_geom) : geom settings

        p_tel: (Param_tel) : telescope settings

    """
    cdef int dim = p_geom._mpupil.shape[0]
    norms = [np.linalg.norm([w.xpos, w.ypos]) for w in p_wfs]
    cdef long patchDiam = long(p_geom.pupdiam + 2 * np.max(norms) * 4.848e-6 *
                               abs(p_dm.alt) / p_geom.pupdiam)
                               
    print "dimkl = %d" %patchDiam

    make_klbas(p_dm,p_dm.nkl,p_tel.cobs,patchDiam,p_dm.kl_type)
    p_dm._i1 = np.zeros((p_dm.nkl), dtype=np.int32) + \
        int((dim - patchDiam) / 2.)
    p_dm._j1 = np.zeros((p_dm.nkl), dtype=np.int32) + \
        int((dim - patchDiam) / 2.)
    p_dm._ntotact = p_dm.nkl


cpdef make_zernike(int nzer,int size,int diameter, float xc=-1, float yc=-1, int ext=0):
    """Compute the zernike modes

    :parameters:
        nzer: (int) : number of modes

        size: (int) : size of the screen

        diameter: (int) : pupil diameter

        xc: (float) : (optional) x-position of the center


        yc: (float) : (optional) y-position of the center

        ext: (int) : (optional) extension
    :return:
        z : (np.ndarray(ndims=3,dtype=np.float64)) : zernikes modes
    """
    cdef int zn, i
    cdef int m = 0
    cdef int n = 0

    if(xc == -1):
        xc = size / 2
    if(yc == -1):
        yc = size / 2

    cdef float radius = (diameter + 1.) / 2.
    cdef np.ndarray[ndim = 2, dtype = np.float32_t] zr = mkP.dist(size, xc, yc).astype(np.float32).T / radius
    cdef np.ndarray[ndim = 3, dtype = np.float32_t] zmask = np.zeros((zr.shape[0], zr.shape[1], nzer),
                                                                     dtype=np.float32)
    cdef np.ndarray[ndim = 3, dtype = np.float32_t] zmaskmod = np.zeros((zr.shape[0], zr.shape[1], nzer),
                                                                        dtype=np.float32)

    zmask[:, :, 0] = (zr <= 1).astype(np.float32)
    zmaskmod[:, :, 0] = (zr <= 1.2).astype(np.float32)

    for i in range(1, nzer):
        zmask[:, :, i] = zmask[:, :, 0]
        zmaskmod[:, :, i] = zmaskmod[:, :, 0]

    cdef np.ndarray[ndim = 2, dtype = np.float32_t] zrmod = zr * zmaskmod[:, :, 0]

    zr = zr * zmask[:, :, 0]

    cdef np.ndarray[ndim = 2, dtype = np.float32_t] x = np.tile(np.linspace(1, size, size).astype(np.float32),
                                                                (size, 1))
    cdef np.ndarray[ndim = 2, dtype = np.float32_t] zteta = np.arctan2(x - yc, x.T - xc).astype(np.float32)

    cdef np.ndarray[ndim = 3, dtype = np.float32_t] z = np.zeros((size, size, nzer), dtype=np.float32)

    for zn in range(nzer):
        zernumero(zn + 1, & n, & m)

        if ext:
            for i in range((n - m) / 2 + 1):
                z[:, :, zn] = z[:, :, zn] + (-1.) ** i * zrmod ** (n - 2. * i) * float(np.math.factorial(n - i)) / \
                    float(np.math.factorial(i) * np.math.factorial((n + m) / 2 - i) *
                          np.math.factorial((n - m) / 2 - i))
        else:
            for i in range((n - m) / 2 + 1):
                z[:, :, zn] = z[:, :, zn] + (-1.) ** i * zr ** (n - 2. * i) * float(np.math.factorial(n - i)) / \
                    float(np.math.factorial(i) * np.math.factorial((n + m) / 2 - i) *
                          np.math.factorial((n - m) / 2 - i))

        if((zn + 1) % 2 == 1):
            if(m == 0):
                z[:, :, zn] = z[:, :, zn] * np.sqrt(n + 1.)
            else:
                z[:, :, zn] = z[:, :, zn] * \
                    np.sqrt(2. * (n + 1)) * np.sin(m * zteta)
        else:
            if(m == 0):
                z[:, :, zn] = z[:, :, zn] * np.sqrt(n + 1.)
            else:
                z[:, :, zn] = z[:, :, zn] * \
                    np.sqrt(2. * (n + 1)) * np.cos(m * zteta)

    if(ext):
        return z * zmaskmod
    else:
        return z * zmask


cdef zernumero(int zn, int * rd, int * an):
    """
    Returns the radial degree and the azimuthal number of zernike
    number zn, according to Noll numbering (Noll, JOSA, 1976)

    :parameters:
        zn: (int) : zernike number

        rd: (int*) : radial degrees

        an: (int*) : azimuthal numbers

    """
    cdef int j, m, n
    j = 0
    for n in range(101):
        for m in range(n + 1):
            if((n - m) % 2 == 0):
                j = j + 1
                if(j == zn):
                    rd[0] = n
                    an[0] = m
                    return
                if(m != 0):
                    j = j + 1
                    if(j == zn):
                        rd[0] = n
                        an[0] = m
                        return


cpdef comp_dmgeom(Param_dm dm, Param_geom geom):
    """Compute the geometry of a DM : positions of actuators and influence functions

    :parameters:
        dm: (Param_dm) : dm settings

        geom: (Param_geom) : geom settings
    """
    cdef int smallsize = dm._influsize
    cdef int nact = dm._ntotact
    cdef int dims = int(dm._n2 - dm._n1 + 1)
    cdef int dim = geom._mpupil.shape[0]
    cdef int offs

    if(dims < dim):
        offs = (dim - dims) / 2
    else:
        offs = 0
        dim = dims

    cdef np.ndarray[ndim = 2, dtype = np.float32_t] mapactu = np.ones((dim, dim), dtype=np.float32)

    if(offs > 0):
        mapactu[:offs, :] = 0
        mapactu[:, :offs] = 0
        mapactu[mapactu.shape[0] - offs:, :] = 0
        mapactu[:, mapactu.shape[1] - offs:] = 0

    cdef np.ndarray[ndim = 2, dtype = np.int32_t] indgen = np.tile(np.arange(smallsize, dtype=np.int32), (smallsize, 1))

    cdef np.ndarray[ndim = 3, dtype = np.int32_t] tmpx = np.tile(indgen, (nact, 1, 1)).T
    cdef np.ndarray[ndim = 3, dtype = np.int32_t] tmpy = np.tile(indgen.T, (nact, 1, 1)).T

    tmpx += offs + np.tile(dm._i1, (smallsize, smallsize, 1))
    tmpy += offs + np.tile(dm._j1, (smallsize, smallsize, 1))

    cdef np.ndarray[ndim = 3, dtype = np.int32_t] tmp = tmpx + dim * tmpy

    # bug in limit of def zone -10 destoe influpos for all actuator
    tmp[tmpx < 0] = -10
    tmp[tmpy < 0] = -10
    tmp[tmpx > dims - 1] = -10
    tmp[tmpy > dims - 1] = -10
    cdef np.ndarray[ndim = 1, dtype = np.int32_t] itmps = np.argsort(tmp.flatten("F")).astype(np.int32)
    cdef np.ndarray[ndim = 1, dtype = np.int32_t] tmps = np.sort(tmp.flatten("F")).astype(np.int32)
    itmps = itmps[np.where(itmps > -1)]


    cdef np.ndarray[ndim = 1, dtype = np.int32_t] istart, npts
    istart = np.zeros((dim * dim), dtype=np.int32)
    npts = np.zeros((dim * dim), dtype=np.int32)

    cdef int cpt, ref
    cpt = 0
    ref = 0

    for i in range(dim * dim):
        if(offs != 0 and mapactu.item(i) == 0):
            npts[i] = 0
        else:
            #print "DEBUG : tmps[cpt], i tmps.size-1 = ",tmps[cpt], i, tmps.size-1
            while(tmps[cpt] == i and cpt < tmps.size - 1):
                cpt += 1
            #print "DEBUG : cpt, ref", cpt, ref
            npts[i] = cpt - ref
            istart[i] = ref
            ref = cpt

    dm._influpos = itmps[:np.sum(npts)].astype(np.int32)
    dm._ninflu = npts.astype(np.int32)
    dm._influstart = istart.astype(np.int32)

    dm._i1 += offs
    dm._j1 += offs


#################################################
# P-Class Dms
#################################################
cdef class Dms:
    def __cinit__(self, long ndm):
        self.dms = new sutra_dms(ndm)

    def __dealloc__(self):
        del self.dms

    cpdef add_dm(self, bytes type_dm, float alt, long dim, long ninflu, long influsize, long ninflupos, long npts, float push4imat, int device=-1):
        """Add a dm into a Dms object

        :parameters:
            type_dm: (str) : dm type to remove,

            alt: (float) : dm conjugaison altitude to remove,

            ninflu: (long) : ,

            influsize: (long) : ,

            ninflupos: (long) : ,

            npts: (long) : ,

            push4imat: (float) : ,

            device: (int) : device where the DM will be create (default=-1):

        """

        cdef carma_context * context = &carma_context.instance()

        if(device > -1):
            context.set_activeDevice(device, 1)
        else:
            device = context.get_activeDevice()

        self.dms.add_dm(context, type_dm, alt, dim, ninflu,
                        influsize, ninflupos, npts, push4imat, device)

    cpdef remove_dm(self, bytes type_dm, float alt):
        """Remove a dm from a Dms object

        :parameters:
            type_dm: (str) : dm type to remove

            alt: (float) : dm conjugaison altitude to remove
        """
        self.dms.remove_dm(type_dm, alt)

    def resetdm(self, bytes type_dm, float alt):
        """Reset the shape of the DM to 0

        :parameters:
            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude
        """
        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            raise StandardError("Error in reset dm ")

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        self.dms.d_dms[inddm].reset_shape()

    def oneactu(self, bytes type_dm, float alt, int nactu, float ampli):
        """Push on on the nactu actuator of the DM with ampli amplitude and compute
        the corresponding shape

        :parameters:
            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude

            nactu: (int) : actuator number

            ampli: (float): amplitude
        """
        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            raise StandardError("One actuator error")

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        self.dms.d_dms[inddm].comp_oneactu(nactu, ampli)

    cpdef load_pzt(self, float alt,
                   np.ndarray[ndim=3, dtype=np.float32_t] influ,
                   np.ndarray[ndim=1, dtype=np.int32_t] influpos,
                   np.ndarray[ndim=1, dtype=np.int32_t] npoints,
                   np.ndarray[ndim=1, dtype=np.int32_t] istart,
                   np.ndarray[ndim=1, dtype=np.int32_t] xoff,
                   np.ndarray[ndim=1, dtype=np.int32_t] yoff,
                   np.ndarray[ndim=2, dtype=np.float32_t] kern):
        """Load all the arrays computed during the initialization
        for a pzt DM in a sutra_dms object

        :parameters:
            alt: (float) : dm conjugaison altitude

            influ: (np.ndarray[ndim=3,dtype=np.float32_t]) : influence functions

            influpos: (np.ndarray[ndim=1,dtype=np.int32_t]) : positions of the IF

            npoints: (np.ndarray[ndim=1,dtype=np.int32_t]) : for each pixel on the DM screen,
                                                            the number of IF which impact on this pixel

            istart: (np.ndarray[ndim=1,dtype=np.int32_t]) :

            xoff: (np.ndarray[ndim=1,dtype=np.int32_t]) : x-offset

            yoff: (np.ndarray[ndim=1,dtype=np.int32_t]) :y-offset

            kern: (np.ndarray[ndim=1,dtype=np.float32_t]) : convoltuon kernel

        """

        cdef np.ndarray[dtype = np.float32_t] influ_F = influ.flatten("F")
        cdef np.ndarray[dtype = np.int32_t] npoints_F = npoints.flatten("F")

        cdef int inddm = self.dms.get_inddm("pzt", alt)
        if(inddm < 0):
            err = "unknown error whith load_pzt\nDM (pzt" + str(
                alt) + ") doesn't exist"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)

        cdef int ntot_influpos = influpos.shape[0]
        cdef int ntot_istart = istart.shape[0]
        cdef int n = self.dms.d_dms[inddm].influsize * \
            self.dms.d_dms[inddm].influsize

        cdef float * influ2
        influ2 = < float * > malloc(ntot_influpos * sizeof(float))
        cdef tuple_t[float] * influ3 = NULL
        cdef int * influpos2
        influpos2 = < int * > malloc(ntot_influpos * sizeof(int))
        cdef int * istart2
        istart2 = < int * > malloc((ntot_istart + 1) * sizeof(int))

        cdef int i
        for i in range(ntot_influpos):
            influ2[i] = influpos[i] / n

        for i in range(ntot_istart):
            istart2[i] = istart[i]
        istart2[ntot_istart] = istart2[
            ntot_istart - 1] + npoints[ntot_istart - 1]

        # TODO preprocessing: COMPN==2 or 3
        """
#if(COMPN == 2)
      influ3 = new struct tuple_t<float>[ntot_influpos];

      for(int  i = 0; i < ntot_influpos; i++)
    	  influ3[i] = {influpos2[i], influ2[i]};

#elif(COMPN == 3)
      influ3 = new struct tuple_t<float>[ntot_istart * MAXSPOT];

      //For each pixel of screen
      int count = 0;
      for(int  i = 0; i < ntot_istart; i++)
      {
    	  int j = 0;
    	  //For each influence function, cpy the value of postition and influ
    	  for(; j < npoints[i]; j++){
    		  influ3[i * MAXSPOT + j] = {influpos2[count], influ2[count]};
    		  count++;
    	  }
    	  //Fill the last element with 0
    	  for(; j < MAXSPOT; j++)
    		  influ3[i * MAXSPOT + j] = {0, 0.0};
      }
#endif
        """

        self.dms.d_dms[inddm].pzt_loadarrays(< float *> influ_F.data,
                                              influ2,
                                              influ3,
                                              < int * > influpos.data,
                                              influpos2,
                                              < int * > npoints_F.data,
                                              istart2,
                                              < int * > xoff.data,
                                              < int * > yoff.data,
                                              < float * > kern.data)
        free(influ2)
        free(influpos2)
        free(istart2)

    cpdef load_kl(self, float alt,
                  np.ndarray[ndim=1, dtype=np.float32_t] rabas,
                  np.ndarray[ndim=1, dtype=np.float32_t] azbas,
                  np.ndarray[ndim=1, dtype=np.int32_t] ord,
                  np.ndarray[ndim=1, dtype=np.float32_t] cr,
                  np.ndarray[ndim=1, dtype=np.float32_t] cp):
        """Load all the arrays computed during the initialization
        for a kl DM in a sutra_dms object

        :parameters:
            alt: (float) : dm conjugaison altitude

            rabas: (np.ndarray[ndim=1,dtype=np.float32_t]) : TODO

            azbas: (np.ndarray[ndim=1,dtype=np.float32_t]) :

            ord: (np.ndarray[ndim=1,dtype=np.int32_t]) :

            cr: (np.ndarray[ndim=1,dtype=np.float32_t]) :

            cp: (np.ndarray[ndim=1,dtype=np.float32_t]) :
        """

        
        cdef int inddm = self.dms.get_inddm("kl", alt)
        if(inddm < 0):
            err = "unknown error whith load_kl\nDM (kl" + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)


        self.dms.d_dms[inddm].kl_loadarrays(< float *> rabas.data,
                                             < float * > azbas.data,
                                             < int * > ord.data,
                                             < float * > cr.data,
                                             < float * > cp.data)

    cpdef load_tt(self, float alt, np.ndarray[ndim=3, dtype=np.float32_t] influ):
        """Load all the arrays computed during the initialization
        for a tt DM in a sutra_dms object

        :parameters:
            alt: (float) : dm conjugaison altitude

            influ: (np.ndarray[ndim=3,dtype=np.float32_t]) : influence functions
        """
        cdef int inddm = self.dms.get_inddm("tt", alt)
        if(inddm < 0):
            err = "unknown error whith load_tt\nDM (tt" + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        cdef np.ndarray[dtype = np.float32_t] influ_F = influ.flatten("F")
        self.dms.d_dms[inddm].d_influ.host2device(< float *> influ_F.data)


    cpdef set_full_comm(self, np.ndarray[ndim=1, dtype=np.float32_t] comm,
                        bool shape_dm=True):
        """Set the voltage command

            comm: (np.ndarray[ndim=1,dtype=np.float32_t]) : voltage vector

            shape_dm: (bool) : perform the dm_shape after the load (default=True)
        """
        nact_total = self.dms.nact_total()
        if nact_total != comm.size:
            raise ValueError("Incorrect size of voltage vector")

        cdef int comm_index = 0
        cdef float * comm_data = < float * > comm.data
        for inddm in range(self.dms.ndm()):
            self.dms.d_dms[inddm].d_comm.host2device(& comm_data[comm_index])
            comm_index += self.dms.d_dms[inddm].nact()
            if shape_dm:
                self.dms.d_dms[inddm].comp_shape()


    cpdef set_comm(self,bytes type_dm,float alt,
                    np.ndarray[ndim=1,dtype=np.float32_t] comm,
                    bool shape_dm=False):
        """Set the voltage command on a sutra_dm

            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude

            comm: (np.ndarray[ndim=1,dtype=np.float32_t]) : voltage vector

            shape_dm: (bool) : perform the dm_shape after the load (default=False)
        """

        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            err = "unknown error whith set_comm\nDM (" + type_dm + ", " + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        self.dms.d_dms[inddm].d_comm.host2device(< float *> comm.data)
        if shape_dm:
            self.dms.d_dms[inddm].comp_shape()
        if shape_dm:
            self.dms.d_dms[inddm].comp_shape()

    cpdef shape_dm(self, bytes type_dm, float alt):
        """Compute the shape of the DM in a sutra_dm object

            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude
        """
        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            err = "unknown error whith set_comm\nDM (" + type_dm + ", " + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)

        self.dms.d_dms[inddm].comp_shape()

    cpdef computeKLbasis(self, bytes type_dm, float alt,
                         np.ndarray[ndim=1, dtype=np.float32_t] xpos, np.ndarray[ndim=1, dtype=np.float32_t] ypos,
                         np.ndarray[ndim=1, dtype=np.int32_t] indx_pup, long dim, float norm, float ampli):
        """Compute a Karhunen-Loeve basis for the dm:
            - compute the phase covariance matrix on the actuators using Kolmogorov
            - compute the geometric covariance matrix
            - double diagonalisation to obtain KL basis

        :parameters:
            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude

            xpos: (np.ndarray[ndim=1,dtype=np.float32_t]) : x-position of actuators

            ypos: (np.ndarray[ndim=1,dtype=np.float32_t]) : y-position of actuators

            indx_pup: (np.ndarray[ndim=1,dtype=np.int32_t]) : indices of where(pup)

            dim: (long) : number of where(pup)

            norm: (float) : normalization factor

            ampli: (float) : amplitude
        """

        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            err = "unknown error with computeKLbasis function\nDM(" + type_dm + "," + \
                str(alt) + ") doesn't exist"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)

        self.dms.d_dms[inddm].compute_KLbasis(< float *> xpos.data, < float *> ypos.data,
                                              < int * > indx_pup.data, dim, norm, ampli)

    cpdef get_KLbasis(self, bytes type_dm, float alt):
        """Return the klbasis computed by computeKLbasis

        :parameters:
            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude
        :return:
            KLbasis : (np.ndarray(dims=2,dtype=np.float32)) : the KL basis
        """

        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            err = "unknown error with get_KLbasis function DM(" + \
                type_dm + "," + alt + ") doesnt exists"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        cdef const long * dims = self.dms.d_dms[inddm].d_KLbasis.getDims()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data = np.zeros((dims[1], dims[2]), dtype=np.float32)

        self.dms.d_dms[inddm].d_KLbasis.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def get_dm(self, bytes type_dm, float alt):
        """Return the shape of the dm

        :parameters:
            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude
        :return:
            data : (np.ndarray(dims=2,dtype=np.float32)) : DM shape
        """

        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            err = "unknown error with get_dm function DM(" + \
                type_dm + "," + alt + ") doesnt exists"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        cdef const long * dims = self.dms.d_dms[inddm].d_shape.d_screen.getDims()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data = np.zeros((dims[1], dims[2]), dtype=np.float32)

        self.dms.d_dms[inddm].d_shape.d_screen.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    cpdef getComm(self, bytes type_dm, float alt):
        """Return the voltage command of the sutra_dm

        :parameters:
            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude
g        :return:
            data : (np.ndarray(dims=1,dtype=np.float32)) : voltage vector
        """

        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            err = "unknown error with get_dm function DM(" + \
                type_dm + "," + alt + ") doesnt exists"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        cdef const long * dims = self.dms.d_dms[inddm].d_comm.getDims()
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data = np.zeros((dims[1]), dtype=np.float32)

        self.dms.d_dms[inddm].d_comm.device2host(< float *> data.data)
        return data

    cpdef getInflu(self, bytes type_dm, float alt):
        """Return the influence functions of the DM

        :parameters:
            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude
        :return:
            data : (np.ndarray(dims=3,dtype=np.float32)) : influence functions
        """

        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            err = "unknown error with get_dm function DM(" + \
                type_dm + "," + alt + ") doesnt exists"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)

        cdef const long * dims = self.dms.d_dms[inddm].d_influ.getDims()
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] data_F = np.zeros((dims[3], dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] data = np.zeros((dims[1], dims[2], dims[3]), dtype=np.float32)

        self.dms.d_dms[inddm].d_influ.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2], dims[3]))
        return data

    cpdef comp_oneactu(self, bytes type_dm, float alt, int nactu, float ampli):
        """Compute the shape of the dm when pushing the nactu actuator

        :parameters:
            type_dm: (str) : dm type

            alt: (float) : dm conjugaison altitude

            nactu: (int) : actuator number pushed

            ampli: (float) : amplitude
        """

        cdef int inddm = self.dms.get_inddm(type_dm, alt)
        if(inddm < 0):
            err = "unknown error with get_dm function DM(" + \
                type_dm + "," + alt + ") doesnt exists"
            raise ValueError(err)

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.dms.d_dms[inddm].device, 1)
        self.dms.d_dms[inddm].comp_oneactu(nactu, ampli)

    def __str__(self):
        info = "DMs object:\n"
        info += "Contains " + str(self.dms.d_dms.size()) + " DMs:\n"
        info += "DM # | Type  |   Alt   | Nact | Dim\n"
        cdef vector[sutra_dm * ].iterator it_dms = self.dms.d_dms.begin()
        cdef vector[type_screen].iterator it_type = self.dms.d_type.begin()
        cdef sutra_dm * dm
        cdef type_screen ts
        cdef int i = 1
        while(it_dms != self.dms.d_dms.end()):
            dm = deref(it_dms)
            ts = deref(it_type)
            info += "%4d" % i + " | " + "%5s" % < bytes > ts.first + " | " + "%7d" % ts.second + " | " + "%4d" % dm.ninflu + " | " + "%4d" % dm.dim + "\n"
            i = i + 1
            inc(it_dms)
            inc(it_type)
        info += "--------------------------------------------------------"
        return info


cpdef compute_klbasis(Dms g_dm, Param_dm p_dm, Param_geom p_geom, Param_atmos p_atmos, Param_tel p_tel):
    """Compute a Karhunen-Loeve basis for the dm:
            - compute the phase covariance matrix on the actuators using Kolmogorov
            - compute the geometric covariance matrix
            - double diagonalisation to obtain KL basis

    :parameters:
        g_dm: (Dms) : Dms object

        p_dm: (Param_dm) : dm settings

        p_geom: (Param_geom) : geom settings

        p_atmos: (Param_atmos) : atmos settings

        p_tel: (Param_tel) : telescope settings
    """

    cdef int tmp
    if(p_dm.type_dm == "pzt"):
        tmp = (p_geom._ipupil.shape[0] - (p_dm._n2 - p_dm._n1 + 1)) / 2
        tmp_e0 = p_geom._ipupil.shape[0] - tmp
        tmp_e1 = p_geom._ipupil.shape[1] - tmp
        pup = p_geom._ipupil[tmp:tmp_e0, tmp:tmp_e1]
        indx_valid = np.where(pup.flatten("F") > 0)[0].astype(np.int32)
        interactp = p_dm._xpos[1] - p_dm._xpos[0]
        interactm = p_tel.diam / (p_dm.nact - 1)
        p2m = interactm / interactp
        norm = -(p2m * p_tel.diam / (2 * p_atmos.r0)) ** (5. / 3)

        g_dm.computeKLbasis(< bytes > "pzt", p_dm.alt, p_dm._xpos, p_dm._ypos, indx_valid, indx_valid.size, norm, 1.0)
        KLbasis = np.fliplr(g_dm.get_KLbasis(< bytes > "pzt", p_dm.alt))
    else:
        raise TypeError("DM must be pzt type")

    return KLbasis

cpdef computeDMbasis(Dms g_dm, Param_dm p_dm, Param_geom p_geom):
    """Compute a the DM basis as a sparse matrix :
            - push on each actuator
            - get the corresponding dm shape
            - apply pupil mask and store in a column

    :parameters:
        g_dm: (Dms) : Dms object

        p_dm: (Param_dm) : dm settings

        p_geom: (Param_geom) : geom settings
    :return:
        IFbasis = (csr_matrix) : DM IF basis
    """
    cdef int tmp
    tmp = (p_geom._ipupil.shape[0] - (p_dm._n2 - p_dm._n1 + 1)) / 2
    tmp_e0 = p_geom._ipupil.shape[0] - tmp
    tmp_e1 = p_geom._ipupil.shape[1] - tmp
    pup = p_geom._ipupil[tmp:tmp_e0, tmp:tmp_e1]
    indx_valid = np.where(pup.flatten("F") > 0)[0].astype(np.int32)

    #IFbasis = np.ndarray((indx_valid.size, p_dm._ntotact), dtype=np.float32)
    for i in range(p_dm._ntotact):
        g_dm.resetdm(p_dm.type_dm, p_dm.alt)
        g_dm.comp_oneactu(p_dm.type_dm, p_dm.alt, i, 1.0)
        shape = g_dm.get_dm(p_dm.type_dm, p_dm.alt)
        IFvec = csr_matrix(shape.flatten("F")[indx_valid])
        if(i==0):
            val = IFvec.data
            col = IFvec.indices
            row = np.append(0,IFvec.getnnz())
        else:
            val = np.append(val,IFvec.data)
            col = np.append(col,IFvec.indices)
            row = np.append(row,row[-1]+IFvec.getnnz())
    g_dm.resetdm(p_dm.type_dm, p_dm.alt)
    IFbasis = csr_matrix((val,col,row))
    return IFbasis

cpdef computeIFsparse(Dms g_dm, list p_dms, Param_geom p_geom):
    """Compute the influence functions of all DMs as a sparse matrix :
            - push on each actuator
            - get the corresponding dm shape
            - apply pupil mask and store in a column

    :parameters:
        g_dm: (Dms) : Dms object

        p_dms: (Param_dms) : dms settings

        p_geom: (Param_geom) : geom settings
    :return:
        IFbasis = (csr_matrix) : DM IF basis
    """
    ndm = len(p_dms)
    for i in range(ndm):
        IFi = computeDMbasis(g_dm,p_dms[i],p_geom)
        if(i==0):
            val = IFi.data
            col = IFi.indices
            row = IFi.indptr
        else:
            val = np.append(val,IFi.data)
            col = np.append(col,IFi.indices)
            row = np.append(row,row[-1]+IFi.indptr[1:])
    IFsparse = csr_matrix((val,col,row))
    return IFsparse
