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

def make_radii(cobs,nr):
    d = (1.-cobs*cobs)/nr
    rad2 = cobs**2 +d/16.+ d * (np.arange(nr,dtype=np.float32))
    radp= np.sqrt(rad2)
    return radp

#__________________________________________________________________________

def make_kernels(cobs,nr,radp,funct,outscl):
    #make kernels
    #DOCUMENT res=make_kernels(cobs,nr,rad,funct=,outscl=)
    #This routine generates the kernel used to find the KL modes.
    #The  kernel constructed here should be simply a discretization
    #of the continuous kernel. It needs rescaling before it is treated
    #as a matrix for finding  the eigen-values. The outer scale
    #should be in units of the diameter of the telescope.
    nth = 5*nr
    kers  = np.zeros((nth,nr, nr),dtype = np.float32)
    cth = np.cos((np.arange(nth,dtype = np.float32))*(2.*np.pi/nth))
    dth = 2.*np.pi/nth
    fnorm = -1./(2*np.pi*(1.-cobs**2))*0.5
    #the 0.5 is to give  the r**2 kernel, not the r kernel
    for i in range(nr):
        for j in range(i+1):
            te = 0.5*np.sqrt(radp[i]**2+radp[j]**2-(2*radp[i]*radp[j])*cth)
            #te in units of the diameter, not the radius
            if (funct=="kolmo"):
                #DOCUMENT var=kolstf(dvec)
                #This routine returns the kolmogorov phase variance at spatial
                #dimension (inverse of the spatial frequency) dvec
                te = 6.88*te**(5./3.)
                
            elif (funct=="karman"):
                #DOCUMENT var=kolstf(dvec)
                #This routine returns the Von Karman phase variance at spatial
                #dimension (inverse of the spatial frequency) dvec. Same as kolstf
                #but with a correcting factor to account for the outter scale.
                #The latter should be in units of telescope diameter
                te = 6.88*te**(5./3.)*(1-1.485*(te/outscl)**(1./3.)+5.383*(te/outscl)**(2)-6.281*(te/outscl)**(7./3.))
                
            else:
                
                raise TypeError("kl funct error")
                
            f = np.fft.fft(te,axis = -1)
            kelt =  fnorm * dth * np.float32(f.real)
            kers [:,i,j] = kers [:,j,i] = kelt
    return kers

#__________________________________________________________________________ 
#__________________________________________________________________________

def piston_orth(nr):
    s = np.zeros((nr,nr),dtype=np.float32)
    for j in range(nr-1):
        rnm = 1./np.sqrt(np.float32((j+1)*(j+2)))
        s[0:j+1,j] = rnm
        s[j+1,j]= -1*(j+1)*rnm
    
    rnm = 1./np.sqrt(nr)
    s[:,nr-1] = rnm
    return s
    
#__________________________________________________________________________
#__________________________________________________________________________
    
def make_azimuth(nord,npp):
    #DOCUMENT piston_orth(nr)

    azbas = np.zeros((np.int32(1+nord), npp), dtype = np.float32)
    th = np.arange(npp,dtype = np.float32)*(2.*np.pi/ npp)

    azbas [0,:] = 1.0
    for i in np.arange(1,nord,2):
        azbas [np.int32(i),:] = np.cos (((i)/2+1) * th)
    for i in np.arange(2,nord,2):
        azbas [np.int32(i),:] = np.sin (((i)/2) * th)
        
    return azbas
    
#__________________________________________________________________________
#__________________________________________________________________________
    
def radii(nr,npp,cobs):   
    #r = radii(nr,npp,cobs)
    #DOCUMENT res=radii(NumberOfR,NumberOfPhi,Dim)
    #This routine generates an nr x npp array with npp copies of the
    #radial coordinate array. Radial coordinate span the range from
    #r=cobs to r=1 with successive annuli having equal areas (ie, the
    #area between cobs and 1 is divided into nr equal rings, and the
    #points are positioned at the half-area mark on each ring). There
    #are no points on the border.
    
    r2 = cobs**2 +(np.arange(nr,dtype = np.float)+0.)/nr*(1.0-cobs**2)
    rs = np.sqrt(r2)
    r= np.transpose(np.tile(rs,(npp,1)))  
    
    return r
    
#__________________________________________________________________________
#__________________________________________________________________________
    
def polang(r):    
    #p = polang(r)
    #DOCUMENT res=polang(RadialCoordArray)
    #This routine generates an array with the same dimensions as r,
    #but containing the azimuthal values for a polar coordinate system.     
    
    s =  r.shape
    nr = s[0]
    np1 = s[1]
    phi1 = np.arange(np1,dtype = np.float)/float(np1)*2.*np.pi
    p1,p2 = np.meshgrid(np.ones(nr),phi1)
    p = np.transpose(p2)
        
    return p

#__________________________________________________________________________
#__________________________________________________________________________

def setpincs(ax,ay,px,py,cobs):
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
    dcar = (ax[nc-1,0] - ax[0,0]) / (nc-1)
    ofcar = ax[0,0]
    rlx = (px - ofcar)/dcar
    rly = (py - ofcar)/dcar
    lx = np.int32(rlx)
    ly = np.int32(rly)
    shx = rlx - lx
    shy = rly - ly
    
    pincx = np.zeros((4,lx.shape[0],lx.shape[1]))
    pincx[[1,2],:,:] = lx+1
    pincx[[0,3],:,:] = lx
    
    pincy = np.zeros((4,ly.shape[0],ly.shape[1]))
    pincy[[0,1],:,:] = ly
    pincy[[2,3],:,:] = ly+1
    
    pincw = np.zeros((4,shx.shape[0],shx.shape[1]))
    pincw[0,:,:] = (1-shx)*(1-shy)
    pincw[1,:,:] = shx*(1-shy)
    pincw[2,:,:] = shx*shy
    pincw[3,:,:] = (1-shx)*shy
    
      
    axy = ax**2 + ay**2
    axyinap = np.clip(axy,cobs**2.+1.e-3,0.999)
    #sizeaxyinap=axyinap.shape[1]# not used
    
    #pincw = pincw*axyinap[pincx+(pincy-1)*sizeaxyinap] --->
    
    for z in range(pincw.shape[0]):
        for i in range(pincw.shape[1]):
            for j in range(pincw.shape[2]):
                pincw[z,i,j] = pincw[z,i,j]*axyinap[np.int32(pincx[z,i,j]),np.int32(pincy[z,i,j])]
    
    pincw = pincw*np.tile(1.0/np.sum(pincw,axis=0),(4,1,1))   

    return pincx,pincy,pincw
    
#__________________________________________________________________________
#__________________________________________________________________________

def pcgeom(nr,npp,cobs,ncp,ncmar):
    #pcgeom,bas,ncp,ncmar;
    #DOCUMENT pcgeom,&geom,ncp,ncmar
    #This routine builds a geom_struct. px and py are the x and y
    #coordinates of points in the polar arrays.  cr and cp are the
    #r and phi coordinates of points in the cartesian grids. ncmar
    #allows the possibility that there is a margin of ncmar points
    #in the cartesian arrays outside the region of interest
    nused = ncp - 2*ncmar
    ff = 0.5*nused
    hw = np.float(ncp-1)/2.
    
    r = radii(nr,npp,cobs)
    p = polang(r)
    
    
    px0 = r * np.cos(p)
    py0 = r * np.sin(p)
    px = ff * px0 + hw
    py = ff * py0 + hw
    ax = np.reshape(np.arange(long(ncp)**2,dtype=np.float)+1,(long(ncp),long(ncp)),order='F')
    ax = np.float32(ax-1) % ncp - 0.5 * (ncp-1)
    ax = ax / (0.5 * nused)
    ay = np.transpose(ax)

    
    pincx,pincy,pincw = setpincs(ax,ay,px0,py0,cobs)
    
    dpi = 2 * np.pi
    cr2 = (ax**2 + ay**2)
    ap = np.clip(cr2,cobs**2 + 1.e-3,0.999)
    #cr = (cr2 - cobs**2) / (1 - cobs**2) * nr - 0.5; 
    cr = (cr2 - cobs**2) / (1 - cobs**2) * nr
    cp = (np.arctan2(ay,ax)+dpi)%dpi
    cp = (npp / dpi) * cp
    
    cr = np.clip(cr,1.e-3,nr-1.001)
    #fudge -----, but one of the less bad ones
    cp = np.clip(cp,1.e-3,npp -1.001)
    #fudge -----  this is the line which
    #gives that step in the cartesian grid
    #at phi = 0.
    return ncp,ncmar,px,py,cr, cp,pincx,pincy,pincw,ap
    
#__________________________________________________________________________
#__________________________________________________________________________
    
def set_pctr(dim,nr,npp,nkl,cobs,nord, ncmar=None, ncp=None):
    #set_pctr,klbasis, ncp= dim
    #DOCUMENT geom=set_pctr(bas, ncp =, ncmar=)
    #This routine calls pcgeom to build a geom_struct with the
    #right initializations. bas is a gkl_basis_struct built with
    #the gkl_bas routine.
    ncp = dim
    if(ncmar==None):ncmar = 2
    if(ncp==None):ncp = 128
    ncp,ncmar,px,py,cr, cp,pincx,pincy,pincw,ap = pcgeom(nr,npp,cobs,ncp,ncmar)
    return ncp,ncmar,px,py,cr, cp,pincx,pincy,pincw,ap
    
#__________________________________________________________________________

#__________________________________________________________________________
    
def gkl_fcom(kers,cobs,nf):
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
    fktom =  (1.-cobs**2)/nr
    #fevtos = np.sqrt(2*nr) #not used
       
    evs = np.zeros((nr, nt), dtype = np.float32)
    #ff isnt used - the normalisation for
    #the eigenvectors is straightforward:
    #integral of surface**2 divided by area = 1,
    #and the cos**2 term gives a factor
    #half, so multiply zero order by
    #sqrt(n) and the rest by sqrt (2n)
    
    #zero order is a special case...
    #need to deflate to eliminate infinite eigenvalue - actually want
    #evals/evecs of zom - b where b is big and negative    
    zom = kers[0,:,:]
    s = piston_orth(nr)
    
    ts = np.transpose(s)
    #b1 = ((ts(,+)*zom(+,))(,+)*s(+,))(1:nr-1, 1:nr-1)
    btemp = (ts.dot(zom).dot(s))[0:nr-1,0:nr-1]
    
    #newev = SVdec(fktom*b1,v0,vt)
    v0, newev, vt = np.linalg.svd(fktom*btemp, full_matrices=True)
    
    v1 = np.zeros((nr,nr),dtype=np.float32)
    v1[0:nr-1,0:nr-1] = v0
    v1[nr-1,nr-1] = 1
    
    vs = s.dot(v1)
    newev = np.concatenate((newev,[0]))
    print np.size(newev)
    evs[:,nxt] = np.float32(newev)
    kers [nxt,:,:] = np.sqrt(nr)*vs
    
    nxt = 1
    while True:
        vs, newev, vt = np.linalg.svd(fktom*kers[nxt,:,:],full_matrices=True) 
        #newev = SVdec(fktom*kers(,,nxt),vs,vt)
        evs[:,nxt] = np.float32(newev)
        kers[nxt,:,:] = np.sqrt(2.*nr)*vs
        mxn = max(np.float32(newev))
        egtmxn = np.floor(evs[:, 0:nxt+1]>mxn)
        nxt = nxt + 1
        if ((2*np.sum(egtmxn)-np.sum(egtmxn[:,0])) >= nkl):
            break
    
    
    nus = nxt - 1
    kers = kers[0:nus+1,:,:]
    
    #evs = reform (evs [:, 1:nus], nr*(nus))
    
    evs = np.reshape(evs[:,0:nus+1],nr*(nus+1),order='F')
    a = np.argsort(-1.*evs)[0:nkl]
    
    #every eigenvalue occurs twice except
    #those for the zeroth order mode. This
    #could be done without the loops, but
    #it isn't the stricking point anyway...
    
    no = 0
    ni = 0
    #oind = array(long,nf+1)
    oind = np.zeros(nkl+1,dtype=np.int32)
    
    
    while True:
        if (a[ni] < nr):
            oind[no] = a[ni]
            no = no + 1
        else:
            oind[no] = a[ni]
            oind[no+1] = a[ni]
            no = no + 2
        
        ni = ni + 1
        if (no >= (nkl)):
            break
    
    
    
    oind = oind[0:nkl]
    tord = (oind)/nr+1
    #odd = ((long(indgen(nf)-1) % 2) == 1)
    odd = np.arange(nkl,dtype=np.int32)%2
    pio = (oind) % nr +1
    
    evals = evs[oind]
    ordd = 2 *(tord-1) - np.floor((tord>1) & (odd))+1
    
    nord = max(ordd)
    
    
    rabas = np.zeros((nr,nkl),dtype = np.float32)
    sizenpo=np.int32(nord)
    npo = np.zeros(sizenpo,dtype=np.int32)
    
      
    for i in range(nkl):
        npo[np.int32(ordd[i])-1] = npo[np.int32(ordd[i])-1] + 1
        rabas[:,i] = kers [tord[i]-1,:,pio[i]-1]
    
    return evals,nord,npo,ordd,rabas

#__________________________________________________________________________


def gkl_sfi(_klbas,i):
    #DOCUMENT 
    #This routine returns the i'th function from the generalised KL
    #basis bas. bas must be generated first with gkl_bas.
    nr = _klbas.nr
    npp = _klbas.npp
    ordd = _klbas.ordd
    rabas = _klbas.rabas
    azbas = _klbas.azbas
    nfunc = _klbas.nfunc

    if (i>nfunc-1):
        raise TypeError("kl funct order it's so big")
        
    else:    
    
        ordi = np.int32(ordd[i])
        rabasi=rabas[:,i]
        azbasi=np.transpose(azbas)
        azbasi = azbasi[ordi,:]
    
        sf1 = np.zeros((nr,npp),dtype=np.float64)
        for j in range(npp):
            sf1[:,j]=rabasi
    
        sf2 = np.zeros((npp,nr),dtype=np.float64)
        for j in range(nr):
            sf2[:,j]=azbasi
        
        sf = sf1*np.transpose(sf2)
        
        return sf
        
def pol2car(pol,_klbas,mask=0):
    # DOCUMENT cart=pol2car(cpgeom, pol, mask=)
    # This routine is used for polar to cartesian conversion.
    # pol is built with gkl_bas and cpgeom with pcgeom.
    # However, points not in the aperture are actually treated
    # as though they were at the first or last radial polar value
    # -- a small fudge, but not serious  ?*******
    #cd = interpolate.interp2d(cr, cp,pol)
    ncp = _klbas.ncp
    cr = _klbas.cr
    cp = _klbas.cp    
    nr = _klbas.nr
    npp = _klbas.npp

    r = np.arange(nr,dtype=np.float64)
    phi = np.arange(npp,dtype=np.float64)
    tab_phi,tab_r = np.meshgrid(phi,r)
    tab_x = (tab_r/(nr))*np.cos((tab_phi/(npp))*2*np.pi)
    tab_y = (tab_r/(nr))*np.sin((tab_phi/(npp))*2*np.pi)
    
    newx = np.linspace(-1,1,ncp)
    newy = np.linspace(-1,1,ncp)
    tx, ty = np.meshgrid(newx,newy)
    
    cd = interpolate.griddata((tab_r.flatten(),tab_phi.flatten()),pol.flatten(),(cr,cp),method='cubic')
    cdf = interpolate.griddata((tab_r.flatten("F"),tab_phi.flatten("F")),pol.flatten("F"),(cr,cp),method='cubic')
    cdxy = interpolate.griddata((tab_y.flatten(),tab_x.flatten()),pol.flatten(),(tx,ty),method='cubic')
    
    if(mask == 1):
        ap = _klbas.ap       
        cd = cd*(ap)
        cdf = cdf*(ap)
        cdxy = cdxy*(ap)
        
    return cd,cdf,cdxy
    
def kl_view(_klbas,mask=1):
    
    nfunc = _klbas.nfunc
    ncp = _klbas.ncp
    
    tab_kl = np.zeros((nfunc,ncp,ncp),dtype=np.float64)
    tab_klf = np.zeros((nfunc,ncp,ncp),dtype=np.float64)
    tab_klxy = np.zeros((nfunc,ncp,ncp),dtype=np.float64)
    
    for i in range(nfunc):
        
        tab_kl[i,:,:],tab_klf[i,:,:],tab_klxy[i,:,:] = pol2car(gkl_sfi(_klbas,i),_klbas,mask)

    
    return tab_kl,tab_klf,tab_klxy
    
    

# --------------------------------------- 
"""
python kl
@author: unknown
Broken translation start
"""


#import numpy as np
#cimport numpy as np
#np.import_array()
#
#cdef class Param_kl_basis:
#    def __cinit__(int dim, long nfunc=500, float cobs=0, long Nr=-1, long Np=-1,
#                bytes funct= < bytes > "kolmo", outscl=None):
#        # TODO check funct=bytes
#        if(Nr == -1):
#            Nr = long(5.0 * np.sqrt(nfunc))
#        if(Np == -1):
#            Np = long(5 * Nr)
#
#        cdef np.ndarray[ndim = 1, dtype = np.float32_t] radp = make_radii(Nr, cobs)
#        cdef np.ndarray[ndim = 3, dtype = np.float32_t] kers = \
#                                    make_kernels(cobs, Nr, radp, funct, outscl)
#
#
#"""
#func make_klbas(nfunc,cobs,dim,nr=,np=,funct=,outscl=)
#  /*DOCUMENT make_klbas(nfunc,cobs,nr=,np=,funct=,outscl=)
#  SEE ALSO : 
#   */
#{
#  if (cobs == []) cobs = 0;
#  if (nfunc == []) nfunc = 500L;
#  
#  if (!is_set(nr)) nr = long(5.0f*sqrt(nfunc));
#  if (!is_set(np)) np = long(5*nr);
#  if (dimsof(funct)==[]) funct="kolmo";
#  
#  radp = make_radii(nr, cobs);
#
#  kers = make_kernels(cobs, nr, radp,funct=funct,outscl=outscl);
# --------------------------------------- 
#  gkl_fcom,kers,cobs,nfunc,evals,nord,npo,ord,rabas;
#
#  azbas = make_azimuth(nord, np);
# 
#  klbasis = kl_basis_struct();
#  klbasis.nr=nr;
#  klbasis.np=np;
#  klbasis.nfunc=nfunc; 
#  klbasis.cobs=cobs;
#  klbasis.radp=&radp;
#  klbasis.evals=&evals;
#  klbasis.nord=nord;
#  klbasis.npo=&npo;
#  klbasis.ord=&ord;
#  klbasis.rabas=&rabas;
#  klbasis.azbas=&transpose(azbas);
#  
#  set_pctr,klbasis, ncp= dim;
#  
#  return klbasis;
#}
#"""
#
#
#
#cdef make_radii(int nr, float cobs):
#    """  This routine generates an nr x np array with np copies of the
#  radial coordinate array. Radial coordinate span the range from
#  r=cobs to r=1 with successive annuli having equal areas (ie, the
#  area between cobs and 1 is divided into nr equal rings, and the
#  points are positioned at the half-area mark on each ring). There
#  are no points on the border.
#    """
#    # TODO check type /nr
#    cdef float d = (1. - cobs * cobs) / nr
#
#    return np.sqrt(cobs ** 2 + d / 16. + d * np.arange(nr, dtype=np.float32))
#
#
#
#cdef make_kernels(float cobs, int nr, np.ndarray[ndim=1, dtype=np.float32_t] rad,
#                    bytes funct , outscl=None):
#    """This routine generates the kernel used to find the KL modes.
#  The  kernel constructed here should be simply a discretization
#  of the continuous kernel. It needs rescaling before it is treated
#  as a matrix for finding  the eigen-values. The outter scale
#  should be in units of the diameter of the telescope.
#    """
#
#    cdef int nth = 5 * nr
#    cdef np.ndarray[ndim = 3, dtype = np.float] kers = np.zeros((nr, nr, nth), dtype=np.float32)
#    cdef np.ndarray[ndim = 1, dtype = np.float] cth = np.cos(np.arange(nth) * 2.*np.pi / nth)
#    cdef float dth = 2.*np.pi / nth
#    cdef float fnorm = -1. / (2 * np.pi * (1. - cobs * cobs)) * 0.5
#    # the 0.5 is to give  the r^2 kernel, not
#    # the r kernel
#    cdef int i, j
#    cdef float te
#
#    for i in range(nr):
#        for j in range(i):
#            te = 0.5 * np.sqrt(rad[i] ** 2 + rad[j] ** 2 - (2 * rad[i] * rad[j]) * cth)
#            # te in units of the diameter, not the radius
#            if (funct == "kolmo"):  te = kolstf(te)
#            if (funct == "karman"): te = karmanstf(te, outscl=outscl)
#            if ((funct != "kolmo") and (funct != "karman")):
#                raise ValueError("The statistics is not known !")
#
#            kelt = fnorm * dth * (np.fft.ifft(te, -1))
#            kers [i, j, :] = kelt
#            kers [j, i, :] = kelt
#
#    return kers
#
#
#
#cdef kolstf(np.ndarray[ndim=1, dtype=np.float32_t] dvec):
#    """This routine returns the kolmogorov phase variance at spatial
#  dimension (inverse of the spatial frequency) dvec"""
#    return 6.88 * dvec ** (5. / 3.)
#
#
#cdef karmanstf(float dvec, int outscl=3):
#    """This routine returns the Von Karman phase variance at spatial
#  dimension (inverse of the spatial frequency) dvec. Same as kolstf
#  but with a correcting factor to account for the outter scale.
#  The latter should be in units of telescope diameter
#  """
#    return  6.88 * dvec ** (5. / 3.) * (1 - 1.485 * (dvec / outscl) ** (1. / 3.) + \
#                              5.383 * (dvec / outscl) ** (2) - 6.281 * \
#                              (dvec / outscl) ** (7. / 3.))
#
#
#
#
#
#cdef gkl_fcom(np.ndarray[ndim=3, dtype=np.float32_t] kers, float cobs, float nf):
#    """This routine does the work : finding the eigenvalues and
#  corresponding eigenvectors. Sort them and select the right
#  one. It returns the KL modes : in polar coordinates : rabas
#  as well as the associated variance : evals. It also returns
#  a bunch of indices used to recover the modes in cartesian
#  coordinates (nord, npo and ord).
#  """
#
#    cdef int nr = kers.shape[0]
#    cdef int nt = kers.shape[3]
#    cdef int nxt = 0
#    cdef float fktom = (1. - cobs ** 2) / nr
#    cdef float fvtos = np.sqrt(2 * nr)
#    cdef np.ndarray[ndim = 2, dtype = np.float32_t] evs = np.zeros((nr, nt), dtype=np.float32)
#    cdef np.ndarray[ndim = 2, dtype = np.float32_t] s = piston_orth(nr)
#
#    cdef np.ndarray[ndim = 2, dtype = np.float32_t] b1, v0, vt, v1, vs
#    cdef np.ndarray[ndim = 2, dtype = np.int32_t] egtmxn
#    cdef np.ndarray[ndim = 1, dtype = np.float32_t] newev
#
#    cdef float do
#
#    # ff isnt used - the normalisation for
#    # the eigenvectors is straightforward:
#    # integral of surface^2 divided by area = 1,
#    # and the cos^2 term gives a factor
#    # half, so multiply zero order by
#    # sqrt(n) and the rest by sqrt (2n)
#
#    # zero order is a special case...
#    # need to deflate to eliminate infinite eigenvalue - actually want
#    # evals/evecs of zom - b where b is big and negative
#
#    b1 = np.dot(np.dot(s.T, kers[:, :, 0]), s)[0:nr - 1, 0:nr - 1]
#
#    v0, newev, vt = np.linalg.svd(fktom * b1)
#
#    v1 = np.zeros((nr, nr), dtype=np.float32)
#    v1[0:nr - 1, 0:nr - 1] = v0
#    v1[nr - 1, nr - 1] = 1
#
#    vs = np.dot(s, v1)
#    evs[:newev.shape[0], nxt] = newev
#    kers[:, :, nxt] = np.sqrt(nr) * vs
#
#    nxt = 1
#    do = nf - 1
#    while(do < nf):
#        vs, newev, vt = np.linalg.svd(fktom * kers[:, :, nxt])
#        evs[:newev.shape[0], nxt] = newev
#        kers[:, :, nxt] = np.sqrt(2 * nr) * vs
#        egtmxn = (evs[:, :nxt] > np.max(newev)).astype(np.int32)
#        do = 2 * np.sum(egtmxn) - sum(egtmxn[:, 0])
#        nxt += 1
#
#    kers = kers[:, :, :nxt]
#    evs = evs[:, :nxt].flatten()
#
#    cdef np.ndarray[ndim = 1, dtype = np.int32_t] a = np.argsort(-evs)[:nf]
#    # every eigenvalue occurs twice except
#    # those for the zeroth order mode. This
#    # could be done without the loops, but
#    # it isn't the stricking point anyway...
#    no = 1
#    ni = 1
#    cdef np.ndarray[ndim = 1, dtype = np.int32_t]oind = np.zeros(nf + 1, dtype=np.int32)
#    while(no < nf + 1):
#        if(a[ni] < nr + 1):
#            oind[no] = a[ni]
#            no += 1
#        else:
#            oind[no] = a[ni]
#            oind[no + 1] = a[ni]
#            no += 2
#        ni += 1
#
#    oind = oind[:nf]
#    cdef np.ndarray[ndim = 1, dtype = np.float32_t] tord = (oind - 1) / nr + 1
#    cdef np.ndarray[ndim = 1, dtype = np.int32_t] odd = np.arange(nf, dtype=np.int32) % 2
#    cdef np.ndarray[ndim = 1, dtype = np.int32_t] pio = (oind - 1) % nr + 1
#
#
#    cdef np.ndarray[ndim = 1, dtype = np.float32_t ] evals = evs[oind]
#    cdef np.ndarray[ndim = 1, dtype = np.int32_t] ord = 2 * (tord - 1) - np.floor((tord > 1) * odd) + 1
#    cdef int nord = np.max(ord)
#    cdef np.ndarray[ndim = 2, dtype = np.float32_t] rabas = np.zeros((nr, nf), dtype=np.float32)
#    cdef np.ndarray[ndim = 1, dtype = np.int32_t] npo = np.zeros(nord, dtype=np.int32)
#
#    cdef int i
#    for i in range(long(nf)):
#        npo[ord[i]] = npo[ord[i]] + 1
#        rabas[:, i] = kers[:, pio[i], tord[i]]
#
#    return evals, nord, npo, ord, rabas
#
#
#cdef make_azimuth(nord, Np):
#    cdef np.ndarray[ndim = 2, dtype = np.float32_t] azi = np.zeros((int(nord + 1), Np),
#                                                    dtype=np.float32)
#    cdef np.ndarray[ndim = 1, dtype = np.float32_t] th = np.arange(Np) * (2.*np.pi / Np)
#
#    cdef int i
#    azi[0, :] = 1.
#
#    for i in range(1, nord + 1, 2):
#        azi[i, :] = np.cos((i / 2 + 1) * th)
#        azi[i + 1, :] = np.sin((i + 1) / 2 * th)
#
#    return azi
#
#
#
#"""
#func make_azimuth(nord, np)
#  /*DOCUMENT piston_orth(nr)
#   */
#{
#  azi = array(float,[2,long(1+nord), np]);
#  th = float(indgen(np)-1)*(2.*pi/ np);
#
#  #TODO check this out
#  azi (1,) = 1.0;
#  for (i = 2; i<=nord;i+=2)  azi (i,) = cos (((i-1)/2+1) * th);
#  for (i = 3; i<=nord;i+=2)  azi (i,) = sin (((i-1)/2) * th);
#  return azi;
#}"""
#
#
#cdef piston_orth(int Nr):
#    cdef np.ndarray[ndim = 2, dtype = np.float32_t] s = np.zeros((), dtype=np.float32)
#    cdef int i
#    cdef float rnm
#
#    for i in range(Nr - 1):
#        nrm = 1. / np.sqrt(i * (i + 1))
#        s[0:i + 1, i] = rnm
#        s[i + 1, i] = -i * rnm
#    rnm = 1. / np.sqrt(Nr)
#    s[:Nr - 1] = rnm
#
#    return s
