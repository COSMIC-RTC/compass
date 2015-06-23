#################################################
# P-Class (parametres) Param_dm
#################################################
cdef class Param_dm:

    def set_type(self,bytes t):
        """set attribute type_dm to t
        t -- string: type of dm
        """
        self.type_dm=t

    def set_nact(self,long n):
        """set attribute nact to n
        n -- long: number of actuators in the diameter
        """
        self.nact=n

    def set_alt(self, float a):
        """set attribute alt to
        a -- float :conjugaison altitude (im m)
        """
        self.alt=a

    def set_thresh(self, float t):
        """set attribute thresh to
        t -- float: threshold on response for selection (<1)
        """
        self.thresh=t

    def set_coupling(self, float c):
        """set attribute coupling to
        c -- float: actuators coupling (<0.3)
        """
        self.coupling=c

    def set_unitpervolt(self,float u ):
        """set attribute unitpervolt to
        u -- float: Influence function sensitivity in unit/volt
        """
        self.unitpervolt=u

    def set_push4imat(self, p):
        """set attribute push4imat to
        p -- float : nominal voltage for imat
        """
        self.push4imat=p

    def set_ntotact(self, long n):
        """set attribute ntotact to
        n -- long :
        """
        self._ntotact=n

    def set_xpos(self,np.ndarray[ndim=1,dtype=np.float32_t] xpos):
        self._xpos=xpos

    def set_ypos(self,np.ndarray[ndim=1,dtype=np.float32_t] ypos):
        self._ypos=ypos

    def set_i1(self,np.ndarray[ndim=1,dtype=np.int32_t] i1):
        self._i1=i1

    def set_j1(self,np.ndarray[ndim=1,dtype=np.int32_t] j1):
        self._j1=j1

    def set_influ(self, np.ndarray[ndim=3,dtype=np.float32_t] influ):
        self._influ=influ

#max_extent signature
cdef _dm_init(Dms dms, Param_dm p_dms, Param_wfs p_wfs, Param_geom p_geom, Param_tel p_tel,int *max_extent):
            """ inits a Dms object on the gpu
            p_dms
            p_wfs
            p_geom
            p_tel


            """

            cdef float patchDiam
            cdef int extent#, max_extent
            cdef long dim
            cdef float disp=0
            cdef long ninflu,influsize,ninflupos,n_npts
            cdef long _nr, _np
            
            cdef float tmp

    #if(len(p_dms)!=0):
    #    dms=Dms(len(p_dms))
    #    for i in range(len(p_dms)):
            if(p_dms.pupoffset is not None):
                p_dms.puppixoffset=p_dms.pupoffset/p_tel.diam*p_geom.pupdiam

            if(p_dms.type_dm=="pzt"):
                #find out the support dimension for the given mirror.
                patchDiam = p_geom.pupdiam#+2*np.max(np.abs(p_wfs.xpos,p_wfs.ypos))* \
                         #4.848e-6*np.abs(p_dms[i].alt)/p_tel.diam*p_geom.pupdiam
                #Patchdiam
                p_dms._pitch=long(patchDiam/(p_dms.nact-1))

                extent=p_dms._pitch*(p_dms.nact+3) # + 1.5 pitch each side
                p_dms._n1=np.floor(p_geom.cent-extent/2)
                p_dms._n2=np.ceil(p_geom.cent+extent/2)
                if( p_dms._n1<1): p_dms._n1=1
                if( p_dms._n2>p_geom.ssize):p_dms._n2=p_geom.ssize

            else:
                # we are dealing with a TT, or kl
                extent=p_geom.pupdiam+16
                p_dms._n1=np.floor(p_geom.cent-extent/2)
                p_dms._n2=np.ceil(p_geom.cent+extent/2)
                if( p_dms._n1<1): p_dms._n1=1
                if( p_dms._n2>p_geom.ssize):p_dms._n2=p_geom.ssize

            #max_extent
            max_extent[0]=max(max_extent[0],p_dms._n2-p_dms._n1+1)
            if(p_dms.alt==0):
                if(p_dms.type_dm=="tt"):
                    # TODO check next line (different max used ?)
                    # (max among all the p_dms, previous p_dms._n modified)
                    # adapted: #max_extent
                    #max_extent=max(p_dms._n2-p_dms._n1+1)
                    extent=int(max_extent[0]*1.05)
                    if (extent %2 != 0): extent +=1
                    p_dms._n1=np.floor(p_geom.cent-extent/2)
                    p_dms._n2=np.ceil(p_geom.cent+extent/2)
                    if( p_dms._n1<1): p_dms._n1=1
                    if( p_dms._n2>p_geom.ssize):p_dms._n2=p_geom.ssize

            if( p_dms.type_dm=="pzt"):
                make_pzt_dm(p_dms,p_geom,disp)
                dim = max(p_dms._n2-p_dms._n1+1, p_geom._mpupil.shape[0])
                ninflu=long(p_dms._ntotact)
                influsize=long(p_dms._influsize)
                ninflupos=long(p_dms._influpos.size)
                n_npts = long(p_dms._ninflu.size)

                dms.add_dm(p_dms.type_dm, p_dms.alt, dim, ninflu, influsize,
                            ninflupos, n_npts,p_dms.push4imat)
                dms.load_pzt(p_dms.alt, p_dms._influ, p_dms._influpos.astype(np.int32),
                            p_dms._ninflu, p_dms._influstart, p_dms._i1,p_dms._j1)

            elif(p_dms.type_dm=="tt"):
                dim = long(p_dms._n2-p_dms._n1+1)
                make_tiptilt_dm(p_dms, p_wfs, p_geom, p_tel)
                dms.add_dm(p_dms.type_dm, p_dms.alt, dim, 2, dim,
                            1, 1, p_dms.push4imat)
                dms.load_tt(p_dms.alt,p_dms._influ)
                

            elif(p_dms.type_dm=="kl"):
                dim=long(p_dms._n2-p_dms._n1+1)
                #TODO make_kl_dm, n;
                ninflu=p_dms.nkl
                influsize=long(p_dms._klbas.ncp)
                _nr=long(p_dms._klbas.nr)
                _np=long(p_dms._klbas.np)
                dms.add_dm(p_dms.type_dm, p_dms.alt, dim, ninflu, influsize,
                            _nr, _np,p_dms.push4imat)
                dms.load_kl(p_dms.alt,p_dms._klbas.rabas,p_dms._klbasazbas,
                            p_dms._klbas.ord,p_dms._klbas.cr,p_dms._klbas.cp)

                #Verif
                #res1 = pol2car(*y_dm(n)._klbas,gkl_sfi(*y_dm(n)._klbas, 1));
                #res2 = yoga_getkl(g_dm,0.,1);


"""
func dm_init(void)
/* DOCUMENT dm_init
   dm_init
   inits a yDM object on the gpu
   no input parameters
   requires 3 externals : 
   y_ geom   : a y_struct for the geometry
   y_ tel    : a y_struct for the telescope
   y_dm    : a y_struct for the dm
   creates 1 external :
   g_dm    : a yDMs object on the gpu
   SEE ALSO:
*/
{
  extern y_dm,g_dm;
  if (y_dm != []) {
    g_dm = yoga_dms(numberof(y_dm));
    for (n=1;n<=numberof(y_dm);n++) {
      if (y_dm(n).type == "pzt") {
      if (y_dm(n).type == "tt") {
      if (y_dm(n).type == "kl") {

      // Florian features
      if (y_dm(n).type == "flo_kl") {
        dim       = long(y_dm(n)._n2-y_dm(n)._n1+1);
        make_flo_kl_dm, n;
        //error;
        
        dim       = long(y_dm(n)._n2-y_dm(n)._n1+1);
        ninflu    = long(y_dm(n).nkl);
        influsize = 0;
        nr        = 0;
        np        = 0;

        yoga_addkl,g_dm,float(y_dm(n).alt),dim,ninflu,influsize,nr,np,float(y_dm(n).push4imat);

        yoga_floloadkl,g_dm,*(*y_dm(n)._klbas).covmat,*(*y_dm(n)._klbas).filter,*(*y_dm(n)._klbas).evals,*(*y_dm(n)._klbas).bas,float(y_dm(n).alt) ;
        //error;
        /*
        // verif :
        res1 = pol2car(*y_dm(n)._klbas,gkl_sfi(*y_dm(n)._klbas, 1));
        res2 = yoga_getkl(g_dm,0.,1);
        */
      }
    }
  }
}
"""

def dm_init(p_dms, Param_wfs p_wfs, Param_geom p_geom, Param_tel p_tel):
    #max_extent
    cdef int max_extent=0
    if(len(p_dms)!=0):
        dms=Dms(len(p_dms))
        for i in range(len(p_dms)):
             #max_extent
            _dm_init(dms, p_dms[i], p_wfs, p_geom, p_tel,&max_extent)
    return dms





cdef make_pzt_dm(Param_dm p_dm,Param_geom geom,disp):

    cdef int i
    #best parameters, as determined by a multi-dimensional fit
    #(see coupling3.i)
    cdef float p1,p2,irc, coupling
    coupling=p_dm.coupling

    p1=4.49469+7.25509*coupling+(-32.1948)*coupling**2+17.9493*coupling**3

    p2=2.49456+(-0.65952)*coupling+8.78886*coupling**2+(-6.23701)*coupling**3

    irc=1.16136+2.97422*coupling+(-13.2381)*coupling**2+20.4395*coupling**3

    #cdef long size=geom.ssize
    cdef long nxact=p_dm.nact
    cdef float cent=geom.cent
    cdef long pitch=p_dm._pitch
    cdef float ir=irc*pitch

    cdef float tmp_c=pitch/np.abs(ir)
    cdef float c = (coupling - 1.+ tmp_c**p1)/(np.log(tmp_c)*tmp_c**p2)

    # compute IF on partial (local) support:
    cdef long smallsize=np.ceil(2*ir+10)
    if(smallsize%2!=0):smallsize+=1
    p_dm._influsize=smallsize

    cdef np.ndarray[ndim=2,dtype=np.float32_t] tmpx, tmp
    tmpx=np.tile(np.arange(smallsize,dtype=np.float32)-smallsize/2+0.5,(smallsize,1))
    tmpx=np.abs(tmpx/ir)
    #clip
    tmpx[tmpx<1e-8]=1e-8
    tmpx[tmpx>2]=2.
    tmp=(1.-tmpx**p1+c*np.log(tmpx)*tmpx**p2)*(1.-tmpx.T**p1+c*np.log(tmpx.T)*tmpx.T**p2)

    #influ   = tmp*(tmpx <= 1.)*(tmpy <= 1.);
    tmp= tmp*(tmpx <= 1.)*(tmpx.T <= 1.)

    # compute location (x,y and i,j) of each actuator:
    cdef np.ndarray[ndim=3,dtype=np.float32_t] cub = np.zeros((nxact,nxact,2),dtype=np.float32)
    # make X and Y indices array:
    cub[:,:,0]=np.tile(np.arange(nxact),(nxact,1))+1
    cub[:,:,1]=np.tile(np.arange(nxact),(nxact,1)).T+1
    # express "centered" coordinate of actuator in pixels
    cub=(cub-1.-(nxact-1.)/2.)*pitch
    #the following determine if an actuator is to be considered or not
    # relative to the pitchmargin parameter.
    cdef np.ndarray dis=np.sqrt(cub[:,:,0]**2+cub[:,:,1]**2)

    cdef float pitchMargin=1.44
    if(p_dm.margin!=0):pitchMargin=p_dm.margin

    cdef rad=((nxact-1.)/2.+pitchMargin)*pitch
    # 1 if valid actuator, 0 if not:
    inbigcirc= np.where(dis < rad)

    # converting to array coordinates:
    cub += cent
    # TODO cub      = cub(*,);
    # cub now has two indices: first one is actuator number
    # second one is: 1:Xcoord, 2:Ycoord

    #filtering actuators outside of a disk radius = rad (see above)
    cdef np.ndarray cubval = cub[inbigcirc]
    p_dm._ntotact=cubval.shape[0]
    p_dm._xpos=cubval[:,0]
    p_dm._ypos=cubval[:,1]
    p_dm._i1=(cubval[:,0]-smallsize/2+0.5-p_dm._n1).astype(np.int32)
    p_dm._j1=(cubval[:,1]-smallsize/2+0.5-p_dm._n1).astype(np.int32)

    # influ = influ(,,-)*array(1.f,y_dm(nm)._ntotact)(-,-,);
    # tmp= tmp*(tmpx <= 1.)*(tmpx.T <= 1.)
    # influ = tmp(,,-)*array(1.f,y_dm(nm)._ntotact)(-,-,);
    cdef np.ndarray[ndim=3,dtype=np.float32_t] influ=np.zeros((tmp.shape[0],tmp.shape[1],p_dm._ntotact),dtype=np.float32)
    for i in range(p_dm._ntotact):
        influ[:,:,i]=tmp

    if(p_dm._puppixoffset is not None):
        p_dm._xpos +=p_dm._puppixoffset[0]
        p_dm._ypos +=p_dm._puppixoffset[1]
        #dm._i1
        #dm._j1

        influ=influ*float(p_dm.unitpervolt/np.max(influ))
        p_dm._influ=influ

        comp_dmgeom(p_dm,geom)

    cdef long dim=max(geom._mpupil.shape[0],p_dm._n2-p_dm._n1+1)
    cdef float off=(dim-p_dm._influsize)/2

    cdef np.ndarray[ndim=2, dtype=np.float32_t] kernconv=np.zeros((dim,dim),dtype=np.float32)
    kernconv [off:off+p_dm._influsize,off:off+p_dm._influsize] = influ[:,:,0]
    kernconv=np.roll(kernconv,kernconv.shape[0]/2,axis=0)
    kernconv=np.roll(kernconv,kernconv.shape[1]/2,axis=1)
    p_dm._influkernel= kernconv

    return influ


cdef make_tiptilt_dm(Param_dm p_dm,Param_wfs p_wfs, Param_geom p_geom, Param_tel p_tel):
    cdef int dim = max(p_dm._n2-p_dm._n1+1,p_geom._mpupil.shape[0])

    cdef int patchDiam = p_geom.pupdiam+\
    2*np.max(np.abs([p_wfs.xpos,p_wfs.ypos]))*4.848e-6*abs(p_dm.alt/p_tel.diam*p_geom.pupdiam)

    cdef nzer=2
    cdef np.ndarray[ndim=3,dtype=np.float32_t] influ =make_zernike(nzer+1, dim,
            patchDiam, p_geom.cent-p_dm._n1+1, p_geom.cent-p_dm._n1+1, 1)[:,:,1:]

    #normalization factor: one unit of tilt gives 1 arcsec:
    cdef float current = influ[dim/2-1,dim/2-1,0]-influ[dim/2-2,dim/2-2,0]
    cdef float fact = p_dm.unitpervolt*p_tel.diam/p_geom.pupdiam*4.848/current

    influ=influ*fact
    p_dm._ntotact=influ.shape[2]
    p_dm._influ=influ

    return influ


cdef make_kl_dm(Param_dm p_dm, Param_wfs p_wfs,Param_geom p_geom, Param_tel p_tel):
    cdef int dim=p_geom._mpupil.shape[0]


    cdef long patchDiam=long( p_geom.pupdiam+2*max(abs(p_wfs.xpos),abs(p_wfs.ypos))*4.848e-6*\
                            abs(p_dm.alt)/p_geom.pupdiam )

    #TODO klbas
    #p_dm._klbas=make_klbas(p_dm.nkl,p_tel.cobs,patchDiam,"kolmo")
    p_dm._i1=np.zeros((p_dm.nkl),dtype=np.int32)+int((dim-patchDiam)/2.)
    p_dm._j1=np.zeros((p_dm.nkl),dtype=np.int32)+int((dim-patchDiam)/2.)
    p_dm._ntotact=p_dm.nkl


cdef make_zernike(int nzer,int size,int diameter, float xc=-1, float yc=-1, int ext=0):

    cdef int zn,i
    cdef int m=0
    cdef int n=0

    if(xc==-1):
        xc=size/2
    if(yc==-1):
        yc=size/2

    cdef float radius=(diameter+1.)/2.
    cdef np.ndarray[ndim=2,dtype=np.float32_t] zr=mkP.dist(size,xc,yc).astype(np.float32)/radius
    cdef np.ndarray[ndim=3,dtype=np.float32_t] zmask=np.zeros((zr.shape[0],zr.shape[1],nzer),
                                                        dtype=np.float32)
     #(zr<=1).astype(np.float32)
    cdef np.ndarray[ndim=3,dtype=np.float32_t] zmaskmod=np.zeros((zr.shape[0],zr.shape[1],nzer),
                                                        dtype=np.float32)
    #(zr<=1.2).astype(np.float32)

    zmask[:,:,0]=(zr<=1).astype(np.float32)
    zmaskmod[:,:,0]=(zr<=1.2).astype(np.float32)

    for i in range(1,nzer):
        zmask[:,:,i]=zmask[:,:,0]
        zmaskmod[:,:,i]=zmaskmod[:,:,0]

    cdef np.ndarray[ndim=2,dtype=np.float32_t] zrmod=zr*zmaskmod[:,:,0]

    
    zr	= zr*zmask[:,:,0]

    cdef np.ndarray[ndim=2,dtype=np.float32_t] x=np.tile(np.linspace(1,size,size,
                                                    dtype=np.float32), (size,1))
    cdef np.ndarray[ndim=2,dtype=np.float32_t] zteta=np.arctan2(x.T-yc,x-xc).astype(np.float32)

    cdef np.ndarray[ndim=3,dtype=np.float32_t] z=np.zeros((size,size,nzer),dtype=np.float32)



    for zn in range(nzer):
        zernumero(zn+1,&n,&m)

        if ext :
            for i in range( (n-m)/2+1):
                z[:,:,zn]=z[:,:,zn]+(-1.)**i*zrmod**(n-2.*i)*float(np.math.factorial(n-i))/\
                            float(np.math.factorial(i)*np.math.factorial((n+m)/2-i)*
                                  np.math.factorial((n-m)/2-i))
        else:
            for i in range( (n-m)/2+1):
                z[:,:,zn]=z[:,:,zn]+(-1.)**i*zr**(n-2.*i)*float(np.math.factorial(n-i))/\
                            float(np.math.factorial(i)*np.math.factorial((n+m)/2-i)*
                                  np.math.factorial((n-m)/2-i))

        if((zn+1)%2==1):
            if(m==0):
                z[:,:,zn]=z[:,:,zn]*np.sqrt(n+1.)
            else:
                z[:,:,zn]=z[:,:,zn]*np.sqrt(2.*(n+1))*np.sin(m*zteta)
        else:
            if(m==0):
                z[:,:,zn]=z[:,:,zn]*np.sqrt(n+1.)
            else:
                z[:,:,zn]=z[:,:,zn]*np.sqrt(2.*(n+1))*np.cos(m*zteta)


    if(ext): 
        return z*zmaskmod#np.tile(zmaskmod,(nzer,1,1))
    else:
        # nezr?
        return z*zmask#np.tile(zmask,(nzer,1,1))




cdef zernumero(int zn, int *rd, int *an):
    """
    Returns the radial degree and the azimuthal number of zernike
    number zn, according to Noll numbering (Noll, JOSA, 1976)
    """
    cdef int j,m,n
    j=0
    for n in range(101):
        for m in range(n+1):
            if( (n-m)%2 == 0):
                j=j+1
                if(j==zn):
                    rd[0]=n
                    an[0]=m
                    return
                if(m!=0):
                    j=j+1
                    if(j==zn):
                        rd[0]=n
                        an[0]=m
                        return




cdef comp_dmgeom(Param_dm dm, Param_geom geom):
    cdef int smallsize = dm._influsize
    cdef int nact=dm._ntotact
    cdef int dims = int(dm._n2-dm._n1+1)
    cdef int dim = geom._mpupil.shape[0]
    cdef int offs

    if(dims < dim):
        offs=(dim - dims)/2
    else:
        offs = 0
        dim  = dims

    cdef np.ndarray[ndim=2,dtype=np.float32_t] mapactu = np.ones((dim,dim),dtype=np.float32)

    if(offs>0):
        mapactu[:offs,:]=0
        mapactu[:,:offs]=0
        mapactu[mapactu.shape[0]-offs:,:]=0
        mapactu[:,mapactu.shape[1]-offs:]=0

    cdef np.ndarray[ndim=2,dtype=np.int32_t] indgen= np.tile(np.arange(smallsize,dtype=np.int32),(smallsize,1))

    cdef np.ndarray[ndim=3,dtype=np.int32_t] tmpx=np.tile(indgen,(nact,1,1)).T
    cdef np.ndarray[ndim=3,dtype=np.int32_t] tmpy=np.tile(indgen.T,(nact,1,1)).T

    tmpx+=offs+np.tile(dm._i1,(smallsize,smallsize,1))
    tmpy+=offs+np.tile(dm._j1,(smallsize,smallsize,1))

    cdef np.ndarray[ndim=3,dtype=np.int32_t] tmp=tmpx+dim*tmpy


    tmp[tmpx<0]=-10
    tmp[tmpy<0]=-10
    tmp[tmpx>dims-1]=-10
    tmp[tmpy>dims-1]=-10
   
    cdef np.ndarray[ndim=1,dtype=np.int32_t] itmps= np.argsort(tmp.flatten("F")).astype(np.int32)
    cdef np.ndarray[ndim=1,dtype=np.int32_t] tmps=np.sort(tmp.flatten("F")).astype(np.int32)
    itmps=itmps[itmps>-1]

    cdef np.ndarray[ndim=2,dtype=np.int32_t] istart,npts
    istart=np.zeros((dim,dim),dtype=np.int32)
    npts=np.zeros((dim,dim),dtype=np.int32) 


    cdef int cpt,ref
    cpt=1
    ref=1
    
    for i in range(dim*dim):
        if( offs!=0 and mapactu[[i]]==0):
            npts[i]=0
        else:
            while(tmps[cpt]==i-1 and cpt<tmp.size-1 ):
                cpt+=1
            npts.itemset(i,cpt-ref)
            istart.itemset(i,ref)
            ref=cpt

    dm._influpos = itmps[:np.sum(npts)].astype(np.int64)
    dm._ninflu = npts.astype(np.int32)
    dm._influstart= istart.astype(np.int32)

    dm._i1+=offs
    dm._j1+=offs
        


#################################################
# P-Class Dms
#################################################
cdef class Dms:
    def __cinit__(self,long ndm):
        cdef carma_context *context=carma_context.instance()
        self.dms= new sutra_dms(ndm)

    cdef add_dm(self, bytes type_dm, float alt, long dim, long ninflu, long influsize, long ninflupos, long npts, float push4imat, int device=-1):

        cdef carma_context *context=carma_context.instance()

        if(device>-1):
            context.set_activeDevice(device,1)
        else:
            device=context.get_activeDevice()

        self.dms.add_dm(context,type_dm, alt, dim,ninflu,
                        influsize, ninflupos ,npts ,push4imat, device)


    cdef remove_dm(self,bytes type_dm,float alt):
        self.dms.remove_dm(type_dm,alt)

    def resetdm(self, bytes type_dm, float alt):
        cdef int inddm = self.dms.get_inddm(type_dm,alt)
        if(inddm<0):
            raise StandardError("Unknown error")
        self.dms.d_dms[inddm].reset_shape()


    def oneactu( self, bytes type_dm, float alt, int nactu, float ampli):
        cdef int inddm = self.dms.get_inddm(type_dm,alt)
        if(inddm<0):
            raise StandardError("Unknown error")
        self.dms.d_dms[inddm].comp_oneactu(nactu,ampli)


    cdef load_pzt(self, float alt,
        np.ndarray[ndim=3,dtype=np.float32_t] influ,
        np.ndarray[ndim=1,dtype=np.int32_t] influpos,
        np.ndarray[ndim=2,dtype=np.int32_t] npoints,
        np.ndarray[ndim=2,dtype=np.int32_t] istart,
        np.ndarray[ndim=1,dtype=np.int32_t] xoff,
        np.ndarray[ndim=1,dtype=np.int32_t] yoff):#,
        #np.ndarray[dtype=np.float32_t] kern):

        cdef np.ndarray[dtype=np.float32_t] influ_F=influ.flatten("F")
        cdef np.ndarray[dtype=np.int32_t] npoints_F=npoints.flatten("F")
        cdef np.ndarray[dtype=np.int32_t] istart_F=istart.flatten("F")

        cdef int inddm=self.dms.get_inddm("pzt",alt)
        if(inddm<0):
            err="unknown error whith load_pzt\nDM (pzt"+str(alt)+") doesn't exist"
            raise ValueError(err)

        #TODO check kern array?
        self.dms.d_dms[inddm].pzt_loadarrays(<float*>influ_F.data,
                                              <int*> influpos.data,
                                              <int*>npoints_F.data,
                                              <int*>istart_F.data,
                                              <int*>xoff.data,
                                              <int*>yoff.data,
                                              NULL) #<float*>kern.data)


    #TODO dims of arrays
    cdef load_kl(self,float alt,
        np.ndarray[ndim=1,dtype=np.float32_t] rabas,
        np.ndarray[ndim=1,dtype=np.float32_t] azbas,
        np.ndarray[ndim=1,dtype=np.int32_t] ord,
        np.ndarray[ndim=1,dtype=np.float32_t] cr,
        np.ndarray[ndim=1,dtype=np.float32_t] cp):

        cdef int inddm=self.dms.get_inddm("kl",alt)
        if(inddm<0):
            err="unknown error whith load_kl\nDM (kl"+str(alt)+") doesn't exist"
            raise ValueError(err)


        self.dms.d_dms[inddm].kl_loadarrays(<float*> rabas.data,
                                            <float*>azbas.data,
                                            <int*>ord.data,
                                            <float*>cr.data,
                                            <float*>cp.data)



    cdef load_tt(self,float alt,np.ndarray[ndim=3,dtype=np.float32_t] influ):
        cdef int inddm=self.dms.get_inddm("tt",alt)
        if(inddm<0):
            err="unknown error whith load_tt\nDM (tt"+str(alt)+") doesn't exist"
            raise ValueError(err)
        cdef const long *dims
        dims= self.dms.d_dms[inddm].d_influ.getDims()
        cdef np.ndarray[dtype=np.float32_t] influ_F=influ.flatten("F")
        self.dms.d_dms[inddm].d_influ.host2device(<float*>influ_F.data)


    cdef set_comm(self,bytes type_dm,float alt,
                    np.ndarray[ndim=1,dtype=np.float32_t] comm):

        cdef int inddm=self.dms.get_inddm(type_dm,alt)
        if(inddm<0):
            err="unknown error whith load_kl\nDM (kl"+str(alt)+") doesn't exist"
            raise ValueError(err)

        self.dms.d_dms[inddm].d_influ.host2device(<float*>comm.data)

    cdef shape_dm(self,bytes type_dm,float alt):
        cdef int inddm=self.dms.get_inddm(type_dm,alt)
        if(inddm<0):
            err="unknown error whith load_kl\nDM (kl"+str(alt)+") doesn't exist"
            raise ValueError(err)

        self.dms.d_dms[inddm].comp_shape()



    cdef computeKLbasis(self, bytes type_dm, float alt, 
        np.ndarray[ndim=1,dtype=np.float32_t] xpos, np.ndarray[ndim=1,dtype=np.float32_t] ypos,
        np.ndarray[ndim=1,dtype=np.int32_t] indx_pup, long dim, float norm, float ampli):

        cdef int inddm=self.dms.get_inddm(type_dm,alt)
        if(inddm<0):
            err= "unknown error with computeKLbasis function\nDM("+type_dm+","+\
                str(alt)+") doesn't exist"
            raise ValueError(err)

        self.dms.d_dms[inddm].compute_KLbasis(<float*>xpos.data, <float*>ypos.data, 
                <int*>indx_pup.data, dim, norm, ampli)


    cdef get_KLbasis(self,bytes type_dm, float alt):

        cdef int inddm=self.dms.get_inddm(type_dm,alt)
        if(inddm<0):
            err="unknown error with yoga_getKLbasis function DM("+type_dm+","+alt+") doesnt exists"
            raise ValueError(err)
        cdef const long *dims=self.dms.d_dms[inddm].d_KLbasis.getDims()
        #TODO check ndim (=2?)
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data_F=np.zeros((dims[1],dims[2]),dtype=np.float32)
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data=np.zeros((dims[1],dims[2]),dtype=np.float32)

        self.dms.d_dms[inddm].d_KLbasis.device2host(<float*>data_F.data)
        data=np.reshape(data_F.flatten("F"),(dims[1],dims[2]))
        return data


    def __str__(self):
        info= "DMs object:\n"
        info+= "Contains "+str(self.dms.d_dms.size())+" DMs:\n"
        info+= "DM # | Type  |   Alt   | Nact | Dim\n"
        cdef vector[sutra_dm *].iterator it_dms = self.dms.d_dms.begin()
        cdef vector[type_screen].iterator it_type= self.dms.d_type.begin()
        cdef sutra_dm * dm
        cdef type_screen ts
        cdef int i=1
        while(it_dms !=self.dms.d_dms.end()):
            dm=deref(it_dms)
            ts=deref(it_type)
            info+= "%4d"%i+" | "+"%5s"%<bytes> ts.first+" | "+"%7d"%ts.second+" | "+"%4d"%dm.ninflu+" | "+"%4d"%dm.dim+"\n"
            i=i+1
            inc(it_dms)
            inc(it_type)
        info+= "--------------------------------------------------------"
        return info


cdef compute_klbasis(Dms g_dm,Param_dm p_dm, Param_geom p_geom,Param_atmos p_atmos,Param_tel p_tel):
    cdef int tmp
    if(p_dm.type_dm=="pzt"):
        tmp=(p_geom._ipupil.shape[0]-(p_dm._n2-p_dm._n1+1))/2
        tmp_e0=p_geom._ipupil.shape[0]-tmp
        tmp_e1=p_geom._ipupil.shape[1]-tmp
        pup=p_geom._ipupil[tmp:tmp_e0,tmp:tmp_e1]
        indx_valid=np.where(pup.flatten("F")>0)
        interactp=p_dm._xpos1[1]-p_dm._xpos1[2]
        interactm=p_tel.diam/(p_dm.nact-1)
        p2m=interactm/interactp
        norm=-(p2m*p_tel.diam/(2*p_atmos.r0))**(5./3)
        # TODO
        #g_dm.computeKLbasis("pzt",p_dm.alt, p_dm._xpos,p_dm._ypos,indx_valid,indx_valid.size,norm,1.0)
        # KLbasis=get_KLbasis(g_dm,"pzt",y_dm(ndm).alt)[,::-1]
        KLbasis=np.array([])
    else:
        KLbasis=np.array([])

    return KLbasis

