cdef class Param_kl_basis:
    def __cinit__(int dim, long nfunc=500,float cobs=0, long Nr=-1, long Np=-1,
                bytes funct=<bytes>"kolmo", outscl=None):
        #TODO check funct=bytes
        if(Nr==-1):
            Nr=long(5.0*np.sqrt(nfunc))
        if(Np==-1):
            Np=long(5*Nr)

        cdef np.ndarray[ndim=1,dtype=np.float32_t] radp=make_radii(Nr,cobs)
        cdef np.ndarray[ndim=3,dtype=np.float32_t] kers=\
                                    make_kernels(cobs,Nr,radp,funct,outscl)


"""
func make_klbas(nfunc,cobs,dim,nr=,np=,funct=,outscl=)
  /*DOCUMENT make_klbas(nfunc,cobs,nr=,np=,funct=,outscl=)
  SEE ALSO : 
   */
{
  if (cobs == []) cobs = 0;
  if (nfunc == []) nfunc = 500L;
  
  if (!is_set(nr)) nr = long(5.0f*sqrt(nfunc));
  if (!is_set(np)) np = long(5*nr);
  if (dimsof(funct)==[]) funct="kolmo";
  
  radp = make_radii(nr, cobs);

  kers = make_kernels(cobs, nr, radp,funct=funct,outscl=outscl);
 --------------------------------------- 
  gkl_fcom,kers,cobs,nfunc,evals,nord,npo,ord,rabas;

  azbas = make_azimuth(nord, np);
 
  klbasis = kl_basis_struct();
  klbasis.nr=nr;
  klbasis.np=np;
  klbasis.nfunc=nfunc; 
  klbasis.cobs=cobs;
  klbasis.radp=&radp;
  klbasis.evals=&evals;
  klbasis.nord=nord;
  klbasis.npo=&npo;
  klbasis.ord=&ord;
  klbasis.rabas=&rabas;
  klbasis.azbas=&transpose(azbas);
  
  set_pctr,klbasis, ncp= dim;
  
  return klbasis;
}
"""



cdef make_radii(int nr,float cobs):
    """  This routine generates an nr x np array with np copies of the
  radial coordinate array. Radial coordinate span the range from
  r=cobs to r=1 with successive annuli having equal areas (ie, the
  area between cobs and 1 is divided into nr equal rings, and the
  points are positioned at the half-area mark on each ring). There
  are no points on the border.
    """
    #TODO check type /nr
    cdef float d=(1.-cobs*cobs)/nr

    return np.sqrt(cobs**2+d/16.+d*np.arange(nr,dtype=np.float32))



cdef make_kernels(float cobs, int nr, np.ndarray[ndim=1,dtype=np.float32_t] rad,
                    bytes funct ,outscl=None):
    """This routine generates the kernel used to find the KL modes.
  The  kernel constructed here should be simply a discretization
  of the continuous kernel. It needs rescaling before it is treated
  as a matrix for finding  the eigen-values. The outter scale
  should be in units of the diameter of the telescope.
    """

    cdef int nth=5*nr
    cdef np.ndarray[ndim=3,dtype=np.float] kers=np.zeros((nr, nr, nth),dtype=np.float32)
    cdef np.ndarray[ndim=1,dtype=np.float] cth =np.cos(np.arange(nth)*2.*np.pi/nth)
    cdef float dth = 2.*np.pi/nth
    cdef float fnorm = -1./(2*np.pi*(1.-cobs*cobs))*0.5
    #the 0.5 is to give  the r^2 kernel, not
    #the r kernel
    cdef int i,j
    cdef float te

    for i in range(nr):
        for j in range(i):
            te = 0.5*np.sqrt(rad[i]**2+rad[j]**2-(2*rad[i]*rad[j])*cth)
            #te in units of the diameter, not the radius
            if (funct=="kolmo"):  te = kolstf(te)
            if (funct=="karman"): te = karmanstf(te,outscl=outscl)
            if ((funct!="kolmo") and (funct!="karman")):
                raise ValueError("The statistics is not known !")

            kelt =  fnorm * dth * (np.fft.ifft(te,-1))
            kers [i, j,:]  = kelt
            kers [j, i,:] = kelt

    return kers



cdef kolstf(np.ndarray[ndim=1,dtype=np.float32_t] dvec):
    """This routine returns the kolmogorov phase variance at spatial
  dimension (inverse of the spatial frequency) dvec"""
    return 6.88 * dvec**(5./3.)


cdef karmanstf(float dvec, int outscl=3):
    """This routine returns the Von Karman phase variance at spatial
  dimension (inverse of the spatial frequency) dvec. Same as kolstf
  but with a correcting factor to account for the outter scale.
  The latter should be in units of telescope diameter
  """
    return  6.88 * dvec**(5./3.)*(1-1.485*(dvec/outscl)**(1./3.)+\
                              5.383*(dvec/outscl)**(2)-6.281*\
                              (dvec/outscl)**(7./3.))





cdef gkl_fcom(np.ndarray[ndim=3,dtype=np.float32_t] kers, float cobs,float nf):
    """This routine does the work : finding the eigenvalues and
  corresponding eigenvectors. Sort them and select the right
  one. It returns the KL modes : in polar coordinates : rabas
  as well as the associated variance : evals. It also returns
  a bunch of indices used to recover the modes in cartesian
  coordinates (nord, npo and ord).
  """

    cdef int nr=kers.shape[0]
    cdef int nt=kers.shape[3]
    cdef int nxt=0
    cdef float fktom=(1.-cobs**2)/nr
    cdef float fvtos=np.sqrt(2*nr)
    cdef np.ndarray[ndim=2,dtype=np.float32_t] evs=np.zeros((nr,nt),dtype=np.float32)
    cdef np.ndarray[ndim=2,dtype=np.float32_t] s=piston_orth(nr)

    cdef np.ndarray[ndim=2,dtype=np.float32_t] b1,v0,vt,v1,vs
    cdef np.ndarray[ndim=2,dtype=np.int32_t] egtmxn
    cdef np.ndarray[ndim=1,dtype=np.float32_t] newev

    cdef float do

    #ff isnt used - the normalisation for
    #the eigenvectors is straightforward:
    #integral of surface^2 divided by area = 1,
    #and the cos^2 term gives a factor
    #half, so multiply zero order by
    #sqrt(n) and the rest by sqrt (2n)

    #zero order is a special case...
    #need to deflate to eliminate infinite eigenvalue - actually want
    #evals/evecs of zom - b where b is big and negative

    b1=np.dot(np.dot(s.T,kers[:,:,0]),s)[0:nr-1,0:nr-1]

    v0,newev,vt=np.linalg.svd(fktom*b1)

    v1=np.zeros((nr,nr),dtype=np.float32)
    v1[0:nr-1,0:nr-1]=v0
    v1[nr-1,nr-1]=1

    vs=np.dot(s,v1)
    evs[:newev.shape[0],nxt]=newev
    kers[:,:,nxt]=np.sqrt(nr)*vs

    nxt=1
    do=nf-1
    while(do<nf):
        vs,newev,vt=np.linalg.svd(fktom*kers[:,:,nxt])
        evs[:newev.shape[0],nxt]=newev
        kers[:,:,nxt]=np.sqrt(2*nr)*vs
        egtmxn = (evs[:,:nxt]>np.max(newev)).astype(np.int32)
        do=2*np.sum(egtmxn)-sum(egtmxn[:,0])
        nxt+=1

    kers=kers[:,:,:nxt]
    evs=evs[:,:nxt].flatten()

    cdef np.ndarray[ndim=1,dtype=np.int32_t] a=np.argsort(-evs)[:nf]
    #every eigenvalue occurs twice except
    #those for the zeroth order mode. This
    #could be done without the loops, but
    #it isn't the stricking point anyway...
    no=1
    ni=1
    cdef np.ndarray[ndim=1,dtype=np.int32_t]oind=np.zeros(nf+1,dtype=np.int32)
    while(no<nf+1):
        if(a[ni]<nr+1):
            oind[no]=a[ni]
            no+=1
        else:
            oind[no]=a[ni]
            oind[no+1]=a[ni]
            no+=2
        ni+=1

    oind=oind[:nf]
    cdef np.ndarray[ndim=1,dtype=np.float32_t] tord=(oind-1)/nr+1
    cdef np.ndarray[ndim=1,dtype=np.int32_t] odd=np.arange(nf,dtype=np.int32)%2
    cdef np.ndarray[ndim=1,dtype=np.int32_t] pio=(oind-1)%nr+1


    cdef np.ndarray[ndim=1,dtype=np.float32_t ] evals=evs[oind]
    cdef np.ndarray[ndim=1,dtype=np.int32_t] ord=2*(tord-1)-np.floor((tord>1)*odd) +1
    cdef int nord=np.max(ord)
    cdef np.ndarray[ndim=2,dtype=np.float32_t] rabas=np.zeros((nr,nf),dtype=np.float32)
    cdef np.ndarray[ndim=1,dtype=np.int32_t] npo=np.zeros(nord,dtype=np.int32)

    cdef int i
    for i in range(long(nf)):
        npo[ord[i]]=npo[ord[i]]+1
        rabas[:,i]=kers[:,pio[i],tord[i]]

    return evals,nord,npo,ord,rabas


cdef make_azimuth(nord,Np):
    cdef np.ndarray[ndim=2,dtype=np.float32_t] azi= np.zeros((int(nord+1),Np),
                                                    dtype=np.float32)
    cdef np.ndarray[ndim=1, dtype=np.float32_t] th=np.arange(Np)*(2.*np.pi/Np)

    cdef int i
    azi[0,:]=1.
    azi[1,:]=np.cos(1.5*th)

    for i in range(2,nord+1):
        azi[i,:]=np.sin((i/2)*th)

    return azi



"""
func make_azimuth(nord, np)
  /*DOCUMENT piston_orth(nr)
   */
{
  azi = array(float,[2,long(1+nord), np]);
  th = float(indgen(np)-1)*(2.*pi/ np);

  #TODO check this out
  azi (1,) = 1.0;
  for (i = 2; i<=nord;i+=2)  azi (i,) = cos (((i-1)/2+1) * th);
  for (i = 3; i<=nord;i+=2)  azi (i,) = sin (((i-1)/2) * th);
  return azi;
}"""


cdef piston_orth(int Nr):
    cdef np.ndarray[ndim=2,dtype=np.float32_t] s=np.zeros((),dtype=np.float32)
    cdef int i
    cdef float rnm

    for i in range(Nr-1):
        nrm=1./np.sqrt(i*(i+1))
        s[0:i+1,i]=rnm
        s[i+1,i]=-i*rnm
    rnm=1./np.sqrt(Nr)
    s[:Nr-1]=rnm

    return s
