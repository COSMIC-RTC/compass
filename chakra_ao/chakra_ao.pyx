import cython
# #cython: profile=True
import numpy as np
cimport numpy as np

import os.path

import iterkolmo as itK
import make_pupil as mkP
import matplotlib.pyplot as pl

from cython.operator cimport dereference as deref, preincrement as inc

import time

cdef float dtor = np.pi/180
cdef long RASC = 180.*3600./np.pi

assert sizeof(int) == sizeof(np.int32_t)
assert sizeof(long) == sizeof(np.int64_t)
assert sizeof(float) == sizeof(np.float32_t)
assert sizeof(double) == sizeof(np.float64_t)


#chakra_ao_dir=os.environ.get('YOGA_AO_DIR')
chakra_ao_dir="/home/ndoucet/workspace/compass/trunk/chakra_ao"
print chakra_ao_dir
chakra_ao_savepath=chakra_ao_dir+"/data/"
print chakra_ao_savepath
#data="./data/"

include "atmos.pyx"
include "geom.pyx"
include "target.pyx"
include "tel.pyx"
include "wfs.pyx"
include "sensors.pyx"
include "loop.pyx"
include "rtc.pyx"
include "dms.pyx"
include "centroider.pyx"
include "controller.pyx"
include "kl.pyx"

def see_atmos_target_disp(int n, Atmos atm, Target tar,Sensors wfs, float alt=0, int n_tar=0,float f=1, int log=0):
    """Display the turbulence of the atmos and the image of the target after a call to the function:
    - move_atmos
    - target.atmos_raytrace
    no mpi involved

    n       -- int      : number of iterations
    atm     -- Atmos    : Atmos used
    tar     -- Target   : Target used
    alt     -- float    : altitude of the turbulence to diplay
    n_tar   -- int      : number of the target 
    f       -- float    : fraction of the image to display (centered on the image center)
    """
    f=f/2
    fig, (turbu,image,sh)=pl.subplots(1,3, figsize=(15,10))

    ph=tar.get_image(n_tar,"se")
    s0=max(0,ph.shape[0]*0.5-ph.shape[0]*f)
    e0=min(ph.shape[0],ph.shape[0]*0.5+ph.shape[0]*f)
    s1=max(0,ph.shape[1]*0.5-ph.shape[1]*f)
    e1=min(ph.shape[1],ph.shape[1]*0.5+ph.shape[1]*f)

   
    pl.ion()
    pl.show()
    cdef double start,end, t1,t2,t3,t4,t5,t6

    screen=atm.get_screen(alt)
    im1=turbu.imshow(screen,cmap='Blues')
    ph=tar.get_image(n_tar,"se")
    im2=image.matshow(ph[s0:e0,s1:e1],cmap='Blues_r')
    ph=np.roll(ph,ph.shape[0]/2,axis=0)
    ph=np.roll(ph,ph.shape[1]/2,axis=1)
    shak=wfs._get_binimg(0)
    im3=sh.matshow(shak,cmap='Blues_r')
    pl.draw()


    start=time.time()
    for i in range(n):
        atm.move_atmos()
        tar.atmos_trace(n_tar, atm)
        wfs.sensors_trace(0,"atmos",atm)
        wfs.sensors_compimg(0)
        shak=wfs._get_binimg(0)
        turbu.clear()
        screen=atm.get_screen(alt)
        im1=turbu.imshow(screen,cmap='Blues')
        ph=tar.get_image(n_tar,"se")
        ph=np.roll(ph,ph.shape[0]/2,axis=0)
        ph=np.roll(ph,ph.shape[1]/2,axis=1)
        image.clear()
        if(log==1):
            ph=np.log(ph[s0:e0,s1:e1])
        im2=image.matshow(ph[s0:e0,s1:e1],cmap='Blues_r')
        sh.clear()
        im3=sh.matshow(shak,cmap='Blues_r')
        pl.draw()
    end=time.time()
    print "time:",end-start


def see_atmos_target_disp_mpi(int n, Atmos atm, Target tar,Sensors wfs, MPI.Intracomm comm, float alt=0, int n_tar=0,float f=1, int log=0):
    """Display the turbulence of the atmos and the image of the target after a call to the function:
    - move_atmos
    - target.atmos_raytrace
    need mpi communicator

    n       -- int      : number of iterations
    atm     -- Atmos    : Atmos used
    tar     -- Target   : Target used
    wfs     -- Sensors  : Sensor used
    comm    --  MPI.Intracomm : mpi communicator
    alt     -- float    : altitude of the turbulence to diplay
    n_tar   -- int      : number of the target 
    f       -- float    : fraction of the image to display (centered on the image center)
    """
    if(wfs.get_rank(0)==0):
        f=f/2
        fig, (turbu,image,sh)=pl.subplots(1,3, figsize=(15,10))

        ph=tar.get_image(n_tar,"se")
        s0=max(0,ph.shape[0]*0.5-ph.shape[0]*f)
        e0=min(ph.shape[0],ph.shape[0]*0.5+ph.shape[0]*f)
        s1=max(0,ph.shape[1]*0.5-ph.shape[1]*f)
        e1=min(ph.shape[1],ph.shape[1]*0.5+ph.shape[1]*f)
  
    pl.ion()
    pl.show()
    cdef double start,end, t1,t2,t3,t4,t5,t6

    start=time.time()
    for i in range(n):
        if(wfs.get_rank(0)==0):
            print i
            atm.move_atmos()
            tar.atmos_trace(n_tar, atm)
            wfs.sensors_trace(0,"atmos",atm)
        wfs.Bcast_dscreen()
        print "Bcast",wfs.get_rank(0),i
        wfs.sensors_compimg(0)
        print "gather",wfs.get_rank(0),i
        wfs.gather_bincube(comm,0)
        if(wfs.get_rank(0)==0):
            print wfs.get_rank(0),"disp"
            turbu.clear()
            screen=atm.get_screen(alt)
            im1=turbu.imshow(screen,cmap='Blues')

            ph=tar.get_image(n_tar,"se")
            ph=np.roll(ph,ph.shape[0]/2,axis=0)
            ph=np.roll(ph,ph.shape[1]/2,axis=1)
            image.clear()
            if(log==1):
                ph=np.log(ph[s0:e0,s1:e1])
            im2=image.matshow(ph[s0:e0,s1:e1],cmap='Blues_r')

            sh.clear()
            shak=wfs._get_binimg(0)
            im3=sh.matshow(shak,cmap='Blues_r')
            pl.draw()

            sh_file="dbg/shak_"+str(i)+"_np_"+str(comm.Get_size())+".npy"
            im_file="dbg/imag_"+str(i)+"_np_"+str(comm.Get_size())+".npy"
            np.save(sh_file,shak)
            np.save(im_file,ph)

    end=time.time()
    print comm.Get_rank(),"time:",end-start

def see_atmos_target_mpi_cu(int n, Atmos atm, Target tar,Sensors wfs, MPI.Intracomm comm, float alt=0, int n_tar=0,float f=1, int log=0):
    """Compute the turbulence of the atmos and the image of the target after a call to the function:
    - move_atmos
    - target.atmos_raytrace
    need mpi communicator
    use mpi cuda_aware

    n       -- int      : number of iterations
    atm     -- Atmos    : Atmos used
    tar     -- Target   : Target used
    wfs     -- Sensors  : Sensor used
    comm    --  MPI.Intracomm : mpi communicator
    alt     -- float    : altitude of the turbulence to diplay
    n_tar   -- int      : number of the target 
    f       -- float    : fraction of the image to display (centered on the image center)
    """

    cdef double start,end, t1,t2,t3,t4,t5,t6

    start=time.time()
    for i in range(n):
        if(wfs.get_rank(0)==0):
            atm.move_atmos()
            tar.atmos_trace(n_tar, atm)
            wfs.sensors_trace(0,"atmos",atm)
        wfs.Bcast_dscreen_cuda_aware()
        wfs.sensors_compimg(0)
        wfs.gather_bincube_cuda_aware(comm,0)
    end=time.time()
    print comm.Get_rank(),"time:",end-start

def see_atmos_target_mpi(int n, Atmos atm, Target tar,Sensors wfs, MPI.Intracomm comm, float alt=0, int n_tar=0,float f=1, int log=0):
    """Compute the turbulence of the atmos and the image of the target after a call to the function:
    - move_atmos
    - target.atmos_raytrace
    need mpi communicator

    n       -- int      : number of iterations
    atm     -- Atmos    : Atmos used
    tar     -- Target   : Target used
    wfs     -- Sensors  : Sensor used
    comm    --  MPI.Intracomm : mpi communicator
    alt     -- float    : altitude of the turbulence to diplay
    n_tar   -- int      : number of the target 
    f       -- float    : fraction of the image to display (centered on the image center)
    """

    cdef double start,end, t1,t2,t3,t4,t5,t6

    start=time.time()
    for i in range(n):
        if(wfs.get_rank(0)==0):
            atm.move_atmos()
            tar.atmos_trace(n_tar, atm)
            wfs.sensors_trace(0,"atmos",atm)
        wfs.Bcast_dscreen()
        wfs.sensors_compimg(0)
        wfs.gather_bincube(comm,0)
    end=time.time()
    print comm.Get_rank(),"time:",end-start


def see_atmos_target(int n, Atmos atm, Target tar,Sensors wfs, float alt=0, int n_tar=0,float f=1, int log=0):
    """Compute the turbulence of the atmos and the image of the target after a call to the function:
    - move_atmos
    - target.atmos_raytrace
    no mpi involved

    n       -- int      : number of iterations
    atm     -- Atmos    : Atmos used
    tar     -- Target   : Target used
    wfs     -- Sensors  : Sensor used
    alt     -- float    : altitude of the turbulence to diplay
    n_tar   -- int      : number of the target 
    f       -- float    : fraction of the image to display (centered on the image center)
    """

    cdef double start,end, t1,t2,t3,t4,t5,t6

    start=time.time()
    for i in range(n):
        if(wfs.get_rank(0)==0):
            atm.move_atmos()
            tar.atmos_trace(n_tar, atm)
            wfs.sensors_trace(0,"atmos",atm)
        wfs.sensors_compimg(0)
    end=time.time()
    print "time:",end-start


cdef bin2d(np.ndarray data_in, int binfact):
    """
   Returns the input 2D array "array", binned with the binning factor "binfact".
   The input array X and/or Y dimensions needs not to be a multiple of
   "binfact"; The final/edge pixels are in effect replicated if needed.
   This routine prepares the parameters and calls the C routine _bin2d.
   The input array can be of type long, float or double.
   Last modified: Dec 15, 2003.
   Author: F.Rigaut
   SEE ALSO: _bin2d
 */
    """
    if(binfact<1):
        raise ValueError("binfact has to be >= 1")

    cdef int nx,ny,fx,fy
    nx=data_in.shape[0]
    ny=data_in.shape[1]
    fx=int(np.ceil(nx/float(binfact)))
    fy=int(np.ceil(ny/float(binfact)))

    

    cdef np.ndarray data_out=np.zeros((fx,fy),dtype=data_in.dtype) 

    cdef int i,j,i1,i2,j1,j2

    for i1 in range(fx):
        for j1 in range(fy):
            for i2 in range(binfact):
                for j2 in range(binfact):
                    i = i1*binfact+i2
                    j = j1*binfact+j2
                    if(i>=nx):
                        i=nx-1
                    if(j>=ny):
                        j=ny-1
                    data_out[i1,j1]+=data_in[i,j]

    return data_out



#cdef  indices(np.ndarray x, np.ndarray y,
#              int dim1, int dim2=-1):
def  indices(int dim1, int dim2=-1):
    """DOCUMENT indices(dim)
  Return a dimxdimx2 array. First plane is the X indices of the pixels
  in the dimxdim array. Second plane contains the Y indices.
  Inspired by the Python scipy routine of the same name.
  New (June 12 2002): dim can either be :
    - a single number N (e.g. 128) in which case the returned array are
      square (NxN)
    - a Yorick array size, e.g. [#dimension,N1,N2], in which case
      the returned array are N1xN2
    - a vector [N1,N2], same result as previous case
  F.Rigaut 2002/04/03
  SEE ALSO: span
    """


    if (dim2<0):
        y =np.tile( (np.arange(dim1)+1),(dim1,1))
        x =np.copy(y.T) 
        return y,x
    else :
        x =np.tile( (np.arange(dim1)+1),(dim2,1))
        y =np.tile( (np.arange(dim2)+1),(dim1,1)).T
        return y,x


cdef fft_goodsize(long s):
    """find best size for a fft from size s"
    long s
    """
    return 2**(long(np.log2(s))+1)



cdef makegaussian(int size, float fwhm, int xc=-1, int yc=-1, int norm=0):
    """makegaussian(size,fwhm,xc,yc)
        Returns a centered gaussian of specified size and fwhm.
        norm returns normalized 2d gaussian
    """
    cdef np.ndarray tmp
    tmp = np.exp(-(mkP.dist(size,xc,yc)/(fwhm/1.66))**2.)
    if (norm>0):
        tmp = tmp/(fwhm**2.*1.140075)
    return tmp




cdef make_apodizer(int dim, int pupd, bytes filename, float angle):

    print "Opening apodizer"
    cdef np.ndarray pup=np.load(filename)
    cdef int A=pup.shape[0]

    if(A>dim):
        raise ValueError("Apodizer dimensions must be smaller.") 
    
    if (A != pupd):
        #use misc.imresize (with bilinear)
        print "TODO pup=bilinear(pup,pupd,pupd)"
    
    if (angle != 0):
        #use ndimage.interpolation.rotate
        print "TODO pup=rotate2(pup,angle)"
    
    reg=np.where(mkP.dist(pupd)>pupd/2.)
    pup[reg]=0.
    
    cdef np.ndarray[ndim=2,dtype=np.float32_t] pupf= np.zeros((dim,dim),dtype=np.float32)
    
    if (dim != pupd):
        if ((dim-pupd)%2 != 0):
            pupf[(dim-pupd+1)/2:(dim+pupd+1)/2,(dim-pupd+1)/2:(dim+pupd+1)/2]=pup
        
        else:
            pupf[(dim-pupd)/2:(dim+pupd)/2,(dim-pupd)/2:(dim+pupd)/2]=pup
       
    else:
        pupf=pup
    
    pupf=np.abs(pupf).astype(np.float32)
    
    return pupf


'''
cdef rotate2(image,angle, xc=-1,yc=-1, splin=0,outside=0):
    """rotate2(image,angle,xc,yc,splin,outside)

    Rotate the input image. Angle is in degrees, CCW.

    KEYWORDS:
    xc, yc: Center for coordinate transform. Note that this is
    compatible with the center defined by dist(), but is
    offset by 0.5 pixels w.r.t what you read on the yorick graphic
    window. I.e. the center of the bottom- left pixel is (1,1) in this
    function's conventions, not (0.5,0.5).

    splin: use spline2() instead of bilinear() for the interpolation

    outside: value for outliers.
    """

    angle *= np.pi/180.

    x,y = indices(image.shape[0],image.shape[1])

    if (xc<0): xc=np.ceil(image.shape[0]/2.+0.5)
    if (yc<0): yc=np.ceil(image.shape[1]/2.+0.5)

    x-=xc
    y-=yc

    x =  np.cos(angle)*x + np.sin(angle)*y
    y = -np.sin(angle)*x + np.cos(angle)*y

    x +=xc
    y +=yc

#    if (splin!=0) return spline2(image,x,y,outside=outside)
#    return bilinear(image,x,y,outside=outside)
'''




cdef Bcast(carma_obj[float] *obj, int root):
    cdef int i
    cdef int size=<int>obj.getNbElem()
    cdef int ndim=<int>obj.getDims(0)


    cdef float *ptr
    ptr=<float*>malloc(size*sizeof(float))

    obj.device2host(ptr)

    mpi.MPI_Bcast(ptr,size,mpi.MPI_FLOAT,root,mpi.MPI_COMM_WORLD)

    obj.host2device(ptr)


    free(ptr)
