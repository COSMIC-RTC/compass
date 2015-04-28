import cython
# #cython: profile=True
import numpy as np

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

data="./data/"

include "atmos.pyx"
include "geom.pyx"
include "target.pyx"
include "tel.pyx"
include "wfs.pyx"
include "sensors.pyx"
include "loop.pyx"

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

    start=time.time()
    for i in range(n):
        atm.move_atmos()
        tar.atmos_trace(n_tar, atm)
        wfs.sensors_trace(0,"atmos",atm,0)
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
    f=f/2
    fig, (turbu,image,sh)=pl.subplots(1,3, figsize=(15,10))

    ph=tar.get_image(n_tar,"se")
    s0=max(0,ph.shape[0]*0.5-ph.shape[0]*f)
    e0=min(ph.shape[0],ph.shape[0]*0.5+ph.shape[0]*f)
    s1=max(0,ph.shape[1]*0.5-ph.shape[1]*f)
    e1=min(ph.shape[1],ph.shape[1]*0.5+ph.shape[1]*f)
   
    if(wfs.get_rank(0)==0):
        pl.ion()
        pl.show()
    cdef double start,end, t1,t2,t3,t4,t5,t6

    start=time.time()
    for i in range(n):
        if(wfs.get_rank(0)==0):
            atm.move_atmos()
            tar.atmos_trace(n_tar, atm)
            wfs.sensors_trace(0,"atmos",atm,0)
        wfs.Bcast_dscreen()
        wfs.sensors_compimg(0)
        wfs.gather_bincube(comm,0)
        if(wfs.get_rank(0)==0):
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
            wfs.sensors_trace(0,"atmos",atm,0)
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
            wfs.sensors_trace(0,"atmos",atm,0)
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
            wfs.sensors_trace(0,"atmos",atm,0)
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
        #return x,y
        return y,x
    else :
        x =np.tile( (np.arange(dim1)+1),(dim2,1))
        y =np.tile( (np.arange(dim2)+1),(dim1,1)).T
        return x,y


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

