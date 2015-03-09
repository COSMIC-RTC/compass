
import numpy as np

import os.path

import iterkolmo as itK
import matplotlib.pyplot as pl

from cython.operator cimport dereference as deref, preincrement as inc


cdef float dtor = np.pi/180

assert sizeof(int) == sizeof(np.int32_t)
assert sizeof(long) == sizeof(np.int64_t)
assert sizeof(float) == sizeof(np.float32_t)
assert sizeof(double) == sizeof(np.float64_t)

data="./data/"

include "atmos.pyx"
include "geom.pyx"
include "target.pyx"
include "tel.pyx"


def see_atmos_target(int n, atmos atm, target tar, float alt=0, int n_tar=0,float f=1, int log=0):
    """Display the turbulence of the atmos and the image of the target after a call to the function:
    - move_atmos
    - target.atmos_raytrace

    n       -- int      : number of iterations
    atm     -- atmos    : atmos used
    tar     -- target   : target used
    alt     -- float    : altitude of the turbulence to diplay
    n_tar   -- int      : number of the target 
    f       -- float    : fraction of the image to display (centered on the image center)
    """
    f=f/2
    fig, (turbu,image)=pl.subplots(1,2, figsize=(15,10))

    ph=tar.get_image(n_tar,"se")
    s0=max(0,ph.shape[0]*0.5-ph.shape[0]*f)
    e0=min(ph.shape[0],ph.shape[0]*0.5+ph.shape[0]*f)
    s1=max(0,ph.shape[1]*0.5-ph.shape[1]*f)
    e1=min(ph.shape[1],ph.shape[1]*0.5+ph.shape[1]*f)
    
    pl.ion()
    pl.show()
    for i in range(n):
        atm.move_atmos()
        tar.atmos_trace(n_tar, atm)
        turbu.clear()
        im1=turbu.imshow(atm.get_screen(alt),cmap='Blues')
        ph=tar.get_image(n_tar,"se")
        ph=np.roll(ph,ph.shape[0]/2,axis=0)
        ph=np.roll(ph,ph.shape[1]/2,axis=1)
        image.clear()
        if(log==1):
            ph=np.log(ph[s0:e0,s1:e1])
        im2=image.matshow(ph[s0:e0,s1:e1],cmap='Blues_r')
        pl.draw()
