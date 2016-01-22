# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 11:19:01 2016

@author: sevin
"""

cdef class Target_brama(Target): # child constructor must have the same prototype (same number of non-optional arguments)
    def __cinit__(self,naga_context ctxt, Telescope telescope, int ntargets, 
                    np.ndarray[ndim=1,dtype=np.float32_t] xpos,
                    np.ndarray[ndim=1,dtype=np.float32_t] ypos,
                    np.ndarray[ndim=1,dtype=np.float32_t] Lambda,
                    np.ndarray[ndim=1,dtype=np.float32_t] mag,
                    float zerop,
                    np.ndarray[ndim=1,dtype=np.int64_t] size,
                    int Npts,
                    int device=-1
                    ):
        del self.target

        cdef carma_context *context =carma_context.instance()
        self.target=new sutra_target_brama(context, "rtc_brama", telescope.telescope, -1, ntargets,
                    <float*>xpos.data,<float*>ypos.data,
                    <float*>Lambda.data,<float*>mag.data,
                    zerop, <long*>size.data,
                    Npts, context.get_activeDevice())

    def __dealloc__(self):
        pass #del self.target

    cpdef publish(self):
        cdef sutra_target_brama* target = <sutra_target_brama*>(self.target)
        target.publish()
