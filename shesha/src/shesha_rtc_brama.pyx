# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:31:50 2016

@author: sevin
"""

from shesha_param import *
from shesha_param cimport *
from shesha_rtc import *
from shesha_rtc cimport *


cdef class Rtc_brama(Rtc): # child constructor must have the same prototype (same number of non-optional arguments)
    def __cinit__(self, Sensors sensor=None, Target target=None, device=-1):
        del self.rtc
        
        cdef carma_context *context =carma_context.instance()
        self.rtc=new sutra_rtc_brama(context, sensor.sensors, target.target, "rtc_brama")

    def __dealloc__(self):
        pass #del self.rtc

    cpdef publish(self):
        cdef sutra_rtc_brama* rtc = <sutra_rtc_brama*>(self.rtc)
        rtc.publish()
        