#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Created on 13 juil. 2017

@author: vdeo
'''
from . import config_setter_utils as csu


#################################################
# P-Class (parametres) Param_geom
#################################################
class Param_geom:

    def __init__(self):
        """ Private members were initialized yet """
        self.isInit = False
        """linear size of full image (in pixels)."""
        self.__ssize = 0
        """observations zenith angle (in deg)."""
        self.__zenithangle = 0.
        """boolean for apodizer"""
        self.__apod = False
        """ File to load an apodizer from """
        self.apodFile = None
        """linear size of total pupil (in pixels)."""
        self.__pupdiam = 0
        """central point of the simulation."""
        self.__cent = 0.

        # Internals
        self._ipupil = None  # total pupil (include full guard band)
        self._mpupil = None  # medium pupil (part of the guard band)
        self._spupil = None  # small pupil (without guard band)
        self._phase_ab_M1 = None  # Phase aberration in the pupil (small size)
        # Phase aberration in the pupil (medium size)
        self._phase_ab_M1_m = None
        self._apodizer = None  # apodizer (same size as small pupil)
        self._p1 = 0  # min x,y for valid points in mpupil
        self._p2 = 0  # max x,y for valid points in mpupil
        self._n = 0  # linear size of mpupil
        self._n1 = 0  # min x,y for valid points in ipupil
        self._n2 = 0  # max x,y for valid points in ipupil

    def set_ssize(self, s):
        """Set linear size of full image

         :param s: (long) : linear size of full image (in pixels)."""
        self.__ssize = csu.enforce_int(s)

    ssize = property(lambda x: x.__ssize, set_ssize)

    def set_zenithangle(self, z):
        """Set observations zenith angle

         :param z: (float) : observations zenith angle (in deg)."""
        self.__zenithangle = csu.enforce_float(z)

    zenithangle = property(lambda x: x.__zenithangle, set_zenithangle)

    def set_pupdiam(self, p):
        """Set the linear size of total pupil

        :param p: (long) : linear size of total pupil (in pixels)."""
        self.__pupdiam = csu.enforce_int(p)

    pupdiam = property(lambda x: x.__pupdiam, set_pupdiam)

    def set_cent(self, c):
        """Set the central point of the simulation

         :param c: (float) : central point of the simulation."""
        self.__cent = csu.enforce_float(c)

    cent = property(lambda x: x.__cent, set_cent)

    def set_apod(self, a):
        """
            Tells if the apodizer is used
            The apodizer is used if a is not 0
        :param a: (int) boolean for apodizer
        """
        self.__apod = csu.enforce_or_cast_bool(a)

    apod = property(lambda x: x.__apod, set_apod)

    def get_ipupil(self):
        """return the full pupil support"""
        if self._ipupil is None:
            raise ValueError("Geometry is unintialized and ipupil is None")
        return self._ipupil

    def get_mpupil(self):
        """return the padded pupil"""
        if self._mpupil is None:
            raise ValueError("Geometry is unintialized and ipupil is None")
        return self._mpupil

    def get_spupil(self):
        """return the small pupil"""
        if self._spupil is None:
            raise ValueError("Geometry is unintialized and ipupil is None")
        return self._spupil

    def get_n(self):
        """Return the linear size of the medium pupil"""
        return self._n

    def get_n1(self):
        """Return the min(x,y) for valid points for the total pupil"""
        return self._n1

    def get_n2(self):
        """Return the max(x,y) for valid points for the total pupil"""
        return self._n2

    def get_p1(self):
        """Return the min(x,y) for valid points for the medium pupil"""
        return self._p1

    def get_p2(self):
        """Return the max(x,y) for valid points for the medium pupil"""
        return self._p2
