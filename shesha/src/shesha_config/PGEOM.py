#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Created on 13 juil. 2017

@author: vdeo
'''
import os
try:
    shesha_dir = os.environ['SHESHA_ROOT']
    os.environ["PATH"] += shesha_dir + '/src'
except KeyError as err:
    raise EnvironmentError(
            "Environment variable 'SHESHA_ROOT' must be defined")
try:
    shesha_db = os.environ['SHESHA_DB_ROOT']
except KeyError as err:
    import warnings
    shesha_db = shesha_dir + "/data/"
    warnings.warn(
            "'SHESHA_DB_ROOT' not defined, using default one: " +
            str(shesha_db))
finally:
    shesha_savepath = shesha_db

import numpy as np

from . import config_setter_utils as csu
from . import shesha_constants as const

import shesha_util.make_pupil as makeP
import shesha_util.make_apodizer as makeA

from shesha_config.PTEL import Param_tel


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
        self._phase_ab_M1_m = None  # Phase aberration in the pupil (medium size)
        self._apodizer = None  # apodizer (same size as small pupil)
        self._p1 = 0  # min x,y for valid points in mpupil
        self._p2 = 0  # max x,y for valid points in mpupil
        self._n = 0  # linear size of mpupil
        self._n1 = 0  # min x,y for valid points in ipupil
        self._n2 = 0  # max x,y for valid points in ipupil

    def geom_init(self, p_tel: Param_tel):
        """
            Initialize the system geometry

        :parameters:
            p_tel: (Param_tel) : telescope settings
        """

        # First power of 2 greater than pupdiam
        self.ssize = int(2**np.ceil(np.log2(self.pupdiam) + 1))
        # Using images centered on 1/2 pixels
        self.cent = self.ssize / 2 + 0.5

        self._p1 = int(np.ceil(self.cent - self.pupdiam / 2.))
        self._p2 = int(np.floor(self.cent + self.pupdiam / 2.))

        self.pupdiam = self._p2 - self._p1 + 1

        self._n = self.pupdiam + 4
        self._n1 = self._p1 - 2
        self._n2 = self._p2 + 2

        cent = self.pupdiam / 2. + 0.5

        # Useful pupil
        self._spupil = makeP.make_pupil(
                self.pupdiam, self.pupdiam, p_tel, cent,
                cent).astype(np.float32)

        self._phase_ab_M1 = makeP.make_phase_ab(
                self.pupdiam, self.pupdiam, p_tel,
                self._spupil).astype(np.float32)

        # large pupil (used for image formation)
        self._ipupil = makeP.pad_array(self._spupil,
                                       self.ssize).astype(np.float32)

        # useful pupil + 4 pixels
        self._mpupil = makeP.pad_array(self._spupil,
                                       self._n).astype(np.float32)

        self._phase_ab_M1_m = makeP.pad_array(self._phase_ab_M1,
                                              self._n).astype(np.float32)

        if self.apod:
            if self.apodFile is None or self.apodFile == '':
                apod_filename = shesha_savepath + \
                    "apodizer/SP_HARMONI_I4_C6_N1024.npy"
            self._apodizer = makeA.make_apodizer(
                    self.pupdiam, self.pupdiam,
                    apod_filename.encode(), 180. / 12.).astype(np.float32)
        else:
            self._apodizer = np.ones(self._spupil.shape, dtype=np.int32)

        self.isInit = True

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
