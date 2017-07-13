import numpy as np
cimport numpy as np

import h5py

import hdf5_utils as h5u
import pandas

import os
import iterkolmo as itK

from cython.operator cimport dereference as deref, preincrement as inc
from subprocess import check_output

import shesha

#################################################
# P-Class atmos
#################################################
cdef class Atmos:
    def __cinit__(self):
        self.context = None

    def __dealloc__(self):
        if(self.s_a != NULL):
            del self.s_a

    cdef realinit(self, naga_context ctxt, int nscreens,
                  np.ndarray[ndim=1, dtype=np.float32_t] r0,
                  np.ndarray[ndim=1, dtype=np.int64_t] size,
                  np.ndarray[ndim=1, dtype=np.float32_t] altitude,
                  np.ndarray[ndim=1, dtype=np.float32_t] windspeed,
                  np.ndarray[ndim=1, dtype=np.float32_t] winddir,
                  np.ndarray[ndim=1, dtype=np.float32_t] deltax,
                  np.ndarray[ndim=1, dtype=np.float32_t] deltay,
                  int device):
        """
            Create a sutra_atmos object

        :parameters:
            c: (naga_context): context
            nscreens: (int): number of turbulent layers
            r0: (np.ndarray[ndim=1,dtype=np.float32_t]): global r0
            size: (np.ndarray[ndim=1,dtype=np.int64_t]): screen size of each layer
            altitude: (np.ndarray[ndim=1,dtype=np.float32_t]): altitude of each layer
            windspeed: (np.ndarray[ndim=1,dtype=np.float32_t]): wind speed of each layer
            winddir: (np.ndarray[ndim=1,dtype=np.float32_t]): wind direction of each layer
            deltax: (np.ndarray[ndim=1,dtype=np.float32_t]): x translation speed
            deltay: (np.ndarray[ndim=1,dtype=np.float32_t]): y translation speed
            device: (int): device index
        """

        cdef np.ndarray[ndim= 1, dtype = np.int64_t]size2
        size2 = compute_size2(size)

        self.s_a = new sutra_atmos(ctxt.c, nscreens,
                                   < np.float32_t * > r0.data,
                                   < long * > size.data,
                                   < long * > size2.data,
                                   < np.float32_t * > altitude.data,
                                   < np.float32_t * > windspeed.data,
                                   < np.float32_t * > winddir.data,
                                   < np.float32_t * > deltax.data,
                                   < np.float32_t * > deltay.data,
                                   device)
        self.context = ctxt

    def get_screen(self, float alt):
        """
            Return a numpy array containing the turbulence at a given altitude

        :param alt: (float) :altitude of the screen to get
        """
        cdef carma_obj[float] * screen = self.s_a.d_screens[alt].d_tscreen.d_screen
        self.context.set_activeDevice(screen.getDevice(), 1)
        cdef const long * dims
        dims = screen.getDims()
        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.float32)
        screen.device2host(< float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def disp(self, float alt):
        """
            Display the screen phase at a given altitude

        :param alt: (float) : altitude of the screen to display
        """
        cdef carma_obj[float] * c_phase = self.s_a.d_screens[alt].d_tscreen.d_screen
        cdef const long * dims = c_phase.getDims()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data = np.zeros((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] phase = np.ndarray((dims[1], dims[2]),
                                                                           dtype=np.float32)

        c_phase.device2host(< float * > data.data)
        phase = np.reshape(data.flatten("F"), (dims[1], dims[2]))

        pl.ion()
        pl.clf()
        pl.imshow(phase, cmap="Blues")
        pl.show()

    def add_screen(self, long size, float amplitude, float altitude,
                   float windspeed, float winddir, float deltax, float deltay, int device):
        """
            Add a screen to the atmos object.

        :parameters:
            size: (float) : dimension of the screen (size x size)
            amplitude: (float) : frac
            altitude: (float) : altitude of the screen in meters
            windspeed: (float) : windspeed of the screen [m/s]
            winddir: (float) : wind direction (deg)
            deltax: (float) : extrude deltax pixels in the x-direction at each iteration
            deltay: (float) : extrude deltay pixels in the y-direction at each iteration
            device: (int) : device number
        """
        cdef long size2 = compute_size2(np.array([size], dtype=np.int64))[0]

        if(self.s_a.d_screens.find(altitude) != self.s_a.d_screens.end()):
            print("There is already a screen at this altitude")
            print("No screen created")
            return

        cdef sutra_tscreen * screen = new sutra_tscreen(self.s_a.current_context, size, size2, amplitude, altitude, windspeed, winddir, deltax, deltay, device)

        cdef pair[float, sutra_tscreen *] p
        p.first, p.second = altitude, screen
        self.s_a.d_screens.insert(p)
        self.s_a.nscreens += 1

    def del_screen(self, float alt):
        """
            Delete a screen from the atmos object

        :param alt: (float) : altitude of the screen to delete
        """
        if(self.s_a.d_screens.find(alt) == self.s_a.d_screens.end()):
            print("No screen at this altitude")
            print("No screen deleted")
            return
        self.s_a.nscreens -= 1
        self.s_a.d_screens.erase(alt)

    def list_alt(self):
        """Display the list of the screens altitude"""

        cdef map[float, sutra_tscreen *].iterator it
        cdef int i = 0
        cdef np.ndarray alt = np.zeros(self.s_a.nscreens, dtype=np.float32)
        it = self.s_a.d_screens.begin()

        while it != self.s_a.d_screens.end():
            alt[i] = deref(it).first
            inc(it)
            i += 1
        print(alt)

    def move_atmos(self):
        """
            Move the turbulence in the atmos screen following previous loaded
            parameters such as windspeed and wind direction
        """

        self.s_a.move_atmos()

    def __str__(self):
        cdef map[float, sutra_tscreen *].iterator it = self.s_a.d_screens.begin()
        cdef sutra_tscreen * screen
        cdef int i = 1
        info = "Atmos obect:\n"
        info += "Contains " + \
            str(self.s_a.nscreens) + " turbulent screen(s):\n"
        info += "Screen # | alt.(m) | speed (m/s) | dir.(deg) | r0 (pix) | deltax | deltay\n"
        while it != self.s_a.d_screens.end():
            screen = deref(it).second
            info += "%8d" % i + " | " + "%7.4f" % screen.altitude + " | " + "%11.4f" % screen.windspeed + \
                " | " + "%9.4f" % screen.winddir + " | " + "%8.4f" % screen.amplitude ** -(6. / 5.) + \
                " | " + "%6.4f" % screen.deltax + \
                " | " + "%6.4f" % screen.deltay + "\n"
            i = i + 1
            inc(it)

        info += "--------------------------------------------------------"

        return info
