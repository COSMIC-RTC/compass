import numpy as np
cimport numpy as np

import shesha_util.iterkolmo as itK

from cython.operator cimport dereference as deref, preincrement as inc


#################################################
# P-Class atmos
#################################################
cdef class Atmos:
    def __cinit__(self, naga_context context, int nscreens,
                  float r0,
                  float pupixsize,
                  np.ndarray[ndim=1, dtype=np.int64_t] dim_screens,
                  np.ndarray[ndim=1, dtype=np.float32_t] frac,
                  np.ndarray[ndim=1, dtype=np.float32_t] alt,
                  np.ndarray[ndim=1, dtype=np.float32_t] windspeed,
                  np.ndarray[ndim=1, dtype=np.float32_t] winddir,
                  np.ndarray[ndim=1, dtype=np.float32_t] deltax,
                  np.ndarray[ndim=1, dtype=np.float32_t] deltay):
        """atmos_create(naga_context c, int nscreens,
                        float r0,
                        float pupixsize,
                        np.ndarray[ndim=1,dtype=np.int64_t] dim_screens,
                        np.ndarray[ndim=1,dtype=np.float32_t] frac,
                        np.ndarray[ndim=1,dtype=np.float32_t] alt,
                        np.ndarray[ndim=1,dtype=np.float32_t] windspeed,
                        np.ndarray[ndim=1,dtype=np.float32_t] winddir,
                        np.ndarray[ndim=1,dtype=np.float32_t] deltax,
                        np.ndarray[ndim=1,dtype=np.float32_t] deltay,
                        np.ndarray[ndim=1,dtype=np.int64_t] seeds,
                        int clean, dict load)

        Create and initialise an atmos object.

        :parameters:
            c: (naga_context) : context
            nscreens: (float) : number of turbulent layers
            r0: (float) : global r0
            L0: (np.ndarray[ndim=1, dtype=np.float32_t]) : L0
            pupixsize: (float) : pixel size [m]
            dim_screens: (np.ndarray[ndim=1,dtype=np.int64_t]) : screens dimensions
            frac: (np.ndarray[ndim=1,dtype=np.float32_t]) : fraction of r0
            alt: (np.ndarray[ndim=1,dtype=np.float32_t]) : altitudes [m]
            windspeed: (np.ndarray[ndim=1,dtype=np.float32_t]) : wind speed [m/s]
            winddir: (np.ndarray[ndim=1,dtype=np.float32_t]) : wind direction [deg]
            deltax: (np.ndarray[ndim=1,dtype=np.float32_t]) : extrude deltax pixels in the x-direction at each iteration
            deltay: (np.ndarray[ndim=1,dtype=np.float32_t]) : extrude deltay pixels in the y-direction at each iteration
        """

        self.context = context

        # get fraction of r0 for corresponding layer
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] r0_layers
        r0_layers = r0 / (frac ** (3. / 5.) * pupixsize)
        # create atmos object on gpu

        cdef np.ndarray[ndim= 1, dtype = np.int64_t] stencil_size
        stencil_size = itK.stencil_size_array(dim_screens)

        self.s_a = new sutra_atmos(self.context.c, nscreens,
                                   < np.float32_t * > r0_layers.data,
                                   < long * > dim_screens.data,
                                   < long * > stencil_size.data,
                                   < np.float32_t * > alt.data,
                                   < np.float32_t * > windspeed.data,
                                   < np.float32_t * > winddir.data,
                                   < np.float32_t * > deltax.data,
                                   < np.float32_t * > deltay.data,
                                   self.context.get_activeDevice())

    def __dealloc__(self):
        if(self.s_a != NULL):
            del self.s_a

    def refresh_screen(self, float alt):
        """
            Refresh the selected screen by extrusion

        :param alt: (float) :altitude of the screen to get
        """
        if alt not in self.list_alt():
            raise ValueError("No screen at this altitude")

        self.s_a.d_screens[alt].refresh_screen()

    def set_seed(self, float alt, int seed):
        """
            Set the seed of the selected screen RNG

        :param alt: (float) :altitude of the screen to get
        :param seed: (int) :new seed

        """
        if alt not in self.list_alt():
            raise ValueError("No screen at this altitude")

        self.s_a.d_screens[alt].set_seed(seed)

    def get_screen(self, float alt):
        """
            Return a numpy array containing the turbulence at a given altitude

        :param alt: (float) :altitude of the screen to get
        """
        if alt not in self.list_alt():
            raise ValueError("No screen at this altitude")

        cdef carma_obj[float] * screen = self.s_a.d_screens[alt].d_tscreen.d_screen
        cdef const long * dims
        dims = screen.getDims()
        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        status = self.s_a.get_screen(alt, < float * > data_F.data)
        return data_F.T.copy()

    def add_screen(self, long size, float amplitude, float altitude,
                   float windspeed, float winddir, float deltax, float deltay, int device):
        """
            Add a screen to the atmos object.

        :parameters:
            size: (long) : dimension of the screen (size x size)
            amplitude: (float) : frac
            altitude: (float) : altitude of the screen in meters
            windspeed: (float) : windspeed of the screen [m/s]
            winddir: (float) : wind direction (deg)
            deltax: (float) : extrude deltax pixels in the x-direction at each iteration
            deltay: (float) : extrude deltay pixels in the y-direction at each iteration
            device: (int) : device number
        """
        cdef long stencil_size = itK.stencil_size_array(
            np.array([size], dtype=np.int64))[0]

        self.s_a.add_screen(altitude, size, stencil_size, amplitude, windspeed,
                            winddir, deltax, deltay, device)

    def init_screen(self, float alt, np.ndarray[ndim=2, dtype=np.float32_t] A,
                    np.ndarray[ndim=2, dtype=np.float32_t] B,
                    np.ndarray[dtype=np.uint32_t] istx,
                    np.ndarray[dtype=np.uint32_t] isty, int seed):
        """
            TODO: docstring
        """
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] A_F = A.T.copy()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] B_F = B.T.copy()

        self.s_a.init_screen(< float > alt, < float * > (A_F.data),
                              < float * > (B_F.data),
                              < unsigned int * > istx.data,
                              < unsigned int * > isty.data, seed)

    def del_screen(self, float alt):
        """
            Delete a screen from the atmos object

        :param alt: (float) : altitude of the screen to delete
        """
        self.s_a.del_screen(alt)

    def list_alt(self):
        """Display the list of the screens altitude"""

        cdef np.ndarray alts = np.zeros(self.s_a.nscreens, dtype=np.float32)
        self.s_a.list_alt(< float * > alts.data)
        return alts

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
