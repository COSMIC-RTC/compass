import numpy as np
cimport numpy as np
# np.import_array()

cdef class Acquisition:

    def __cinit__(self, int wfs_num, Sensors sensors):

        self.acquisition = new sutra_acquisim(sensors.sensors, wfs_num)

    def __dealloc__(self):
        del self.acquisition



    def comp_image(self, np.ndarray[ndim=2, dtype=np.float32_t] bimage):
        """Return the bincube

        :return:
           bincube : (np.ndarray[ndim=2,dtype=np.float32_t]) : cube of sub-apertures

        """

        cdef long dims[3] 
        dims[0] = 2
        dims[1] = bimage.shape[0]
        dims[2] = bimage.shape[1]

        self.acquisition.comp_image(dims, < float *> bimage.data)
