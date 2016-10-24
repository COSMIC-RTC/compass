import numpy as np
cimport numpy as np
np.import_array()

from libcpp.string cimport string

"""
#################################################
# C-Class carma_context
#################################################
cdef extern from "carma_context.h":
    cdef cppclass carma_context "carma_context":
        carma_context()
        int get_ndevice()
        int set_activeDevice(int newDevice, int silent)
        int set_activeDeviceForce(int newDevice, int silent)
        int set_activeDeviceForCpy(int newDevice, int silent)
        int get_activeDevice()
"""
#################################################
# P-Class naga_context
#################################################
cdef class naga_context:

    def __cinit__(self, device=None, np.ndarray[ndim=1, dtype=np.int32_t] devices=None):
        if device is not None :
            self.c = &carma_context.instance_1gpu(device)
        elif devices is not None :
            self.c = &carma_context.instance_ngpu(devices.size, <np.int32_t*>devices.data)
        else :
            self.c = &carma_context.instance()

    def get_ndevice(self):
        """Return number of device."""
        return self.c.get_ndevice()

    def get_device_names(self):
        """Return names of devices."""
        cdef int nb_dev = self.c.get_ndevice()
        names = ['']*nb_dev
        for num_dev in range(nb_dev):
            names[num_dev]=self.c.get_device(num_dev).getName()
        return names

    def set_activeDevice(self, int newDevice, int silent=1):
        """Activate a device.

        newDevice -- int device to activate
        silent    -- int (default=1)
        """
        return self.c.set_activeDevice(newDevice, silent)

    def set_activeDeviceForce(self, newDevice, int silent=1):
        """Activate a device.

        newDevice -- int device to activate
        silent    -- int (default=1)
        """
        return self.c.set_activeDeviceForce(newDevice, silent)

    def set_activeDeviceForCpy(self, int newDevice, int silent=1):
        """Activate a device.

        newDevice -- int device to activate
        silent    -- int (default=1)
        """
        return self.c.set_activeDeviceForCpy(newDevice, silent)

    def get_activeDevice(self):
        """Return the index of actual activated device"""
        return self.c.get_activeDevice()

    def get_cudaRuntimeGetVersion(self):
        """Return the version number of the installed CUDA Runtime"""
        return self.c.get_cudaRuntimeGetVersion() / 1000.

    def get_cudaDriverGetVersion(self):
        """Return the version number of the installed CUDA driver"""
        return self.c.get_cudaDriverGetVersion() / 1000.

    def get_magma_info(self):
        """Return the information of the installed MAGMA"""
        return self.c.magma_info()
