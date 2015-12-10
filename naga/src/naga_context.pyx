
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

    def __cinit__(self):
        self.c = carma_context.instance()

    def get_ndevice(self):
        """Return number of device."""
        return self.c.get_ndevice()

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


