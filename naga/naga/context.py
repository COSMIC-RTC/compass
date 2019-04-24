''' @package naga.context

More details.
'''

from carmaWrap import context
import numpy as np


class Context:
    """
    Python class for wrapping a carma_context
    """

    def __init__(self, devices=None):
        """ Initializes a Context object

        Parameters
        ------------
        devices: (list): (optionnal) list of devices to use. Default is 0
        """
        if devices is None:
            self.context = context.get_instance_1gpu(0)
            self.devices = 0
        else:
            if isinstance(devices, list):
                devices = np.array(devices, dtype=np.int32)
            else:
                raise TypeError("Devices must be a list of integer")
            self.devices = devices
            self.context = context.get_instance_ngpu(devices)

    def getActiveDevice(self):
        """ Return the index of the current active device """
        return self.context.activeDevice

    def setActiveDevice(self, index: int):
        """ Set the device index as the active device """
        if index in self.devices:
            self.context.set_activeDevice(index)
        else:
            raise ValueError("Index given is not valid")

    def enableTensorCores(self):
        """ Enable the tensor cores math mode """
        self.context.activate_tensor_cores(True)

    def disableTensorCores(self):
        """ Disable the tensor cores math mode """
        self.context.activate_tensor_cores(True)
