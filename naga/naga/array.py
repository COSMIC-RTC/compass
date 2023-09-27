## @package   naga.array
## @brief     Documentation for naga
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   5.5.0
## @date      2022/01/24
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
#  All rights reserved.
#  Distributed under GNU - LGPL
#
#  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
#  General Public License as published by the Free Software Foundation, either version 3 of the License,
#  or any later version.
#
#  COMPASS: End-to-end AO simulation tool using GPU acceleration
#  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
#
#  The final product includes a software package for simulating all the critical subcomponents of AO,
#  particularly in the context of the ELT and a real-time core based on several control approaches,
#  with performances consistent with its integration into an instrument. Taking advantage of the specific
#  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
#  conduct large simulation campaigns called to the ELT.
#
#  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
#  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
#  various systems configurations such as multi-conjugate AO.
#
#  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
#  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
#  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.

import numpy as np
from naga.context import Context

from carmaWrap import obj_float, obj_double, obj_int, obj_float_complex, obj_double_complex, obj_uint16
try:
    from carmaWrap import obj_half
    USE_HALF = 1
except:
    USE_HALF = 0
context = Context()


class DimensionsError(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Array():

    def __init__(self, data=None, shape=None, dtype=None):
        if shape is not None:
            if dtype is None:
                dtype = np.float32
            if isinstance(shape, tuple):
                data = np.zeros(shape, dtype=dtype)
            else:
                raise TypeError("shape must be a tuple")

        if data is not None:
            if isinstance(data, list):
                data = np.array(data)
            if isinstance(data, np.ndarray):
                if dtype is not None:
                    data = data.astype(dtype)
                if data.dtype == np.int64 or data.dtype == np.int32:
                    self.__data = obj_int(context.context, data)
                elif data.dtype == np.float32:
                    self.__data = obj_float(context.context, data)
                elif data.dtype == np.float64:
                    self.__data = obj_double(context.context, data)
                elif USE_HALF and data.dtype == np.float16:
                    self.__data = obj_half(context.context, data)
                elif data.dtype == np.complex64:
                    self.__data = obj_float_complex(context.context, data)
                elif data.dtype == np.complex128:
                    self.__data = obj_double_complex(context.context, data)
                else:
                    raise TypeError("Data type not implemented")
                self.__dtype = data.dtype
                self.__shape = data.shape
            elif isinstance(data, obj_int):
                self.__data = data
                self.__dtype = np.int64
                self.__shape = tuple(data.shape[k] for k in range(len(data.shape)))
            elif isinstance(data, obj_float):
                self.__data = data
                self.__dtype = np.float32
                self.__shape = tuple(data.shape[k] for k in range(len(data.shape)))
            elif isinstance(data, obj_double):
                self.__data = data
                self.__dtype = np.float64
                self._shape = tuple(data.shape[k] for k in range(len(data.shape)))
            elif USE_HALF and isinstance(data, obj_half):
                self.__data = data
                self.__dtype = np.float16
                self.__shape = tuple(data.shape[k] for k in range(len(data.shape)))
            elif isinstance(data, obj_float_complex):
                self.__data = data
                self.__dtype = np.complex64
                self.__shape = tuple(data.shape[k] for k in range(len(data.shape)))
            elif isinstance(data, obj_double_complex):
                self.__data = data
                self.__dtype = np.complex128
                self.__shape = tuple(data.shape[k] for k in range(len(data.shape)))
            elif isinstance(data, obj_uint16):
                self.__data = data
                self.__dtype = np.uint16
                self.__shape = tuple(data.shape[k] for k in range(len(data.shape)))
            else:
                raise TypeError("Data must be a list, a numpy array or a carmaWrap.obj")
            self.__size = self.__data.nbElem
        else:
            raise AttributeError("You must provide data or shape at least")

    shape = property(lambda x: x.__shape)
    dtype = property(lambda x: x.__dtype)
    data = property(lambda x: x.__data)

    def __repr__(self):
        return "GPU" + self.toarray().__repr__()

    def __add__(self, idata):
        tmp = self.copy()
        if isinstance(idata, Array):
            tmp.data.axpy(1, idata.data)
        elif isinstance(idata, np.ndarray):
            tmp.data.axpy(1, Array(idata).data)
        else:
            raise TypeError("operator + is defined only between Arrays and np.arrays")
        return tmp

    def __sub__(self, idata):
        tmp = self.copy()
        if isinstance(idata, Array):
            tmp.data.axpy(-1, idata.data)
        elif isinstance(idata, np.ndarray):
            tmp.data.axpy(-1, Array(idata).data)
        else:
            raise TypeError("operator + is defined only between Arrays and np.arrays")
        return tmp

    def __mul__(self, idata):
        if isinstance(idata, float) or isinstance(idata, int):
            tmp = self.copy()
            tmp.data.scale(idata)
            return tmp
        else:
            raise NotImplementedError("Operator not implemented yet")

    def __getitem__(self, indices):
        return self.toarray().__getitem__(indices)

    def copy(self):
        tmp = Array(shape=self.shape, dtype=self.dtype)
        tmp.data.copy_from(self.data)
        return tmp

    def dot(self, idata):
        if isinstance(idata, np.ndarray):
            if idata.dtype == self.dtype:
                idata = Array(idata)
            else:
                raise TypeError("Data types must be the same for both arrays")
        if isinstance(idata, Array):
            if len(self.shape) == 1:
                if len(idata.shape) == 1:
                    if idata.shape == self.shape:
                        result = self.data.dot(idata.data, 1, 1)
                    else:
                        raise DimensionsError("Dimensions mismatch")
                elif len(idata.shape) == 2:
                    if idata.shape[0] == self.shape[0]:
                        result = Array(idata.data.gemv(self.data, op='T'))
                    else:
                        raise DimensionsError("Dimensions mismatch")
                else:
                    raise DimensionsError("Arrays must be 1D or 2D max")
            elif len(self.shape) == 2:
                if len(idata.shape) == 1:
                    if idata.shape[0] == self.shape[1]:
                        result = Array(self.data.gemv(idata.data))
                    else:
                        raise DimensionsError("Dimensions mismatch")
                elif len(idata.shape) == 2:
                    if self.shape[1] == idata.shape[0]:
                        result = Array(self.data.gemm(idata.data))
                    else:
                        raise DimensionsError("Dimensions mismatch")
            else:
                raise DimensionsError("Arrays must be 1D or 2D max")

        return result

    def argmax(self):
        return self.data.aimax()

    def max(self):
        return self.toarray()[self.argmax()]

    def argmin(self):
        return self.data.aimin()

    def min(self):
        return self.toarray()[self.argmin()]

    def sum(self):
        return self.data.sum()

    def toarray(self):
        tmp = np.array(self.data)
        return tmp


def ones(shape, dtype=np.float32):
    return Array(np.ones(shape, dtype=dtype))


def zeros(shape, dtype=np.float32):
    return Array(np.zeros(shape, dtype=dtype))
