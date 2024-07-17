#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
#
# COMPASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COMPASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COMPASS. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2024 COSMIC Team


import numpy as np
from typing import Union, Type, Tuple


def enforce_int(number: Union[int, np.int32, np.int64]) -> np.int32:
    """
    Enforce that the input number is an integer.

    Args:
        number: The number to enforce as an integer.

    Returns:
        The input number as an integer.

    Raises:
        TypeError: If the input number is not an integer.
    """
    if not isinstance(number, (int, np.int32, np.int64)):
        raise TypeError("Value should be an integer.")
    return np.int32(number)


def enforce_float(number: Union[float, np.float32, np.float64]) -> np.float32:
    """
    Enforce that the input number is a float.

    Args:
        number: The number to enforce as a float.

    Returns:
        The input number as a float.

    Raises:
        TypeError: If the input number is not a float.
    """
    if isinstance(number, (int, np.int32, np.int64)):
        number = np.float32(number)
    if not isinstance(number, (float, np.float32, np.float64)):
        raise TypeError("Value should be a float.")
    return np.float32(number)


def enforce_or_cast_bool(
    number: Union[bool, int, np.int32, np.int64, float, np.float32, np.float64],
) -> bool:
    """
    Enforce that the input number is a boolean or cast it to a boolean if possible.

    Args:
        number: The number to enforce or cast as a boolean.

    Returns:
        The input number as a boolean.

    Raises:
        ValueError: If the input number is not 0 or 1.
        TypeError: If the input number cannot be cast to a boolean.
    """
    if isinstance(number, bool):
        return number
    if isinstance(number, (int, np.int32, np.int64, float, np.float32, np.float64)):
        if number == 0:
            return False
        elif number == 1:
            return True
        else:
            raise ValueError("Will not cast non 0/1 int or float to boolean.")
    raise TypeError("Will only cast int and float to booleans.")


def enforce_array(
    data: Union[
        int,
        np.int32,
        np.int64,
        float,
        np.float32,
        np.float64,
        complex,
        np.complex64,
        np.complex128,
        list,
        np.ndarray,
    ],
    size: int,
    dtype: Type[np.float32] = np.float32,
    scalar_expand: bool = False,
) -> np.ndarray:
    """
    Enforce that the input data is an array of a specific size.

    Args:
        data: The data to enforce as an array.
        size: The expected size of the array.
        dtype: The data type of the array (default: np.float32).
        scalar_expand: Whether to expand a scalar input to an array of the specified size (default: False).

    Returns:
        The input data as a numpy array.

    Raises:
        TypeError: If the input data is not a valid array.
    """
    if isinstance(
        data,
        (
            int,
            np.int32,
            np.int64,
            float,
            np.float32,
            np.float64,
            complex,
            np.complex64,
            np.complex128,
        ),
    ):
        data = [data]
    if len(data) == 1:
        if scalar_expand or size == 1:
            return np.full(size, data[0], dtype=dtype)
        else:
            raise TypeError("This non-singleton array cannot be initialized with a scalar.")
    if len(data) != size:
        raise TypeError("Input argument has wrong number of elements.")
    if isinstance(data, np.ndarray) and len(data.shape) > 1:
        raise TypeError("Multidimensional ndarray input is not allowed")
    if isinstance(data, list) and not all([isinstance(x, (float, int, complex)) for x in data]):
        raise TypeError("Input list may only contain numerical values.")
    return np.array(data, dtype=dtype)


def enforce_arrayMultiDim(
    data: np.ndarray, shape: Tuple[int], dtype: Type[np.float32] = np.float32
) -> np.ndarray:
    """
    Enforce that the input data is a multi-dimensional array with a specific shape.

    Args:
        data: The data to enforce as a multi-dimensional array.
        shape: The expected shape of the array.
        dtype: The data type of the array (default: np.float32).

    Returns:
        The input data as a numpy array.

    Raises:
        TypeError: If the input data is not a valid multi-dimensional array.
    """
    if not isinstance(data, np.ndarray):
        raise TypeError("Input argument must be a np.ndarray")
    else:
        doRaise = False
        if len(data.shape) != len(shape):
            doRaise = True
        for i, j in zip(data.shape, shape):
            if j != -1 and i != j:
                doRaise = True
        if doRaise:
            raise TypeError("Input has wrong dimensions, expect multi-dimensional arrays")
    return np.array(data, dtype=dtype)
