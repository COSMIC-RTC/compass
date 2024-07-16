## @package   carma.test
## @brief     Unit tests for carma
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @date      2022/01/24
## @copyright 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>
#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>

# COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation, either version 3 of the 
# License, or any later version.

# COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
# If not, see <https://www.gnu.org/licenses/>

# Copyright (C) 2011-2024 COSMIC Team <https//://github.com/COSMIC-RTC/compass>

import carma as ch
import numpy as np
import numpy.testing as npt
import time

size = 1024
dec = 4
prec = 10**-dec

print("Test cublas 1")
print("precision: ", prec)

seed = np.int32(time.perf_counter() * 1e3)
c = ch.context.get_instance()

# aa = np.random.random((size * size)).astype(np.float32)
# bb = ch.obj_float(c, aa)

# print(aa)
# print(bb)
# print(np.array(bb))


def test_float_aimax():
    Vect = ch.obj_float(c, np.random.randn(size))
    # Vect = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed, 'U')
    v = np.array(Vect)

    #imax return the index in column major of the maximum absolute value element
    imaxC = np.abs(v).argmax()
    imaxG = Vect.aimax()

    print(v[imaxC], v[imaxG])
    npt.assert_equal(imaxC, imaxG)


def test_float_aimin():
    Vect = ch.obj_float(c, np.random.randn(size))
    # Vect = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*2, 'U')
    v = np.array(Vect)

    #imax return the index in column major of the minimum obsolute value element
    iminC = np.abs(v).argmin()
    iminG = Vect.aimin()

    print(v[iminC], v[iminG])
    npt.assert_equal(iminC, iminG)


def test_float_asum():
    Vect = ch.obj_float(c, np.random.randn(size))
    # Vect = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*3, 'U')
    v = np.array(Vect)

    sumC = np.sum(np.abs(v))
    sumG = Vect.asum()
    print(sumC, sumG)
    npt.assert_almost_equal(sumC, sumG, decimal=dec)

    # M = np.max(np.abs(sumC))
    # d = 1
    # if (M > 0):
    #     d = 10**np.ceil(np.log10(M))
    # npt.assert_almost_equal(sumC / d, sumG / d, decimal=dec)


def test_float_nrm2():
    Vect = ch.obj_float(c, np.random.randn(size))
    # Vect = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*4, 'U')
    v = np.array(Vect)

    nC = np.linalg.norm(v)
    nG = Vect.nrm2(1)

    npt.assert_almost_equal(nC, nG, decimal=dec)

    # M = np.max(np.abs(nC))
    # d = 1
    # if (M > 0):
    #     d = 10**np.ceil(np.log10(M))
    # npt.assert_almost_equal(nC / d, nG / d, decimal=dec)


def test_float_scale():
    Vect = ch.obj_float(c, np.random.randn(size))
    # Vect = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*5, 'U')
    v = np.array(Vect)

    scale = 1.67
    sC = v * scale
    Vect.scale(scale, 1)
    sG = np.array(Vect)
    npt.assert_array_almost_equal(sC, sG, decimal=dec)


def test_float_swap():
    Vect = ch.obj_float(c, np.random.randn(size))
    # Vect = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*6, 'U')
    h_Vect = np.array(Vect)

    Vect2 = ch.obj_float(c, np.random.randn(size))
    # Vect2 = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect2.random_host(seed*7, 'U')
    h_Vect2 = np.array(Vect2)

    Vect.swap(Vect2, 1, 1)
    npt.assert_equal(h_Vect2, np.array(Vect))
    npt.assert_equal(h_Vect, np.array(Vect2))


def test_float_copy():
    Vect = ch.obj_float(c, np.random.randn(size))
    # Vect = ch.obj_float(c, np.empty((size)))
    # Vect.random_host(seed*8, 'U')

    Vect2 = ch.obj_float(c, np.random.randn(size))
    # Vect2 = ch.obj_float(c, np.empty((size), dtype=np.float32))
    Vect2.copy(Vect, 1, 1)
    npt.assert_array_equal(np.array(Vect), np.array(Vect2))


def test_float_axpy():
    Vect = ch.obj_float(c, np.random.randn(size))
    # Vect = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*9, 'U')
    h_Vect = np.array(Vect)

    alpha = 1.4

    Vect2 = ch.obj_float(c, np.random.randn(size))
    # Vect2 = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect2.random_host(seed*10, 'U')
    h_Vect2 = np.array(Vect2)

    Vect.axpy(alpha, Vect2, 1, 1)

    h_Vect = h_Vect + alpha * h_Vect2

    npt.assert_almost_equal(h_Vect, np.array(Vect), decimal=dec)


def test_float_dot():
    Vect = ch.obj_float(c, np.random.randn(size))
    # Vect = ch.obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*11, 'U')
    h_Vect = np.array(Vect)

    dotC = np.dot(h_Vect, h_Vect)
    dotG = Vect.dot(Vect, 1, 1)

    M = np.max(np.abs(dotC))
    d = 1
    if (M > 0):
        d = 10**np.ceil(np.log10(M))
    npt.assert_almost_equal(dotC / d, dotG / d, decimal=dec)


def test_double_aimax():
    Vect = ch.obj_double(c, np.random.randn(size))
    # Vect = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed, 'U')
    v = np.array(Vect)

    #imax return the index in column major of the maximum absolute value element
    imaxC = np.abs(v).argmax()
    imaxG = Vect.aimax()

    print(v[imaxC], v[imaxG])
    npt.assert_equal(imaxC, imaxG)


def test_double_aimin():
    Vect = ch.obj_double(c, np.random.randn(size))
    # Vect = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*2, 'U')
    v = np.array(Vect)

    #imax return the index in column major of the minimum obsolute value element
    iminC = np.abs(v).argmin()
    iminG = Vect.aimin()

    print(v[iminC], v[iminG])
    npt.assert_equal(iminC, iminG)


def test_double_asum():
    Vect = ch.obj_double(c, np.random.randn(size))
    # Vect = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*3, 'U')
    v = np.array(Vect)

    sumC = np.sum(np.abs(v))
    sumG = Vect.asum()
    print(sumC, sumG)
    npt.assert_almost_equal(sumC, sumG, decimal=2 * dec)

    # M = np.max(np.abs(sumC))
    # d = 1
    # if (M > 0):
    #     d = 10**np.ceil(np.log10(M))
    # npt.assert_almost_equal(sumC / d, sumG / d, decimal=2*dec)


def test_double_nrm2():
    Vect = ch.obj_double(c, np.random.randn(size))
    # Vect = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*4, 'U')
    v = np.array(Vect)

    nC = np.linalg.norm(v)
    nG = Vect.nrm2(1)

    npt.assert_almost_equal(nC, nG, decimal=2 * dec)

    # M = np.max(np.abs(nC))
    # d = 1
    # if (M > 0):
    #     d = 10**np.ceil(np.log10(M))
    # npt.assert_almost_equal(nC / d, nG / d, decimal=2*dec)


def test_double_scale():
    Vect = ch.obj_double(c, np.random.randn(size))
    # Vect = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*5, 'U')
    v = np.array(Vect)

    scale = 1.67
    sC = v * scale
    Vect.scale(scale, 1)
    sG = np.array(Vect)
    npt.assert_array_almost_equal(sC, sG, decimal=2 * dec)


def test_double_swap():
    Vect = ch.obj_double(c, np.random.randn(size))
    # Vect = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*6, 'U')
    h_Vect = np.array(Vect)

    Vect2 = ch.obj_double(c, np.random.randn(size))
    # Vect2 = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect2.random_host(seed*7, 'U')
    h_Vect2 = np.array(Vect2)

    Vect.swap(Vect2, 1, 1)
    npt.assert_equal(h_Vect2, np.array(Vect))
    npt.assert_equal(h_Vect, np.array(Vect2))


def test_double_copy():
    Vect = ch.obj_double(c, np.random.randn(size))
    # Vect = ch.obj_double(c, np.empty((size)))
    # Vect.random_host(seed*8, 'U')

    Vect2 = ch.obj_double(c, np.random.randn(size))
    # Vect2 = ch.obj_double(c, np.empty((size), dtype=np.float32))
    Vect2.copy(Vect, 1, 1)
    npt.assert_array_equal(np.array(Vect), np.array(Vect2))


def test_double_axpy():
    Vect = ch.obj_double(c, np.random.randn(size))
    # Vect = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*9, 'U')
    h_Vect = np.array(Vect)

    alpha = 1.4

    Vect2 = ch.obj_double(c, np.random.randn(size))
    # Vect2 = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect2.random_host(seed*10, 'U')
    h_Vect2 = np.array(Vect2)

    Vect.axpy(alpha, Vect2, 1, 1)

    h_Vect = h_Vect + alpha * h_Vect2

    npt.assert_almost_equal(h_Vect, np.array(Vect), decimal=2 * dec)


def test_double_dot():
    Vect = ch.obj_double(c, np.random.randn(size))
    # Vect = ch.obj_double(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*11, 'U')
    h_Vect = np.array(Vect)

    dotC = np.dot(h_Vect, h_Vect)
    dotG = Vect.dot(Vect, 1, 1)

    M = np.max(np.abs(dotC))
    d = 1
    if (M > 0):
        d = 10**np.ceil(np.log10(M))
    npt.assert_almost_equal(dotC / d, dotG / d, decimal=2 * dec)
