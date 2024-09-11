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
import carma as ch
import numpy.testing as npt
import time

dec = 4
prec = 10**-dec

sizem = 512
sizen = 1024

seed = np.int32(time.perf_counter())

print("")
print("Test cublas 2")
print("Precision: ", prec)

c = ch.context.get_instance()

# generatig random context.obj 2d

# generating random symetric context.obj 2d

# generating 3 random context.obj 1d


def test_float_gemv():
    # function gemv
    # testing: y=A.x
    # x and y are vector, A a matrix

    alpha = 2
    beta = 1

    Mat = ch.obj_float(c, np.random.randn(sizem, sizen))
    MatT = ch.obj_float(c, np.random.randn(sizen, sizem))
    # Mat.random_host(seed, 'U')
    # MatT.random_host(seed + 2, 'U')

    Vectx = ch.obj_float(c, np.random.randn(sizen))
    Vecty = ch.obj_float(c, np.random.randn(sizem))
    # Vectx.random_host(seed + 3, 'U')
    # Vecty.random_host(seed + 4, 'U')

    A = np.array(Mat)
    AT = np.array(MatT)
    x = np.array(Vectx)
    y = np.array(Vecty)

    y = alpha * A.dot(x) + beta * y
    y2 = alpha * A.dot(x)
    y3 = alpha * AT.T.dot(x)

    # Vecty = ch.obj_float(c, np.random.randn((sizem)))

    Mat.gemv(Vectx, alpha, "N", Vecty, beta)
    Vecty_2 = Mat.gemv(Vectx, alpha, "N")
    Vecty_3 = MatT.gemv(Vectx, alpha, "T")

    npt.assert_array_almost_equal(y, np.array(Vecty), decimal=dec - 1)
    npt.assert_array_almost_equal(y2, np.array(Vecty_2), decimal=dec - 1)
    npt.assert_array_almost_equal(y3, np.array(Vecty_3), decimal=dec - 1)


def test_float_ger():
    # function ger
    # testing: A= x.y
    #   and  : A= x.y+ A
    # x and y are vectors, A a matrix
    Mat = ch.obj_float(c, np.random.randn(sizem, sizen))
    # Mat.random_host(seed + 2, 'U')

    Vectx = ch.obj_float(c, np.random.randn(sizen))
    # Vectx.random_host(seed + 3, 'U')

    Vecty = ch.obj_float(c, np.random.randn(sizem))
    # Vecty.random_host(seed + 4, 'U')

    x = np.array(Vectx)
    A = np.array(Mat)
    y = np.array(Vecty)

    caOresA = Vecty.ger(Vectx, Mat)
    caOresB = Vecty.ger(Vectx)

    A = np.outer(y, x) + A
    B = np.outer(y, x)

    # npt.assert_array_almost_equal(A, np.array(caOresA), decimal=dec)
    npt.assert_array_almost_equal(B, np.array(caOresB), decimal=dec)


def test_float_symv():
    # function symv
    # testing: y=A.x
    # x and y are vector, A a symetric matrix

    MatSym = ch.obj_float(c, np.random.randn(sizem, sizem))
    # MatSym.random_host(seed + 2, 'U')
    data_R = np.array(MatSym)
    data_R = data_R + data_R.T
    MatSym = ch.obj_float(c, data_R)

    Vectx = ch.obj_float(c, np.random.randn((sizem)))
    # Vectx.random_host(seed + 5, 'U')

    Vecty = ch.obj_float(c, np.random.randn((sizem)))

    A = np.array(MatSym)

    x2 = np.array(Vectx)

    y = A.dot(x2)

    MatSym.symv(Vectx, vecty=Vecty)
    Vecty_2 = MatSym.symv(Vectx)

    npt.assert_array_almost_equal(y, np.array(Vecty), decimal=dec)
    npt.assert_array_almost_equal(y, np.array(Vecty_2), decimal=dec)


def test_double_gemv():
    # function gemv
    # testing: y=A.x
    # x and y are vector, A a matrix

    alpha = 2
    beta = 1

    Mat = ch.obj_double(c, np.random.randn(sizem, sizen))
    MatT = ch.obj_double(c, np.random.randn(sizen, sizem))
    # Mat.random_host(seed, 'U')
    # MatT.random_host(seed + 2, 'U')

    Vectx = ch.obj_double(c, np.random.randn(sizen))
    Vecty = ch.obj_double(c, np.random.randn(sizem))
    # Vectx.random_host(seed + 3, 'U')
    # Vecty.random_host(seed + 4, 'U')

    A = np.array(Mat)
    AT = np.array(MatT)
    x = np.array(Vectx)
    y = np.array(Vecty)

    y = alpha * A.dot(x) + beta * y
    y2 = alpha * A.dot(x)
    y3 = alpha * AT.T.dot(x)

    # Vecty = ch.obj_double(c, np.random.randn((sizem)))

    Mat.gemv(Vectx, alpha, "N", Vecty, beta)
    Vecty_2 = Mat.gemv(Vectx, alpha, "N")
    Vecty_3 = MatT.gemv(Vectx, alpha, "T")

    npt.assert_array_almost_equal(y, np.array(Vecty), decimal=dec - 1)
    npt.assert_array_almost_equal(y2, np.array(Vecty_2), decimal=dec - 1)
    npt.assert_array_almost_equal(y3, np.array(Vecty_3), decimal=dec - 1)


def test_double_ger():
    # function ger
    # testing: A= x.y
    #   and  : A= x.y+ A
    # x and y are vectors, A a matrix
    Mat = ch.obj_double(c, np.random.randn(sizem, sizen))
    # Mat.random_host(seed + 2, 'U')

    Vectx = ch.obj_double(c, np.random.randn(sizen))
    # Vectx.random_host(seed + 3, 'U')

    Vecty = ch.obj_double(c, np.random.randn(sizem))
    # Vecty.random_host(seed + 4, 'U')

    x = np.array(Vectx)
    A = np.array(Mat)
    y = np.array(Vecty)

    caOresA = Vecty.ger(Vectx, Mat)
    caOresB = Vecty.ger(Vectx)

    A = np.outer(y, x) + A
    B = np.outer(y, x)

    # npt.assert_array_almost_equal(A, np.array(caOresA), decimal=dec)
    npt.assert_array_almost_equal(B, np.array(caOresB), decimal=dec)


def test_double_symv():
    # function symv
    # testing: y=A.x
    # x and y are vector, A a symetric matrix

    MatSym = ch.obj_double(c, np.random.randn(sizem, sizem))
    # MatSym.random_host(seed + 2, 'U')
    data_R = np.array(MatSym)
    data_R = data_R + data_R.T
    MatSym = ch.obj_double(c, data_R)

    Vectx = ch.obj_double(c, np.random.randn((sizem)))
    # Vectx.random_host(seed + 5, 'U')

    Vecty = ch.obj_double(c, np.random.randn((sizem)))

    A = np.array(MatSym)

    x2 = np.array(Vectx)

    y = A.dot(x2)

    MatSym.symv(Vectx, vecty=Vecty)
    Vecty_2 = MatSym.symv(Vectx)

    npt.assert_array_almost_equal(y, np.array(Vecty), decimal=dec)
    npt.assert_array_almost_equal(y, np.array(Vecty_2), decimal=dec)
