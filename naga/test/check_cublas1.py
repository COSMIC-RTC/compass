import naga as ch
import numpy as np
import numpy.testing as npt
import time

size = 1024
dec = 4
prec = 10**-dec

print("Test cublas 1")
print("precision: ", prec)

seed = 1234  # int(time.clock()*10**6)
c = ch.naga_context.get_instance()

# aa = np.random.random((size * size)).astype(np.float32)
# bb = ch.naga_obj_float(c, aa)

# print(aa)
# print(bb)
# print(np.array(bb))


def test_imax():
    Vect = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect.init_prng(seed)
    Vect.prng('U')

    #imax return the index in column major of the maximum absolute value element
    imaxC = np.array(Vect).argmax() + 1
    imaxG = Vect.imax(1)

    npt.assert_equal(imaxC, imaxG)


def test_imin():
    Vect = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect.init_prng(seed * 2)
    Vect.prng('U')

    #imax return the index in column major of the minimum obsolute value element
    iminC = np.array(Vect).argmin() + 1
    iminG = Vect.imin(1)

    npt.assert_equal(iminC, iminG)


def test_asum():
    Vect = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect.init_prng(seed * 3)
    Vect.prng('U')

    sumC = np.abs(np.array(Vect)).sum()
    sumG = Vect.asum(1)

    M = np.max(np.abs(sumC))
    d = 1
    if (M > 0):
        d = 10**np.ceil(np.log10(M))
    npt.assert_almost_equal(sumC / d, sumG / d, decimal=dec)


def test_nrm2():
    Vect = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect.init_prng(seed * 4)
    Vect.prng('U')

    nC = np.linalg.norm(np.array(Vect))
    nG = Vect.nrm2(1)

    M = np.max(np.abs(nC))
    d = 1
    if (M > 0):
        d = 10**np.ceil(np.log10(M))
    npt.assert_almost_equal(nC / d, nG / d, decimal=dec)


def test_scale():
    Vect = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect.init_prng(seed * 5)
    Vect.prng('U')

    scale = 1.67
    sC = np.array(Vect) * scale
    Vect.scale(scale, 1)
    sG = np.array(Vect)
    npt.assert_array_equal(sC, sG)


def test_swap():
    Vect = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect.init_prng(seed * 6)
    Vect.prng('U')
    h_Vect = np.array(Vect)

    Vect2 = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect2.init_prng(seed * 7)
    Vect2.prng('U')
    h_Vect2 = np.array(Vect2)

    Vect.swap(Vect2, 1, 1)
    npt.assert_array_equal(h_Vect2, np.array(Vect))
    npt.assert_array_equal(h_Vect, np.array(Vect2))


def test_copy():
    Vect = ch.naga_obj_float(c, np.zeros((size * size)).astype(np.float32))
    Vect.init_prng(seed * 8)
    Vect.prng('U')

    Vect2 = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect2.copy(Vect, 1, 1)
    npt.assert_array_equal(np.array(Vect), np.array(Vect2))


def test_axpy():
    Vect = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect.init_prng(seed * 9)
    Vect.prng('U')

    h_Vect = np.array(Vect)
    alpha = 1.4

    Vect2 = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect2.init_prng(seed * 10)
    Vect2.prng('U')
    h_Vect2 = np.array(Vect2)

    Vect.axpy(alpha, Vect2, 1, 1)

    h_Vect = h_Vect + alpha * h_Vect2

    npt.assert_almost_equal(h_Vect, np.array(Vect), decimal=dec)


def test_dot():
    Vect = ch.naga_obj_float(c, np.zeros((size * size), dtype=np.float32))
    Vect.init_prng(seed * 11)
    Vect.prng('U')

    h_Vect = np.array(Vect)
    dotC = np.dot(h_Vect, h_Vect)
    dotG = Vect.dot(Vect, 1, 1)

    M = np.max(np.abs(dotC))
    d = 1
    if (M > 0):
        d = 10**np.ceil(np.log10(M))
    npt.assert_almost_equal(dotC / d, dotG / d, decimal=dec)
