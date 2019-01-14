import carmaWrap as ch
import numpy as np
import numpy.testing as npt
import time

size = 1024
dec = 4
prec = 10**-dec

print("Test cublas 1")
print("precision: ", prec)

seed = 1234  # int(time.clock()*10**6)
c = ch.carmaWrap_context.get_instance()

# aa = np.random.random((size * size)).astype(np.float32)
# bb = ch.carmaWrap_obj_float(c, aa)

# print(aa)
# print(bb)
# print(np.array(bb))


def test_aimax():
    v = np.random.randn(size)
    Vect = ch.carmaWrap_obj_float(c, v)
    # Vect = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed, 'U')
    v = np.array(Vect)

    #imax return the index in column major of the maximum absolute value element
    imaxC = np.abs(v).argmax()
    imaxG = Vect.aimax()

    print(v[imaxC], v[imaxG])
    npt.assert_equal(imaxC, imaxG)


def test_aimin():
    v = np.random.randn(size)
    Vect = ch.carmaWrap_obj_float(c, v)
    # Vect = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*2, 'U')
    v = np.array(Vect)

    #imax return the index in column major of the minimum obsolute value element
    iminC = np.abs(v).argmin()
    iminG = Vect.aimin()

    print(v[iminC], v[iminG])
    npt.assert_equal(iminC, iminG)


def test_asum():
    v = np.random.randn(size)
    Vect = ch.carmaWrap_obj_float(c, v)
    # Vect = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
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


def test_nrm2():
    v = np.random.randn(size)
    Vect = ch.carmaWrap_obj_float(c, v)
    # Vect = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
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


def test_scale():
    v = np.random.randn(size)
    Vect = ch.carmaWrap_obj_float(c, v)
    # Vect = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*5, 'U')
    v = np.array(Vect)

    scale = 1.67
    sC = v * scale
    Vect.scale(scale, 1)
    sG = np.array(Vect)
    npt.assert_array_almost_equal(sC, sG, decimal=dec)


def test_swap():
    h_Vect = np.random.randn(size)
    Vect = ch.carmaWrap_obj_float(c, h_Vect)
    # Vect = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*6, 'U')
    h_Vect = np.array(Vect)

    h_Vect2 = np.random.randn(size)
    Vect2 = ch.carmaWrap_obj_float(c, h_Vect2)
    # Vect2 = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
    # Vect2.random_host(seed*7, 'U')
    h_Vect2 = np.array(Vect2)

    Vect.swap(Vect2, 1, 1)
    npt.assert_array_equal(h_Vect2, np.array(Vect))
    npt.assert_array_equal(h_Vect, np.array(Vect2))


def test_copy():
    Vect = ch.carmaWrap_obj_float(c, np.random.randn(size))
    # Vect = ch.carmaWrap_obj_float(c, np.empty((size)).astype(np.float32))
    # Vect.random_host(seed*8, 'U')

    Vect2 = ch.carmaWrap_obj_float(c, np.random.randn(size))
    # Vect2 = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
    Vect2.copy(Vect, 1, 1)
    npt.assert_array_equal(np.array(Vect), np.array(Vect2))


def test_axpy():
    Vect = ch.carmaWrap_obj_float(c, np.random.randn(size))
    # Vect = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*9, 'U')

    h_Vect = np.array(Vect)
    alpha = 1.4

    Vect2 = ch.carmaWrap_obj_float(c, np.random.randn(size))
    # Vect2 = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
    # Vect2.random_host(seed*10, 'U')
    h_Vect2 = np.array(Vect2)

    Vect.axpy(alpha, Vect2, 1, 1)

    h_Vect = h_Vect + alpha * h_Vect2

    npt.assert_almost_equal(h_Vect, np.array(Vect), decimal=dec)


def test_dot():
    Vect = ch.carmaWrap_obj_float(c, np.random.randn(size))
    # Vect = ch.carmaWrap_obj_float(c, np.empty((size), dtype=np.float32))
    # Vect.random_host(seed*11, 'U')

    h_Vect = np.array(Vect)
    dotC = np.dot(h_Vect, h_Vect)
    dotG = Vect.dot(Vect, 1, 1)

    M = np.max(np.abs(dotC))
    d = 1
    if (M > 0):
        d = 10**np.ceil(np.log10(M))
    npt.assert_almost_equal(dotC / d, dotG / d, decimal=dec)
