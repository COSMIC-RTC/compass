import numpy as np
import naga as ch
import numpy.testing as npt
import time

dec = 4
prec = 10**-dec

sizem = 512
sizen = 1024

seed = 1234  # int(time.clock()*10**6)

print("")
print("Test cublas 2")
print("Precision: ", prec)

c = ch.naga_context.get_instance()

#generatig random naga_obj 2d

#generating random symetric naga_obj 2d

#generating 3 random naga_obj 1d


def test_gemv():

    #function gemv
    # testing: y=A.x
    # x and y are vector, A a matrix

    A = np.random.randn(sizem, sizen)
    x = np.random.randn(sizen)

    Mat = ch.naga_obj_float(c, A)
    #Mat.random(seed)

    Vectx = ch.naga_obj_float(c, x)
    #Vectx.random(seed * 3)

    #A = np.array(Mat)
    #x = np.array(Vectx)
    y = A.dot(x)

    # Vecty = ch.naga_obj_float(c, np.zeros((sizem), dtype=np.float32))

    # Mat.gemv(Vectx, Vecty)  # , "N", 1, Vecty, sizem, 0)
    Vecty_2 = Mat.gemv(Vectx)

    # npt.assert_array_almost_equal(y, np.array(Vecty), decimal=dec)
    npt.assert_array_almost_equal(y, np.array(Vecty_2), decimal=dec)


def test_ger():
    # function ger
    # testing: A= x.y
    #   and  : A= x.y+ A
    # x and y are vectors, A a matrix
    Mat = ch.naga_obj_float(c, np.zeros((sizem, sizen), dtype=np.float32))
    # Mat.random(seed)

    Vectx = ch.naga_obj_float(c, np.zeros((sizen), dtype=np.float32))
    Vectx.random(seed * 3)

    Vecty = ch.naga_obj_float(c, np.zeros((sizem), dtype=np.float32))
    Vecty.random(seed * 4)

    x = np.array(Vectx)
    A = np.array(Mat)
    y = np.array(Vecty)

    caOresA = Vecty.ger(Vectx, Mat)
    caOresB = Vecty.ger(Vectx)

    A = np.outer(y, x) + A
    B = np.outer(y, x)

    # npt.assert_array_almost_equal(A, np.array(caOresA), decimal=dec)
    npt.assert_array_almost_equal(B, np.array(caOresB), decimal=dec)


def test_symv():
    #function symv
    # testing: y=A.x
    # x and y are vector, A a symetric matrix

    MatSym = ch.naga_obj_float(c, np.zeros((sizem, sizem), dtype=np.float32))
    MatSym.random(seed * 2)
    data_R = np.array(MatSym)
    data_R = data_R + data_R.T
    MatSym = ch.naga_obj_float(c, data_R)

    Vectx = ch.naga_obj_float(c, np.zeros((sizem), dtype=np.float32))
    Vectx.random(seed * 5)

    Vecty = ch.naga_obj_float(c, np.zeros((sizem), dtype=np.float32))

    A = np.array(MatSym)

    x2 = np.array(Vectx)

    y = A.dot(x2)

    MatSym.symv(Vectx, Vecty)
    Vecty_2 = MatSym.symv(Vectx)

    npt.assert_array_almost_equal(y, np.array(Vecty), decimal=dec)
    npt.assert_array_almost_equal(y, np.array(Vecty_2), decimal=dec)
