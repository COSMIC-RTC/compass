import numpy as np
import naga as ch
import numpy.testing as npt
import time

from naga.context import context
from naga.obj import obj_Float1D, obj_Float2D

dec = 4
prec = 10**-dec

sizem = 64
sizen = 128
sizek = 256

print("")
print("Test cublas 3")
print("precision: ", prec)

c = context()


def test_gemm():

    # function gemm
    #testing: C=A.B+C
    #A,B,C matrices

    #generating random matrices A,B,C and associated carma_obj
    matA = obj_Float2D(c, dims=np.array([sizem, sizek], dtype=np.int64))
    matB = obj_Float2D(c, dims=np.array([sizek, sizen], dtype=np.int64))
    matC = obj_Float2D(c, dims=np.array([sizem, sizen], dtype=np.int64))

    matA.random(time.clock() * 10**6)
    matB.random(time.clock() * 10**6)
    matC.random(time.clock() * 10**6)

    A = matA.device2host()
    B = matB.device2host()
    C = matC.device2host()

    #matrices product
    matA.gemm(matB, alpha=1, C=matC, beta=1)
    C = A.dot(B) + C

    #test results
    res = matC.device2host()

    M = np.argmax(np.abs(res - C))
    d = 5
    if (0 < np.abs(C.item(M))):
        d = 10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


def test_symm():
    #function symm
    #testing: C=A.B+C
    #A ssymetric matrix, B,C matrices

    #generating random matrices and associated carma_obj
    matA = obj_Float2D(c, dims=np.array([sizek, sizek], dtype=np.int64))
    matB = obj_Float2D(c, dims=np.array([sizek, sizen], dtype=np.int64))
    matC = obj_Float2D(c, dims=np.array([sizek, sizen], dtype=np.int64))

    matA.random(time.clock() * 10**6)
    matB.random(time.clock() * 10**6)
    matC.random(time.clock() * 10**6)

    #A symetric
    A = matA.device2host()
    A = (A + A.T) / 2
    matA.host2device(A)
    B = matB.device2host()
    C = matC.device2host()

    #matrices multiplication
    t1 = time.clock()
    matA.symm(matB, alpha=1, C=matC, beta=1)
    t2 = time.clock()
    C = A.dot(B) + C
    t3 = time.clock()

    print("")
    print("test symm:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #test results
    res = matC.device2host()

    M = np.argmax(np.abs(res - C))
    d = 5
    if (0 < np.abs(C.item(M))):
        d = 10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


def test_dgmm():
    #function dgmm
    #testing: C=A.d
    # C,A matrices, d vector (diagonal matrix as a vector)

    #generating random matrices and associated carma_obj
    matA = obj_Float2D(c, dims=np.array([sizek, sizek], dtype=np.int64))
    Vectx = obj_Float1D(c, dims=np.array([sizek], dtype=np.int64))

    matA.random(time.clock() * 10**6)
    A = matA.device2host()
    d = Vectx.device2host()

    #matrices product
    t1 = time.clock()
    matC = matA.dgmm(Vectx)
    t2 = time.clock()
    C = A * d
    t3 = time.clock()

    print("")
    print("test dgmm:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #test results
    res = matC.device2host()

    M = np.argmax(np.abs(res - C))
    d = 5
    if (0 < np.abs(C.item(M))):
        d = 10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


def test_syrk():
    #function syrk
    #testing: C=A.transpose(A)+C
    #A matrix, C symetric matrix

    #generating random matrices and associated carma_obj
    matA = obj_Float2D(c, dims=np.array([sizen, sizek], dtype=np.int64))
    matC = obj_Float2D(c, dims=np.array([sizen, sizen], dtype=np.int64))

    matA.random(time.clock() * 10**6)
    matC.random(time.clock() * 10**6)

    A = matA.device2host()
    C = matC.device2host()
    #syrk: C matrix is symetric
    C = (C + C.T) / 2

    matC.host2device(C)

    #matrices product
    t1 = time.clock()
    matA.syrk(C=matC, beta=1)
    t2 = time.clock()
    C = A.dot(A.T) + C
    t3 = time.clock()

    print("")
    print("test syrk:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #test results
    res = matC.device2host()

    #only upper triangle is computed
    C = np.triu(C).flatten()
    res = np.triu(res).flatten()

    M = np.argmax(np.abs(res - C))
    d = 5
    if (0 < np.abs(C.item(M))):
        d = 10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


def test_syrkx():
    #function syrkx
    #testing: C=A.transpose(B)+C
    #A matrix, C symetric matrix

    #generating random matrices and associated carma_obj
    matA = obj_Float2D(c, dims=np.array([sizen, sizek], dtype=np.int64))
    matB = obj_Float2D(c, dims=np.array([sizen, sizek], dtype=np.int64))
    matC = obj_Float2D(c, dims=np.array([sizen, sizen], dtype=np.int64))

    matA.random(time.clock() * 10**6)
    matB.copyFrom(matA)
    matC.random(time.clock() * 10**6)

    A = matA.device2host()
    B = matB.device2host()
    C = matC.device2host()

    #C is symetric
    C = np.dot(C, C.T)

    matC.host2device(C)

    #matrices product
    t1 = time.clock()
    matA.syrkx(matB, alpha=1, C=matC, beta=1)
    t2 = time.clock()
    C = A.dot(B.T) + C.T
    t3 = time.clock()

    print("")
    print("test syrkx:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #test results
    res = matC.device2host()

    #only upper triangle is computed
    res = np.triu(res).flatten()
    C = np.triu(C).flatten()

    M = np.argmax(np.abs(res - C))
    d = 5
    if (0 < np.abs(C.item(M))):
        d = 10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


def test_geam():

    #function geam
    #testing: C=A.B
    #A,B matrices

    #generating random matrices and associated carma_obj
    matA = obj_Float2D(c, dims=np.array([sizen, sizek], dtype=np.int64))
    matB = obj_Float2D(c, dims=np.array([sizen, sizek], dtype=np.int64))

    matA.random(time.clock() * 10**6)
    matB.random(time.clock() * 10**6)

    A = matA.device2host()
    B = matB.device2host()

    #matrices product
    t1 = time.clock()
    C = A + B
    t2 = time.clock()
    matC = matA.geam(matB, beta=1)
    t3 = time.clock()

    print("")
    print("test geam:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #testing result

    npt.assert_array_almost_equal(C, matC.device2host(), dec)
