import numpy as np
import carmaWrap as ch
import numpy.testing as npt
import time

dec = 4
prec = 10**-dec

sizem = 128
sizen = 256
sizek = 512

seed = np.int32(time.perf_counter() * 1e3)

print("")
print("Test cublas 3")
print("precision: ", prec)

c = ch.context.get_instance()


def test_gemm():

    # function gemm
    #testing: C=A.B+C
    #A,B,C matrices

    #generating random matrices A,B,C and associated carma_obj

    # np.random.seed(seed)
    A = np.empty((sizem, sizek), dtype=np.float32)
    AT = np.empty((sizek, sizem), dtype=np.float32)
    B = np.empty((sizek, sizen), dtype=np.float32)
    BT = np.empty((sizen, sizek), dtype=np.float32)
    C = np.empty((sizem, sizen), dtype=np.float32)
    C2 = np.empty((sizem, sizen), dtype=np.float32)
    C3 = np.empty((sizem, sizen), dtype=np.float32)

    # A = A.dot(A.T)
    # B = B.dot(B.T)

    matA = ch.obj_float(c, A)
    matAT = ch.obj_float(c, AT)
    matB = ch.obj_float(c, B)
    matBT = ch.obj_float(c, BT)
    matC = ch.obj_float(c, C)
    matC2 = ch.obj_float(c, C2)
    matC3 = ch.obj_float(c, C3)

    matA.random_host(seed, 'U')
    matAT.random_host(seed, 'U')
    matB.random_host(seed + 2, 'U')
    matBT.random_host(seed + 2, 'U')
    matC.random_host(seed + 3, 'U')
    matC2.random_host(seed + 3, 'U')
    matC3.random_host(seed + 3, 'U')

    A = np.array(matA)
    AT = np.array(matAT)
    B = np.array(matB)
    BT = np.array(matBT)
    C = np.array(matC)
    C2 = np.array(matC2)
    C3 = np.array(matC3)

    alpha = 2
    beta = 1

    #matrices product
    matA.gemm(matB, 'n', 'n', alpha, matC, beta)
    matAT.gemm(matB, 't', 'n', alpha, matC2, beta)
    matAT.gemm(matBT, 't', 't', alpha, matC3, beta)
    matC4 = matA.gemm(matB, 'n', 'n', alpha)
    matC5 = matAT.gemm(matB, 't', 'n', alpha)
    matC6 = matAT.gemm(matBT, 't', 't', alpha)

    C = alpha * A.dot(B) + beta * C
    C2 = alpha * AT.T.dot(B) + beta * C2
    C3 = alpha * AT.T.dot(BT.T) + beta * C3
    C4 = alpha * A.dot(B)
    C5 = alpha * AT.T.dot(B)
    C6 = alpha * AT.T.dot(BT.T)

    npt.assert_array_almost_equal(C, np.array(matC), decimal=dec - 1)
    npt.assert_array_almost_equal(C2, np.array(matC2), decimal=dec - 1)
    npt.assert_array_almost_equal(C3, np.array(matC3), decimal=dec - 1)
    npt.assert_array_almost_equal(C4, np.array(matC4), decimal=dec - 1)
    npt.assert_array_almost_equal(C5, np.array(matC5), decimal=dec - 1)
    npt.assert_array_almost_equal(C6, np.array(matC6), decimal=dec - 1)


def test_symm():
    #function symm
    #testing: C=A.B+C
    #A ssymetric matrix, B,C matrices
    A = np.empty((sizek, sizek), dtype=np.float32)
    B = np.empty((sizek, sizen), dtype=np.float32)
    C = np.empty((sizek, sizen), dtype=np.float32)

    #generating random matrices and associated carma_obj
    matA = ch.obj_float(c, A)  #np.zeros((sizek, sizek)))
    matB = ch.obj_float(c, B)  #np.zeros((sizek, sizen)))
    matC = ch.obj_float(c, C)  #np.zeros((sizek, sizen)))

    matA.random_host(seed, 'U')
    matB.random_host(seed + 2, 'U')
    matC.random_host(seed + 3, 'U')

    #A symetric
    A = np.array(matA)
    A = (A + A.T) / 2
    matA.host2device(A)
    B = np.array(matB)
    C = np.array(matC)

    #matrices multiplication
    t1 = time.perf_counter()
    matA.symm(matB, 1, matC, 1)
    t2 = time.perf_counter()
    C = A.dot(B) + C
    t3 = time.perf_counter()

    matC2 = matA.symm(matB)
    C2 = A.dot(B)

    print("")
    print("test symm:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #test results
    res = np.array(matC)

    npt.assert_almost_equal(C, res, decimal=dec - 1)
    npt.assert_almost_equal(C2, np.array(matC2), decimal=dec - 1)

    # M = np.argmax(np.abs(res - C))
    # d = 5
    # if (0 < np.abs(C.item(M))):


# d = 10**np.ceil(np.log10(np.abs(C.item(M))))

# print(res.item(M))
# print(C.item(M))

# npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


def test_dgmm():
    #function dgmm
    #testing: C=A.d
    # C,A matrices, d vector (diagonal matrix as a vector)

    #generating random matrices and associated carma_obj
    matA = ch.obj_float(c, np.zeros((sizek, sizek), dtype=np.float32))
    Vectx = ch.obj_float(c, np.zeros((sizek), dtype=np.float32))

    matA.random_host(seed, 'U')
    Vectx.random_host(seed + 2, 'U')
    A = np.array(matA)
    x = np.array(Vectx)

    #matrices product
    t1 = time.perf_counter()
    matC = matA.dgmm(Vectx)
    t2 = time.perf_counter()
    C = A * x
    t3 = time.perf_counter()

    print("")
    print("test dgmm:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #test results
    res = np.array(matC)

    npt.assert_almost_equal(C, res, decimal=dec)

    # M = np.argmax(np.abs(res - C))
    # d = 5
    # if (0 < np.abs(C.item(M))):
    #     d = 10**np.ceil(np.log10(np.abs(C.item(M))))

    # print(res.item(M))
    # print(C.item(M))

    # npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


def test_syrk():
    #function syrk
    #testing: C=A.transpose(A)+C
    #A matrix, C symetric matrix

    #generating random matrices and associated carma_obj
    matA = ch.obj_float(c, np.zeros((sizen, sizek)))
    matC = ch.obj_float(c, np.zeros((sizen, sizen)))

    matA.random_host(seed, 'U')
    matC.random_host(seed + 2, 'U')

    A = np.array(matA)
    C = np.array(matC)
    #syrk: C matrix is symetric
    C = (C + C.T) / 2

    matC.host2device(C)

    #matrices product
    t1 = time.perf_counter()
    matA.syrk(matC=matC, beta=1)
    t2 = time.perf_counter()
    C = A.dot(A.T) + C
    t3 = time.perf_counter()

    print("")
    print("test syrk:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #test results
    res = np.array(matC)

    #only upper triangle is computed
    C = np.triu(C).flatten()
    res = np.triu(res).flatten()
    npt.assert_almost_equal(C, res, decimal=dec - 1)

    # M = np.argmax(np.abs(res - C))
    # d = 5
    # if (0 < np.abs(C.item(M))):
    #     d = 10**np.ceil(np.log10(np.abs(C.item(M))))

    # print(res.item(M))
    # print(C.item(M))

    # npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


def test_syrkx():
    #function syrkx
    #testing: C=A.transpose(B)+C
    #A matrix, C symetric matrix

    #generating random matrices and associated carma_obj
    matA = ch.obj_float(c, np.zeros((sizen, sizek)))
    matB = ch.obj_float(c, np.zeros((sizen, sizek)))
    matC = ch.obj_float(c, np.zeros((sizen, sizen)))

    matA.random_host(seed, 'U')
    matB.random_host(seed + 2, 'U')
    matC.random_host(seed + 3, 'U')

    A = np.array(matA)
    B = np.array(matB)
    C = np.array(matC)

    #C is symetric
    C = np.dot(C, C.T)

    matC.host2device(C)

    #matrices product
    t1 = time.perf_counter()
    matA.syrkx(matB, alpha=1, matC=matC, beta=1)
    t2 = time.perf_counter()
    C = A.dot(B.T) + C.T
    t3 = time.perf_counter()

    print("")
    print("test syrkx:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #test results
    res = np.array(matC)

    #only upper triangle is computed
    res = np.triu(res)
    C = np.triu(C)
    npt.assert_almost_equal(C, res, decimal=dec - 1)

    # M = np.argmax(np.abs(res - C))
    # d = 5
    # if (0 < np.abs(C.item(M))):
    #     d = 10**np.ceil(np.log10(np.abs(C.item(M))))

    # print(res.item(M))
    # print(C.item(M))

    # npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


def test_geam():

    #function geam
    #testing: C=A.B
    #A,B matrices

    #generating random matrices and associated carma_obj
    matA = ch.obj_float(c, np.empty((sizem, sizen)))
    matB = ch.obj_float(c, np.empty((sizem, sizen)))

    matA.random_host(seed, 'U')
    matB.random_host(seed + 2, 'U')

    A = np.array(matA)
    B = np.array(matB)

    #matrices product
    t1 = time.perf_counter()
    C = A + B
    t2 = time.perf_counter()
    matC = matA.geam(matB, beta=1)
    t3 = time.perf_counter()

    print("")
    print("test geam:")
    print("execution time (s)")
    print("python: ", t3 - t2)
    print("carma : ", t2 - t1)

    #testing result

    npt.assert_array_almost_equal(C, np.array(matC), dec)
