import numpy as np
import naga as ch
import numpy.testing as npt
import time

dec = 4
prec = 10**-dec

sizem = 128
sizen = 256
sizek = 512

seed = 1234  # int(time.clock()*10**6)

print("")
print("Test cublas 3")
print("precision: ", prec)

c = ch.naga_context.get_instance()


def test_gemm():

    # function gemm
    #testing: C=A.B+C
    #A,B,C matrices

    #generating random matrices A,B,C and associated carma_obj

    A = np.random.randn(sizem, sizek)
    B = np.random.randn(sizek, sizen)
    C = np.random.randn(sizem, sizen)

    # A = A.dot(A.T)
    # B = B.dot(B.T)

    matA = ch.naga_obj_float(c, A)
    matB = ch.naga_obj_float(c, B)
    matC = ch.naga_obj_float(c, C)

    matA.random(seed)
    matB.random(seed * 2)
    matC.random(seed * 3)

    A = np.array(matA)
    B = np.array(matB)
    C = np.array(matC)

    alpha = 1
    beta = 0

    #matrices product
    matA.gemm(matB, 'n', 'n', alpha, matC, beta)
    #test results
    res = np.array(matC)

    C = alpha * A.dot(B) + beta * C

    M = np.argmax(np.abs(res - C))
    d = 5
    if (0 < np.abs(C.item(M))):
        d = 10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M) / d, res.item(M) / d, decimal=dec)


'''
def test_symm():
    #function symm
    #testing: C=A.B+C
    #A ssymetric matrix, B,C matrices

    #generating random matrices and associated carma_obj
    matA=ch.naga_obj_float(c,np.zeros((sizek,sizek)))
    matB=ch.naga_obj_float(c,np.zeros((sizek,sizen)))
    matC=ch.naga_obj_float(c,np.zeros((sizek,sizen)))

    matA.random(seed)
    matB.random(seed*2)
    matC.random(seed*3)

    #A symetric
    A=np.array(matA)
    A=(A+A.T)/2
    matA.host2device(A)
    B=np.array(matB)
    C=np.array(matC)

    #matrices multiplication
    t1=time.clock()
    matC.symm(b"l", b"l", 1, matA, sizek, matB, sizek,1, sizek)
    t2=time.clock()
    C=A.dot(B)+C
    t3=time.clock()

    print("")
    print("test symm:")
    print("execution time (s)")
    print("python: ",t3-t2)
    print("carma : ",t2-t1)

    #test results
    res=np.array(matC)


    M=np.argmax(np.abs(res-C))
    d=5
    if(0<np.abs(C.item(M))):
        d=10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M)/d,res.item(M)/d,decimal=dec)


def test_dgmm():
    #function dgmm
    #testing: C=A.d
    # C,A matrices, d vector (diagonal matrix as a vector)

     #generating random matrices and associated carma_obj
    matA=ch.naga_obj_float(c,np.zeros((sizek,sizek)))
    Vectx=ch.naga_obj_Float1D(c,np.zeros((sizek)))

    matA.random(seed)
    A=np.array(matA)
    d=np.array(Vectx)

    #matrices product
    t1=time.clock()
    matC=matA.dgmm(Vectx)
    t2=time.clock()
    C=A*d
    t3=time.clock()

    print("")
    print("test dgmm:")
    print("execution time (s)")
    print("python: ",t3-t2)
    print("carma : ",t2-t1)

    #test results
    res=np.array(matC)

    M=np.argmax(np.abs(res-C))
    d=5
    if(0<np.abs(C.item(M))):
        d=10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M)/d,res.item(M)/d,decimal=dec)


def test_syrk():
    #function syrk
    #testing: C=A.transpose(A)+C
    #A matrix, C symetric matrix

     #generating random matrices and associated carma_obj
    matA=ch.naga_obj_float(c,np.zeros((sizen,sizek)))
    matC=ch.naga_obj_float(c,np.zeros((sizen,sizen)))

    matA.random(seed)
    matC.random(seed*2)

    A=np.array(matA)
    C=np.array(matC)
    #syrk: C matrix is symetric
    C=(C+C.T)/2

    matC.host2device(C)

    #matrices product
    t1=time.clock()
    matA.syrk(C=matC,beta=1)
    t2=time.clock()
    C=A.dot(A.T)+C
    t3=time.clock()

    print("")
    print("test syrk:")
    print("execution time (s)")
    print("python: ",t3-t2)
    print("carma : ",t2-t1)


    #test results
    res=np.array(matC)

    #only upper triangle is computed
    C=np.triu(C).flatten()
    res=np.triu(res).flatten()

    M=np.argmax(np.abs(res-C))
    d=5
    if(0<np.abs(C.item(M))):
        d=10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M)/d,res.item(M)/d,decimal=dec)

def test_syrkx():
    #function syrkx
    #testing: C=A.transpose(B)+C
    #A matrix, C symetric matrix

     #generating random matrices and associated carma_obj
    matA=ch.naga_obj_float(c,np.zeros((sizen,sizek)))
    matB=ch.naga_obj_float(c,np.zeros((sizen,sizek)))
    matC=ch.naga_obj_float(c,np.zeros((sizen,sizen)))

    matA.random(seed)
    matB.copyFrom(matA)
    matC.random(seed*2)

    A=np.array(matA)
    B=np.array(matB)
    C=np.array(matC)

    #C is symetric
    C=np.dot(C,C.T)

    matC.host2device(C)

    #matrices product
    t1=time.clock()
    matA.syrkx(matB,alpha=1,C=matC,beta=1)
    t2=time.clock()
    C=A.dot(B.T)+C.T
    t3=time.clock()

    print("")
    print("test syrkx:")
    print("execution time (s)")
    print("python: ",t3-t2)
    print("carma : ",t2-t1)

    #test results
    res=np.array(matC)

    #only upper triangle is computed
    res=np.triu(res).flatten()
    C=np.triu(C).flatten()

    M=np.argmax(np.abs(res-C))
    d=5
    if(0<np.abs(C.item(M))):
        d=10**np.ceil(np.log10(np.abs(C.item(M))))

    print(res.item(M))
    print(C.item(M))

    npt.assert_almost_equal(C.item(M)/d,res.item(M)/d,decimal=dec)

def test_geam():

    #function geam
    #testing: C=A.B
    #A,B matrices

    #generating random matrices and associated carma_obj
    matA=ch.naga_obj_float(c,np.zeros((sizen,sizek)))
    matB=ch.naga_obj_float(c,np.zeros((sizen,sizek)))

    matA.random(seed)
    matB.random(seed*2)

    A=np.array(matA)
    B=np.array(matB)

    #matrices product
    t1=time.clock()
    C=A+B
    t2=time.clock()
    matC=matA.geam(matB,beta=1)
    t3=time.clock()

    print("")
    print("test geam:")
    print("execution time (s)")
    print("python: ",t3-t2)
    print("carma : ",t2-t1)

    #testing result

    npt.assert_array_almost_equal(C,matCnp.array(),dec)
'''
