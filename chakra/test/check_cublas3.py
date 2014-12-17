import numpy as np
import chakra as ch
import numpy.testing as npt
import time

dec=6
prec=10**-dec

sizem=4
sizen=5
sizek=6


print ""
print "Test cublas 3"
print "precision: ",prec


c=ch.chakra_context()

def test_gemm():

    matA=ch.chakra_obj_Float2D(c,dims=np.array((2,sizem,sizek),dtype=np.int64))
    matB=ch.chakra_obj_Float2D(c,dims=np.array((2,sizek,sizen),dtype=np.int64))
    matC=ch.chakra_obj_Float2D(c,dims=np.array((2,sizem,sizen),dtype=np.int64))

    matA.random(time.clock()*10**6)
    matB.random(time.clock()*10**6)
    matC.random(time.clock()*10**6)

    A=matA.device2host()
    B=matB.device2host() 
    C=matC.device2host() 
    C_gpu=C

    #F-order
    matA.host2device(np.asfortranarray(A))
    matB.host2device(np.asfortranarray(B))
    matC.host2device(np.asfortranarray(C))


    matA.gemm(matB,alpha=1,C=matC,beta=1) 
    C=A.dot(B)+C

    res=matC.device2host()

    #F-order
    res=res.reshape(res.T.shape).T

    res=res.reshape(res.size)
    C=C.reshape(C.size)

    M=np.argmax(res-C)
    m=np.argmin(res-C)
  
  
    C=C.reshape(C.size)
    res=res.reshape(res.size)
   
    npt.assert_array_almost_equal(C,res,dec)
    npt.assert_approx_equal(C.reshape(C.size)[M],res.reshape(res.size)[M],dec)
    npt.assert_approx_equal(C.reshape(C.size)[m],res.reshape(res.size)[m],dec)


def test_symm():
    matA=ch.chakra_obj_Float2D(c,dims=np.array((2,sizek,sizek),dtype=np.int64))
    matB=ch.chakra_obj_Float2D(c,dims=np.array((2,sizek,sizen),dtype=np.int64))
    matC=ch.chakra_obj_Float2D(c,dims=np.array((2,sizek,sizen),dtype=np.int64))

    matA.random(time.clock()*10**6)
    matB.random(time.clock()*10**6)
    matC.random(time.clock()*10**6)

    A=matA.device2host()
    A=(A+A.T)/2
    B=matB.device2host() 
    C=matC.device2host() 

    #order-F 
    matA.host2device(np.asfortranarray(A))
    matB.host2device(np.asfortranarray(B))
    matC.host2device(np.asfortranarray(C))

    t1=time.clock()
    matA.symm(matB,alpha=1,C=matC,beta=1)
    t2=time.clock()    
    C=A.dot(B)+C
    t3=time.clock()

    print ""
    print "test symm:"
    print "execution time (s)"
    print "python: ",t3-t2
    print "carma : ",t2-t1

    res=matC.device2host()
    #F-order
    res=res.reshape((res.T).shape).T

    C=C.reshape(C.size)
    res=res.reshape(res.size)


    M=np.argmax(res)
    m=np.argmin(res)

    npt.assert_array_almost_equal(C,res,dec)
    npt.assert_approx_equal(C.reshape(C.size)[M],res.reshape(res.size)[M],dec)
    npt.assert_approx_equal(C.reshape(C.size)[m],res.reshape(res.size)[m],dec)


def test_dgmm():
    matA=ch.chakra_obj_Float2D(c,dims=np.array((2,sizek,sizek),dtype=np.int64))
    Vectx=ch.chakra_obj_Float1D(c,dims=np.array((1,sizek),dtype=np.int64))

    matA.random(time.clock()*10**6)
    A=matA.device2host()
    d=Vectx.device2host()

    #F-order
    matA.host2device(np.asfortranarray(A))
    Vectx.host2device(np.asfortranarray(d))

    t1=time.clock()
    matC=matA.dgmm(Vectx)
    t2=time.clock()
    C=A*d
    t3=time.clock()

    print ""
    print "test dgmm:"
    print "execution time (s)"
    print "python: ",t3-t2
    print "carma : ",t2-t1

    res=matC.device2host()

    M=np.argmax(res)
    m=np.argmin(res)

    npt.assert_approx_equal(C.reshape(C.size)[M],res.reshape(res.size)[M],dec)
    npt.assert_approx_equal(C.reshape(C.size)[m],res.reshape(res.size)[m],dec)
    


def test_syrk():
    matA=ch.chakra_obj_Float2D(c,dims=np.array((2,sizen,sizek),dtype=np.int64))
    matC=ch.chakra_obj_Float2D(c,dims=np.array((2,sizen,sizen),dtype=np.int64))

    matA.random(time.clock()*10**6)
    matC.random(time.clock()*10**6)

    A=matA.device2host()
    C=matC.device2host()
    #syrk: C matrix is symetric
    C=(C+C.T)/2

    #F-order
    matA.host2device(np.asfortranarray(A))
    matC.host2device(np.asfortranarray(C))


    t1=time.clock()
    matA.syrk(C=matC,beta=1)
    t2=time.clock()
    C=A.dot(A.T)+C
    t3=time.clock()

    print ""
    print "test syrk:"
    print "execution time (s)"
    print "python: ",t3-t2
    print "carma : ",t2-t1


    res=matC.device2host()
    
    #res=(res.T).reshape(res.shape)

    C=np.tril(C)
    res=np.tril(res)

    M=np.argmax(res-C)
    m=np.argmin(res-C)


    npt.assert_array_almost_equal(C,res,dec)
    npt.assert_approx_equal(C.reshape(C.size)[M],res.reshape(res.size)[M],dec)
    npt.assert_approx_equal(C.reshape(C.size)[m],res.reshape(res.size)[m],dec)


def test_syrkx():
    matA=ch.chakra_obj_Float2D(c,dims=np.array((2,sizen,sizek),dtype=np.int64))
    matB=ch.chakra_obj_Float2D(c,dims=np.array((2,sizen,sizek),dtype=np.int64))
    matC=ch.chakra_obj_Float2D(c,dims=np.array((2,sizen,sizen),dtype=np.int64))

    matA.random(time.clock()*10**6)
    matB.copyFrom(matA)
    matC.random(time.clock()*10**6)

    A=matA.device2host()
    B=matB.device2host()
    C=matC.device2host()


    matA.host2device(np.asfortranarray(A))
    matB.host2device(np.asfortranarray(B))
    matC.host2device(np.asfortranarray(C))

    t1=time.clock()
    matA.syrkx(matB,alpha=1,C=matC,beta=1)
    t2=time.clock()
    C=A.dot(B.T)+C.T
    t3=time.clock()

    print ""
    print "test syrkx:"
    print "execution time (s)"
    print "python: ",t3-t2
    print "carma : ",t2-t1

    res=matC.device2host()

    res=np.tril(res)
    C=np.tril(C)

    res=res.reshape(res.size)
    C=C.reshape(C.size)

    M=np.argmax(res-C)
    m=np.argmin(res-C)

    npt.assert_array_almost_equal(C,res,dec)
    npt.assert_approx_equal(C.reshape(C.size)[M],res.reshape(res.size)[M],dec)
    npt.assert_approx_equal(C.reshape(C.size)[m],res.reshape(res.size)[m],dec)


def test_geam():
    matA=ch.chakra_obj_Float2D(c,dims=np.array((2,sizen,sizek),dtype=np.int64))
    matB=ch.chakra_obj_Float2D(c,dims=np.array((2,sizen,sizek),dtype=np.int64))

    matA.random(time.clock()*10**6)
    matB.random(time.clock()*10**6)

    A=matA.device2host()
    B=matB.device2host()

    matA.host2device(np.asfortranarray(A))
    matB.host2device(np.asfortranarray(B))
    
    t1=time.clock()
    C=A+B
    t2=time.clock()
    matC=matA.geam(matB,beta=1)
    t3=time.clock()

    print ""
    print "test geam:"
    print "execution time (s)"
    print "python: ",t3-t2
    print "carma : ",t2-t1

    C=C.reshape((C.T).shape,order='F')

    npt.assert_array_almost_equal(C.T,matC.device2host(),dec)

