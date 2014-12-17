import numpy as np
import chakra as ch
import numpy.testing as npt
import time


dec=6
prec=10**-dec

sizem=512
sizen=1024


print ""
print "Test cublas 2"
print "Precision: ", prec

c=ch.chakra_context()

shm=np.array((2,sizem,sizen),dtype=np.int64)
Mat=ch.chakra_obj_Float2D(c,dims=shm)
Mat.random(time.clock())

shsym=np.array((2,sizem,sizem),dtype=np.int64)
MatSym=ch.chakra_obj_Float2D(c,dims=shsym)
MatSym.random(time.clock())
data_R=MatSym.device2host()
data_R=(data_R+data_R.T)/2
#F-order
MatSym.host2device(np.asfortranarray(data_R))


shx=np.array((1,sizen),dtype=np.int64)
Vectx=ch.chakra_obj_Float1D(c,dims=shx)
Vectx.random(time.clock())

shy=np.array((1,sizem),dtype=np.int64)
Vecty=ch.chakra_obj_Float1D(c,dims=shy)
Vecty.random(time.clock())

Vectx2=ch.chakra_obj_Float1D(c,dims=shy)
Vectx2.random(time.clock())



def test_gemv():

    A=Mat.device2host()
    x=Vectx.device2host()
    y=Vecty.device2host()

    y=A.dot(x)

    #F-order
    Mat.host2device(np.asfortranarray(A)) 

    Mat.gemv(Vectx,1,Vecty,0)

    d=10**np.ceil(np.log10(y))
    npt.assert_array_almost_equal(y/d,Vecty.device2host()/d,dec)

def test_ger():


    x=Vectx.device2host()
    A=Mat.device2host()
    y=Vecty.device2host()

    #F-order
    MatT=ch.chakra_obj_Float2D(c,data=np.asfortranarray(A.T))
   
    caOres=Vecty.ger(Vectx,A=MatT)

    A=np.outer(y,x)+A
    
    res=caOres.device2host()
    #order-F
    res=np.reshape(res,(res.shape[1],res.shape[0]))

    d=10**np.ceil(np.log10(A))

    npt.assert_array_almost_equal(A/d,res/d,dec)


def test_symv():

    A=MatSym.device2host()

    x2=Vectx2.device2host()
    y=Vecty.device2host()

    y=A.dot(x2)

    #F-order
    MatSym.host2device(np.asfortranarray(A))

    MatSym.symv(Vectx2,1,Vecty,0)

    d=10**np.ceil(np.log10(y))
    npt.assert_array_almost_equal(y/d,Vecty.device2host()/d,dec)

