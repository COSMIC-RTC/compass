import numpy as np
import naga as ch
import numpy.testing as npt
import time


dec=6
prec=10**-dec

sizem=512
sizen=1024


print ""
print "Test cublas 2"
print "Precision: ", prec

c=ch.naga_context()

#generatig random naga_obj 2d
shm=np.array([sizem,sizen],dtype=np.int64)
Mat=ch.naga_obj_Float2D(c,dims=shm)
Mat.random(time.clock())

#generating random symetric naga_obj 2d
shsym=np.array([sizem,sizem],dtype=np.int64)
MatSym=ch.naga_obj_Float2D(c,dims=shsym)
MatSym.random(time.clock())
data_R=MatSym.device2host()
data_R=(data_R+data_R.T)/2
MatSym.host2device(data_R)


#generating 3 random naga_obj 1d
shx=np.array([sizen],dtype=np.int64)
Vectx=ch.naga_obj_Float1D(c,dims=shx)
Vectx.random(time.clock())

shy=np.array([sizem],dtype=np.int64)
Vecty=ch.naga_obj_Float1D(c,dims=shy)
Vecty.random(time.clock())

Vectx2=ch.naga_obj_Float1D(c,dims=shy)
Vectx2.random(time.clock())



def test_gemv():

    #function gemv
    # testing: y=A.x
    # x and y are vector, A a matrix

    A=Mat.device2host()
    x=Vectx.device2host()
    y=Vecty.device2host()

    y=A.dot(x)

    Mat.gemv(Vectx,1,Vecty,0)
    Vecty_2=Mat.gemv(Vectx)
    
    yG=Vecty.device2host()

    M=np.argmax(np.abs(y-yG))

    d=1
    if(0<np.abs(y.item(M))):
        d=10**np.ceil(np.log10(np.abs(y.item(M))))
    npt.assert_almost_equal(y.item(M)/d,yG.item(M)/d,decimal=dec)
    npt.assert_array_almost_equal(yG,Vecty_2.device2host(),decimal=dec)

def test_ger():

    # function ger
    # testing: A= x.y
    #   and  : A= x.y+ A
    # x and y are vectors, A a matrix

    x=Vectx.device2host()
    A=Mat.device2host()
    y=Vecty.device2host()

    MatT=ch.naga_obj_Float2D(obj=Mat)

    caOres=Vectx.ger(Vecty,A=MatT)
    caOresB=Vectx.ger(Vecty)

    B=np.outer(y,x)
    A=np.outer(y,x)+A
    
    res=caOres.device2host()
    resB=caOresB.device2host()

    Ma=np.argmax(np.abs(A-res))
    Mb=np.argmax(np.abs(B-resB))

    d=1
    if(0<np.abs(A.item(Ma))):
        d=10**np.ceil(np.log10(np.abs(A.item(Ma))))
    dB=1
    if(0<np.abs(B.item(Mb))): 
        dB=10**np.ceil(np.log10(np.abs(B.item(Mb))))

    npt.assert_array_almost_equal(A/d,res/d,decimal=dec)
    npt.assert_array_almost_equal(B/d,resB/d,decimal=dec)


def test_symv():

    #function symv
    # testing: y=A.x
    # x and y are vector, A a symetric matrix

    A=MatSym.device2host()

    x2=Vectx2.device2host()
    y=Vecty.device2host()

    y=A.dot(x2)

    MatSym.symv(Vectx2,1,Vecty,0)

    yG=Vecty.device2host()

    M=np.argmax(np.abs(y-yG))
    d=1
    if(0<np.abs(y.item(M))):
        d=10**np.ceil(np.log10(np.abs(y.item(M))))

    npt.assert_array_almost_equal(y/d,yG/d,decimal=dec)

