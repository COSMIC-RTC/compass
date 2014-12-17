import chakra as ch
import numpy as np
import numpy.testing as npt
import time

size=128
dec=5
prec=10**-dec

print "Test cublas 1"
print "precision: ", prec 

c=ch.chakra_context()
#sh2=np.array((2,size,size),dtype=np.int64)
sh2=np.ndarray((3),dtype=np.int64)

sh2[0]=2
sh2[1]=size
sh2[2]=size
caF2D=ch.chakra_obj_Float2D(c,dims=sh2)
caF2D.random(time.clock()*10**6)

sh1=np.ndarray((2),dtype=np.int64)
sh1[0]=1
sh1[1]=size*size
Vect=ch.chakra_obj_Float1D(c,dims=sh1)
Vect.random(time.clock()*10**6)


def test_imax():
    imaxC=caF2D.device2host().argmax()+1 
    imaxG=caF2D.imax()
    
    npt.assert_equal(imaxG,imaxC)

def test_imin():
    iminC=caF2D.device2host().argmin()+1 
    iminG=caF2D.imin()
    
    npt.assert_equal(iminG,iminC)

def test_asum():
    sumC=np.abs(caF2D.device2host()).sum()
    sumG=caF2D.asum()
    
    npt.assert_approx_equal(sumG,sumC, dec)

def test_nrm2():
    nC=np.linalg.norm(caF2D.device2host())
    nG=caF2D.nrm2()
    npt.assert_approx_equal(nG,nC, dec)

def test_scale():
    sC=caF2D.device2host()*10.0
    caF2D.scale(10.0)
    sG=caF2D.device2host()
    npt.assert_array_equal(sG,sC)

def test_swap():
    d_caF2D=caF2D.device2host()
    ca2=ch.chakra_obj_Float2D(c,dims=sh2)
    d_ca2=ca2.device2host()
    caF2D.swap(ca2)
    npt.assert_array_equal(d_ca2,caF2D.device2host())
    npt.assert_array_equal(d_caF2D,ca2.device2host())

def test_copy():
    ca2=ch.chakra_obj_Float2D(c,dims=sh2)
    ca2.copy(caF2D)
    npt.assert_array_equal(caF2D.device2host(),ca2.device2host())
   
def test_axpy():
    alpha=10
    beta=1.4
    ca_Vy=ch.chakra_obj_Float1D(obj=Vect)
    ca_Vy.scale(beta) 
    Vect.axpy(alpha,ca_Vy)
    
    d_ca=Vect.device2host()*(alpha+beta)
    
    npt.assert_almost_equal(ca_Vy.device2host(),d_ca,decimal=dec)


def test_dot():
    d_ca=Vect.device2host()
    dotC=np.dot(d_ca,d_ca)
    dotG=Vect.dot(Vect)

    d=10**np.ceil(np.log10(dotC))
    npt.assert_almost_equal(dotC/d,dotG/d,decimal=dec)
