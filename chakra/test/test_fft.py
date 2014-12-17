import numpy as np
import numpy.fft as npf
import numpy.testing as npt
import chakra as ch
import time

c=ch.chakra_context()

sizex=1024
sizey=512
sizez=10


dec=10
prec=10**-dec


def test_fft_C2C():

    m=sizex
    n=sizey

    nElem=m*n

    C1=ch.chakra_obj_ComplexD2D(c,dims=np.array((2,m,n)))
    C2=ch.chakra_obj_ComplexD2D(c,dims=np.array((2,m,n)))

    C1.random(time.clock()*10**6)
    C1_data=C1.device2host()


    t1=time.clock()
    cpu_F=npf.fft2(C1_data)
    t2=time.clock()
    C1.fft(C2)
    t3=time.clock()

    gpu_F=C2.device2host()

  
    data=np.reshape(cpu_F,(n*m))
    data_gpu=np.reshape(gpu_F,(n*m))

    err=np.abs(data-data_gpu)
    MR=np.argmax(err.real)
    MI=np.argmax(err.imag)

    print ""
    print "Test FFT forward C2C"
    print "Precision: ", prec
    print "Execution time:"
    print "Python: ",t2-t1
    print "Carma : ",t3-t2


    npt.assert_approx_equal(data[MR].real,data_gpu[MR].real,dec)
    npt.assert_approx_equal(data[MI].imag,data_gpu[MI].imag,dec)
    npt.assert_array_almost_equal(gpu_F,cpu_F,dec)


    t1=time.clock()
    cpu_B=npf.ifft2(cpu_F)
    t2=time.clock()
    C2.fft(C1,direction=-1)
    t3=time.clock()

    gpu_B=C1.device2host()

  
    data=np.reshape(cpu_B,(nElem))
    data_gpu=np.reshape(gpu_B,(nElem))
    err=data-data_gpu/nElem
    MR=np.argmax(np.abs(err.real))
    MI=np.argmax(np.abs(err.imag))

    print ""
    print "Test FFT backward C2C"
    print "Precision: ", prec
    print "Execution time:"
    print "Python: ",t2-t1
    print "Carma : ",t3-t2
   

    npt.assert_approx_equal(data[MR].real,data_gpu[MR].real/nElem,dec)
    npt.assert_approx_equal(data[MR].real,data_gpu[MR].real/nElem,dec)
    npt.assert_array_almost_equal(C1_data,cpu_B,dec)
    npt.assert_array_almost_equal(C1_data,gpu_B/nElem,dec)



def test_fft_R2C_C2R():

    m=sizex
    n=sizey

    nElem=m*n
    nc=n/2+1
    ncElem=(n+1)/2*m

    R1=ch.chakra_obj_Double2D(c,dims=np.array((2,m,n)))
    C1=ch.chakra_obj_ComplexD2D(c,dims=np.array((2,m,nc)))

    R1.random(time.clock()*10**6)
    R1_data=R1.device2host()

    t1=time.clock()
    cpu_R2C=npf.rfft2(R1_data)
    t2=time.clock()
    R1.fft(C1)
    t3=time.clock()


    gpu_R2C=C1.device2host()

    sh=C1.get_Dims()
  
  
    data=np.reshape(cpu_R2C[:,:nc],nc*m)
    data_gpu=np.reshape(gpu_R2C[:,:nc],(nc*m))

    err=np.abs(data[:ncElem]-data_gpu[:ncElem])
    MR=np.argmax(err.real)

    print ""
    print "Test FFT R2C"
    print "Precision: ", prec
    print "Execution time:"
    print "Python: ",t2-t1
    print "Carma : ",t3-t2


    npt.assert_approx_equal(data[MR].real,data_gpu[MR].real,dec)
    npt.assert_array_almost_equal(data,data_gpu,dec)


    t1=time.clock()
    cpu_C2R=npf.irfft2(cpu_R2C,s=(m,n))
    t2=time.clock()
    C1.fft(R1,direction=-1)
    t3=time.clock()

    gpu_C2R=R1.device2host()/nElem

    err=np.abs(cpu_C2R-gpu_C2R)
    MR=np.argmax(err)
  
  
    data=np.reshape(cpu_C2R,(nElem))
    data_gpu=np.reshape(gpu_C2R,(nElem))
    

    print ""
    print "Test FFT C2R"
    print "Precision: ", prec
    print "Execution time:"
    print "Python: ",t2-t1
    print "Carma : ",t3-t2
   

    npt.assert_approx_equal(data[MR],data_gpu[MR],dec)
    npt.assert_array_almost_equal(R1_data,cpu_C2R,dec)
    npt.assert_array_almost_equal(R1_data,gpu_C2R,dec)


def test_fft_multi():

    m=sizex/2
    n=sizey/2
    l=sizez

    nElem=m*n

    C1=ch.chakra_obj_ComplexD3D(c,dims=np.array((3,m,n,l)))
    C2=ch.chakra_obj_ComplexD3D(c,dims=np.array((3,m,n,l)))
    C3=ch.chakra_obj_ComplexD3D(c,dims=np.array((3,m,n,l)))

    cpu_F=np.ones((m,n,l),dtype=np.complex128)
    cpu_B=np.ones((m,n,l),dtype=np.complex128)

    C1.random(time.clock()*10**6)
    R1_data=C1.device2host()


    #cufftManyPlan: l successive 2D plan ( != 3D fft)
    R1_plan=np.ones((l,m,n),dtype=np.complex128)
    for i in range(l):
        R1_plan[i,:,:]=R1_data[:,:,i]

    
    C1.host2device(R1_plan)

    t1=time.clock()
    for i in range(l):
        cpu_F[:,:,i]=npf.fft2(R1_data[:,:,i])
    t2=time.clock()
    C1.fft(C2)
    t3=time.clock()

    print ""
    print "Test FFT Multi C2C"
    print "Precision: ", prec
    print "Execution time:"
    print "Python: ",t2-t1
    print "Carma : ",t3-t2


    gpu_F=C2.device2host()

    # rearange layout for 2D plan to match
    gpu_plan= C2.device2host().reshape((m*n*l))
    for i in range(l):
        gpu_F[:,:,i]=gpu_plan[(m*n*i):(m*n*(i+1))].reshape((m,n))


    for i in range(l):
        npt.assert_almost_equal(gpu_F[:,:,i],cpu_F[:,:,i],dec)


    t1=time.clock()
    for i in range(l):
        cpu_B[:,:,i]=npf.ifft2(cpu_F[:,:,i])
    t2=time.clock()
    C2.fft(C3,-1)
    t3=time.clock()

    print ""
    print "Test FFT Multi backward C2C"
    print "Precision: ", prec
    print "Execution time:"
    print "Python: ",t2-t1
    print "Carma : ",t3-t2


    gpu_B=C3.device2host()/nElem


    npt.assert_almost_equal(gpu_B,C1.device2host(),dec)
    npt.assert_almost_equal(cpu_B,R1_data,dec)


#def test_fft_C2C_3D():
#
#    i=l/2
#    j=m/2
#    k=n
#
#    C1=ch.chakra_obj_ComplexD3D(c,dims=np.array((3,i,j,k)))
#    C2=ch.chakra_obj_ComplexD3D(c,dims=np.array((3,i,j,k)))
#
#    C1.random(time.clock()*10**6)
#    C1_data=C1.device2host()
#
#    t1=time.clock()
#    cpu_F=npf.fftn(C1_data)
#    t2=time.clock()
#    C1.fft(C2)
#    t3=time.clock()
#
#
#    gpu_F=C2.device2host()
#
#  
#    data=np.reshape(cpu_F,(i*j*k))
#    data_gpu=np.reshape(gpu_F,(i*j*k))
#
#    err=np.abs(data-data_gpu)
#    MR=np.argmax(err.real)
#    MI=np.argmax(err.imag)
#
#    print ""
#    print "Test FFT 3D forward C2C"
#    print "Precision: ", prec
#    print "Execution time:"
#    print "Python: ",t2-t1
#    print "Carma : ",t3-t2
#
#
#    npt.assert_approx_equal(data[MR].real,data_gpu[MR].real,dec)
#    npt.assert_approx_equal(data[MI].imag,data_gpu[MI].imag,dec)
#    #npt.assert_array_almost_equal(gpu_F,cpu_F,dec)
#
#
#    t1=time.clock()
#    cpu_B=npf.ifftn(cpu_F)
#    t2=time.clock()
#    C2.fft(C1,direction=-1)
#    t3=time.clock()
#
#    gpu_B=C1.device2host()
#
#    err=np.abs(cpu_B-gpu_B)
#    MR=np.argmax(err.real)
#  
#    data=np.reshape(cpu_B,(i*j*k))
#    data_gpu=np.reshape(gpu_B,(i*j*k))
#    
#
#    print ""
#    print "Test FFT 3D backward C2C "
#    print "Precision: ", prec
#    print "Execution time:"
#    print "Python: ",t2-t1
#    print "Carma : ",t3-t2
#   
#
#    npt.assert_approx_equal(data[MR].real,data_gpu[MR].real/(i*j*k),dec)
#    npt.assert_approx_equal(data[MI].imag,data_gpu[MI].imag/(i*j*k),dec)
#    npt.assert_array_almost_equal(C1_data,cpu_B,dec)
#    npt.assert_array_almost_equal(C1_data,gpu_B/(i*j*k),dec)




#def test_fft_R2C_C2R_3D():
#
#    i=l/2
#    j=m/2
#    k=n
#
#    print i," ",j," ",k
#
#    R1=ch.chakra_obj_Double3D(c,dims=np.array((3,i,j,k)))
#    C1=ch.chakra_obj_ComplexD3D(c,dims=np.array((3,i,j,nc)))
#
#    R1.random(time.clock()*10**6)
#    R1_data=R1.device2host()
#
#    #R1_data=np.ones((m/2,m/2,n),dtype=np.float64)
#    #R1.host2device(R1_data)
#
#    t1=time.clock()
#    cpu_R2C=npf.rfftn(R1_data)
#    t2=time.clock()
#    R1.fft(C1)
#    t3=time.clock()
#
#
#    gpu_R2C=C1.device2host()
#
#    sh=C1.get_Dims()
#  
#  
#    data=np.reshape(cpu_R2C[:,:,:nc],nc*i*j)
#    data_gpu=np.reshape(gpu_R2C[:,:,:nc],(nc*i*j))
#
#    err=np.abs(data[:i*j*nc]-data_gpu[:(i*j*nc)])
#    MR=np.argmax(err.real)
#    MI=np.argmax(err.imag)
#
#    print ""
#    print "Test FFT 3D R2C"
#    print "Precision: ", prec
#    print "Execution time:"
#    print "Python: ",t2-t1
#    print "Carma : ",t3-t2
#
#
#    npt.assert_approx_equal(data[MR].real,data_gpu[MR].real,dec)
#    npt.assert_approx_equal(data[MI].imag,data_gpu[MI].imag,dec)
#    #npt.assert_array_almost_equal(data,data_gpu,dec)
#
#
#    t1=time.clock()
#    cpu_C2R=npf.irfftn(cpu_R2C)
#    t2=time.clock()
#    C1.fft(R1,direction=-1)
#    t3=time.clock()
#
#
#    gpu_C2R=R1.device2host()/nElem
#
#    err=np.abs(cpu_C2R-gpu_C2R)
#    MR=np.argmax(err)
#  
#  
#    data=np.reshape(cpu_C2R,(i*j*k))
#    data_gpu=np.reshape(gpu_C2R,(i,j,k))
#    
#
#    print ""
#    print "Test FFT 3D C2R"
#    print "Precision: ", prec
#    print "Execution time:"
#    print "Python: ",t2-t1
#    print "Carma : ",t3-t2
#   
#
#    npt.assert_array_almost_equal(R1_data,cpu_C2R,dec)
#    print C1.get_Dims()
#    print i," ",j," ",k
#    print "py ok  ",float(i)/2,float(j)/2,float(k)/2 
#    #npt.assert_approx_equal(data[MR],data_gpu[MR],dec)#/(m*m*n)*4,dec)
#    npt.assert_array_almost_equal(R1_data,gpu_C2R/R1_data,dec)
#    npt.assert_array_almost_equal(R1_data,gpu_C2R/i*2,dec)
#    #npt.assert_array_almost_equal(R1_data,gpu_C2R,dec)
