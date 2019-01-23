import numpy as np
import numpy.fft as npf
import numpy.testing as npt
import carmaWrap as ch
import time

c = ch.context.get_instance()

sizex = 256
sizey = 512
sizez = 10

dec = 5
prec = 10**-dec


def test_fft_C2C():

    #testing FFT C2C

    #matrices sizes
    m = sizex
    n = sizey
    nElem = m * n

    #generating random carma_obj
    C1 = ch.obj_ComplexD2D(c, dims=np.array((m, n), dtype=np.int64))
    C1.random(time.clock() * 10**6)
    C2 = ch.obj_ComplexD2D(c, dims=np.array((m, n), dtype=np.int64))

    #matrix associated to carma_obj
    C1_data = C1.device2host()
    #switch from data layout
    C1_data = np.reshape(C1_data.flatten("F"), (m, n), order="C")

    #FFT
    t1 = time.clock()
    cpu_F = npf.fft2(C1_data)
    t2 = time.clock()
    C1.fft(dest=C2)
    t3 = time.clock()

    gpu_F = C2.device2host()
    #switch from data layout for testing
    data_gpu = gpu_F.flatten("F")
    data = cpu_F.flatten("C")

    #test results
    err = np.abs(data - data_gpu)
    MR = np.argmax(err.real)
    MI = np.argmax(err.imag)

    dR = 1
    if (0 < np.abs(data.item(MR).real)):
        dR = 10**np.ceil(np.log10(np.abs(data.item(MR).real)))
    dI = 1
    if (0 < np.abs(data.item(MI).imag)):
        dI = 10**np.ceil(np.log10(np.abs(data.item(MI).imag)))

    print("")
    print("Test FFT forward C2C")
    print("Precision: ", prec)
    print("Execution time:")
    print("Python: ", t2 - t1)
    print("Carma : ", t3 - t2)

    npt.assert_almost_equal(
            data.item(MR).real / dR,
            data_gpu.item(MR).real / dR, decimal=dec)
    npt.assert_almost_equal(
            data.item(MI).imag / dI,
            data_gpu.item(MI).imag / dI, decimal=dec)

    #IFFT
    t1 = time.clock()
    cpu_B = npf.ifft2(cpu_F)
    t2 = time.clock()
    C2.fft(dest=C1, direction=-1)
    t3 = time.clock()

    gpu_B = C1.device2host() / nElem

    #switch layout
    data = cpu_B.flatten("C")
    data_gpu = gpu_B.flatten("F")

    #test results
    err = data - data_gpu
    MR = np.argmax(np.abs(err.real))
    MI = np.argmax(np.abs(err.imag))

    dR = 1
    if (0 < np.abs(data.item(MR).real)):
        dR = 10**np.ceil(np.log10(np.abs(data.item(MR).real)))
    dI = 1
    if (0 < np.abs(data.item(MI).imag)):
        dI = 10**np.ceil(np.log10(np.abs(data.item(MI).imag)))

    print("")
    print("Test FFT backward C2C")
    print("Precision: ", prec)
    print("Execution time:")
    print("Python: ", t2 - t1)
    print("Carma : ", t3 - t2)

    npt.assert_almost_equal(
            data.item(MR).real / dR,
            data_gpu.item(MR).real / dR, decimal=dec)
    npt.assert_almost_equal(
            data.item(MI).imag / dI,
            data_gpu.item(MI).imag / dI, decimal=dec)

    npt.assert_array_almost_equal(C1_data, cpu_B, decimal=dec)


def test_fft_R2C_C2R():
    #testing FFT R2C and C2R

    #matrices sizes
    m = sizex
    n = sizey

    nElem = m * n
    nc = n / 2 + 1
    ncElem = (n + 1) / 2 * m

    #generating random carma_obj
    R1 = ch.obj_Double2D(c, dims=np.array((m, n), dtype=np.int64))
    R1.random(time.clock() * 10**6)

    #matrix associated to carma_obj
    R1_data = R1.device2host()
    #switch from data layout
    R2_data = np.reshape(R1_data.flatten("F"), (m, n), order="C")

    #FFT R2C
    t1 = time.clock()
    cpu_R2C = npf.rfft2(R2_data)
    t2 = time.clock()
    C1 = R1.fft()
    t3 = time.clock()

    gpu_R2C = C1.device2host()

    #switch layout
    data = cpu_R2C.flatten("C")
    data_gpu = gpu_R2C.flatten("F")

    #test results
    err = np.abs(data - data_gpu)
    MR = np.argmax(err.real)
    MI = np.argmax(err.imag)

    dR = 1
    if (0 < np.abs(data.item(MR).real)):
        dR = 10**np.ceil(np.log10(np.abs(data.item(MR).real)))
    dI = 1
    if (0 < np.abs(data.item(MI)).imag):
        dI = 10**np.ceil(np.log10(np.abs(data.item(MI).imag)))

    print("")
    print("Test FFT R2C")
    print("Precision: ", prec)
    print("Execution time:")
    print("Python: ", t2 - t1)
    print("Carma : ", t3 - t2)

    npt.assert_almost_equal(
            data.item(MR).real / dR,
            data_gpu.item(MR).real / dR, decimal=dec)
    npt.assert_almost_equal(
            data.item(MI).imag / dI,
            data_gpu.item(MI).imag / dI, decimal=dec)

    #IFFT
    t1 = time.clock()
    cpu_C2R = npf.irfft2(cpu_R2C)  # ,s=(m,n))
    t2 = time.clock()
    R1 = C1.fft(C2R=True, direction=-1)
    t3 = time.clock()

    gpu_C2R = R1.device2host() / nElem

    #switch layout
    data = cpu_R2C.flatten("C")
    data_gpu = gpu_R2C.flatten("F")

    #test results
    err = np.abs(data - data_gpu)
    MR = np.argmax(err.real)

    dR = 1
    if (0 < np.abs(data.item(MR).real)):
        dR = 10**np.ceil(np.log10(np.abs(data.item(MR).real)))

    print("")
    print("Test FFT C2R")
    print("Precision: ", prec)
    print("Execution time:")
    print("Python: ", t2 - t1)
    print("Carma : ", t3 - t2)

    npt.assert_almost_equal(
            data.item(MR).real / dR,
            data_gpu.item(MR).real / dR, decimal=dec)
    npt.assert_array_almost_equal(R2_data, cpu_C2R, decimal=dec)


#this test fails:
#   for each plan, the data sould be reorganized
#   to fit the row major requirement
#   as for the 2 test before
def Ntest_fft_multi():

    m = sizex / 2
    n = sizey / 2
    l = sizez

    nElem = m * n

    C1 = ch.obj_ComplexD3D(c, dims=np.array((m, n, l), dtype=np.int64))
    C2 = ch.obj_ComplexD3D(c, dims=np.array((m, n, l), dtype=np.int64))
    C3 = ch.obj_ComplexD3D(c, dims=np.array((m, n, l), dtype=np.int64))

    cpu_F = np.ones((m, n, l), dtype=np.complex128)
    cpu_B = np.ones((m, n, l), dtype=np.complex128)

    C1.random(time.clock() * 10**6)
    R1_data = C1.device2host()

    #cufftManyPlan: l successive 2D plan ( != 3D fft)
    R1_plan = np.ones((l, m, n), dtype=np.complex128)
    for i in range(l):
        R1_plan[i, :, :] = R1_data[:, :, i]

    C1.host2device(R1_plan)

    t1 = time.clock()
    for i in range(l):
        cpu_F[:, :, i] = npf.fft2(R1_data[:, :, i])
    t2 = time.clock()
    C1.fft(dest=C2)
    t3 = time.clock()

    print("")
    print("Test FFT Multi C2C")
    print("Precision: ", prec)
    print("Execution time:")
    print("Python: ", t2 - t1)
    print("Carma : ", t3 - t2)

    gpu_F = C2.device2host()

    # rearange layout for 2D plan to match
    gpu_plan = C2.device2host().reshape((m * n * l))
    for i in range(l):
        gpu_F[:, :, i] = gpu_plan[(m * n * i):(m * n * (i + 1))].reshape((m, n))

    for i in range(l):
        npt.assert_almost_equal(gpu_F[:, :, i], cpu_F[:, :, i], decimal=dec)

    t1 = time.clock()
    for i in range(l):
        cpu_B[:, :, i] = npf.ifft2(cpu_F[:, :, i])
    t2 = time.clock()
    C2.fft(dest=C3, direction=-1)
    t3 = time.clock()

    print("")
    print("Test FFT Multi backward C2C")
    print("Precision: ", prec)
    print("Execution time:")
    print("Python: ", t2 - t1)
    print("Carma : ", t3 - t2)

    gpu_B = C3.device2host() / nElem

    npt.assert_almost_equal(gpu_B, C1.device2host(), decimal=dec)
    npt.assert_almost_equal(cpu_B, R1_data, decimal=dec)
