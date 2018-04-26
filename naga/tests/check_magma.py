import numpy as np
import naga as ch
import numpy.testing as npt
import time

from naga.context import context
from naga.obj import obj_Float1D, obj_Float2D
from naga.host_obj import host_obj_Float1D, host_obj_Float2D
from naga.magma import svd_host_Float, getri_host_Float, potri_host_Float, getri_Float, potri_Float, syevd_Float

m = 256
n = 256
min_mn = min(m, n)

dec = 4
prec = 10**(-dec)

c = context()

print("precision: ", prec)


def test_svd():

    a = np.random.rand(m, n).astype(np.float32)

    h_mat = host_obj_Float2D(data=a, mallocType="pagelock")
    h_eig = host_obj_Float1D(data=np.zeros([min_mn], dtype=np.float32),
                             mallocType="pagelock")
    h_U = host_obj_Float2D(data=np.zeros((m, m), dtype=np.float32),
                           mallocType="pagelock")
    h_VT = host_obj_Float2D(data=np.zeros((n, n), dtype=np.float32),
                            mallocType="pagelock")

    Mat = a

    npt.assert_array_equal(a, h_mat.getData())

    svd_host_Float(h_mat, h_eig, h_U, h_VT)

    # expected: U.S.V=Mat
    U = h_U.getData()
    S = np.zeros((m, n), dtype=np.float32)
    S[np.diag_indices_from(S[:min_mn, :min_mn])] = h_eig.getData()
    V = h_VT.getData()

    res = np.dot(U, S)
    res = np.dot(res, V)

    iErr = np.argmax(np.abs(res - Mat))
    err = np.abs(np.abs(res.item(iErr) - Mat.item(iErr)))
    if (Mat.item(iErr) != 0):
        err = err / Mat.item(iErr)

    npt.assert_almost_equal(err, 0., decimal=dec)
    print("")
    print(err)

    p_U, p_S, p_V = np.linalg.svd(a)
    npt.assert_array_almost_equal(S.diagonal(), p_S, decimal=dec)


def test_getri_cpu():

    a = np.random.rand(m, m).astype(np.float32)
    a = np.dot(a, a.T)

    a = +np.identity(m, dtype=np.float32)

    h_mat = host_obj_Float2D(data=a, mallocType="pagelock")
    h_mat2 = host_obj_Float2D(obj=h_mat, mallocType="pagelock")

    getri_host_Float(h_mat2)

    a_1 = h_mat2.getData()
    #expected: a.a_1=Id
    res = np.dot(a, a_1)

    err = np.amax(np.abs(res - np.identity(m)))
    print("")
    print(err)

    npt.assert_almost_equal(err, 0., decimal=dec)


def test_potri_cpu():

    a = np.random.rand(m, m).astype(np.float32)
    a = np.dot(a, a.T)
    a = +np.identity(m, dtype=np.float32)

    h_mat = host_obj_Float2D(data=a, mallocType="pagelock")
    h_mat2 = host_obj_Float2D(obj=h_mat, mallocType="pagelock")

    potri_host_Float(h_mat2)

    a_1 = h_mat2.getData()
    #expected: a.a_1=Id
    res = np.dot(a, a_1)

    err = np.amax(np.abs(res - np.identity(m)))
    print("")
    print(err)

    npt.assert_almost_equal(err, 0., decimal=dec)


def test_getri_gpu():

    d_mat = obj_Float2D(c, dims=np.array([m, m], dtype=np.int64))
    d_mat.random(time.clock() * 10**6)
    a = d_mat.device2host()
    a = a.T
    a.reshape(a.T.shape)

    identity = obj_Float2D(c, data=np.identity(m, dtype=np.float32))

    d_res = d_mat.gemm(d_mat, opA='n', opB='t', beta=1, C=identity)

    b = np.dot(a, a.T)
    mat = d_res.device2host()

    d_mat = obj_Float2D(obj=d_res)
    getri_Float(d_res)

    d_id = d_mat.gemm(d_res)
    res_id = d_id.device2host()

    err = np.amax(np.abs(res_id - np.identity(m)))
    print("")
    print(err)

    npt.assert_almost_equal(err, 0., decimal=dec)


def test_potri_gpu():

    d_mat = obj_Float2D(c, dims=np.array([m, m], dtype=np.int64))
    d_mat.random(time.clock() * 10**6)
    a = d_mat.device2host()
    a = a.T
    a.reshape(a.T.shape)

    identity = obj_Float2D(c, data=np.identity(m, dtype=np.float32))

    d_res = d_mat.gemm(d_mat, opA='n', opB='t', beta=1, C=identity)

    b = np.dot(a, a.T)
    mat = d_res.device2host()

    d_mat = obj_Float2D(obj=d_res)
    potri_Float(d_res)

    d_id = d_mat.gemm(d_res)
    res_id = d_id.device2host()

    err = np.amax(np.abs(res_id - np.identity(m)))
    print("")
    print(err)

    npt.assert_almost_equal(err, 0, decimal=dec)


def test_syevd():

    d_mat = obj_Float2D(c, dims=np.array([m, m], dtype=np.int64))
    d_U = obj_Float2D(c, dims=np.array([m, m], dtype=np.int64))
    h_EV = np.zeros(m, dtype=np.float32)
    h_EV2 = np.zeros(m, dtype=np.float32)

    d_res = d_mat.gemm(d_mat, opA='n', opB='t')

    syevd_Float(d_res, h_EV, d_U)

    U = d_U.device2host().T
    Mat = d_mat.device2host().T
    EV = np.diag(h_EV)

    err = np.amax(np.abs(Mat - np.dot(np.dot(U, EV), U.T)))

    print("")

    print("out of place, compute U")
    print(err)

    npt.assert_almost_equal(err, 0., decimal=dec)

    d_res = d_mat.gemm(d_mat, opA='n', opB='t')
    syevd_Float(d_res, h_EV2, computeU=False)

    err = np.amax(np.abs(h_EV - h_EV2))

    print("in place, U not computed")
    print(err)
    npt.assert_array_equal(h_EV, h_EV2)