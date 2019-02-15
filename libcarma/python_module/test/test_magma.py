import numpy as np
import carmaWrap as ch
import numpy.testing as npt
import time

m = 256
n = 256
min_mn = min(m, n)

dec = 4
prec = 10**(-dec)

c = ch.context.get_instance()

print("precision: ", prec)


def test_svd():

    mat = np.random.rand(m, n).astype(np.float32)

    h_mat = ch.host_obj_float(mat, ch.MA_PAGELOCK)
    h_eig = ch.host_obj_float(np.zeros([min_mn], dtype=np.float32), ch.MA_PAGELOCK)
    h_U = ch.host_obj_float(np.zeros((m, m), dtype=np.float32), ch.MA_PAGELOCK)
    h_VT = ch.host_obj_float(np.zeros((n, n), dtype=np.float32), ch.MA_PAGELOCK)

    npt.assert_array_equal(mat, np.array(h_mat))

    ch.svd_cpu_float(h_mat, h_eig, h_U, h_VT)

    # expected: U.S.V=mat
    #     U = np.array(h_U)
    #     S = np.zeros((m, n), dtype=np.float32)
    #     S = np.diag(h_eig)
    #     VT = np.array(h_VT)

    #     res = np.dot(U, S)
    #     res = np.dot(res, VT.T)

    #     iErr = np.argmax(np.abs(res - mat))
    #     err = np.abs(np.abs(res.item(iErr) - mat.item(iErr)))
    #     if (mat.item(iErr) != 0):
    #         err = err / mat.item(iErr)

    #     npt.assert_almost_equal(err, 0., decimal=dec)
    #     print("")
    #     print(err)

    p_U, p_S, p_V = np.linalg.svd(mat)
    npt.assert_array_almost_equal(np.array(h_eig), p_S, decimal=dec)


def test_getri_cpu():

    a = np.random.rand(m, m).astype(np.float32)
    a = np.dot(a, a.T)

    a += np.identity(m, dtype=np.float32)

    h_mat = ch.host_obj_float(a, ch.MA_PAGELOCK)
    h_mat2 = ch.host_obj_float(a, ch.MA_PAGELOCK)

    ch.getri_cpu_float(h_mat2)

    a_1 = np.array(h_mat2)
    #expected: a.a_1=Id
    res = np.dot(a, a_1)

    err = np.amax(np.abs(res - np.identity(m)))

    npt.assert_almost_equal(err, 0., decimal=dec)


def test_potri_cpu():

    a = np.random.rand(m, m).astype(np.float32)
    a = np.dot(a, a.T)
    a += np.identity(m, dtype=np.float32)

    h_mat = ch.host_obj_float(a, ch.MA_PAGELOCK)
    h_mat2 = ch.host_obj_float(a, ch.MA_PAGELOCK)

    ch.potri_cpu_float(h_mat2)

    a_1 = np.array(h_mat2)
    #expected: a.a_1=Id
    res = np.dot(a, a_1)

    err = np.amax(np.abs(res - np.identity(m)))
    print("")
    print(err)

    npt.assert_almost_equal(err, 0., decimal=dec)


def test_getri_gpu():

    d_mat = ch.obj_float(c, np.zeros([m, m], dtype=np.float32))
    d_mat.random(int(time.perf_counter() * 1e3))
    a = np.array(d_mat)
    a = a.T
    a.reshape(a.T.shape)

    identity = ch.obj_float(c, np.identity(m, dtype=np.float32))

    d_res = d_mat.gemm(d_mat, 'n', 't', 1, identity, 1)

    b = np.dot(a, a.T)
    mat = np.array(d_res)

    d_mat = ch.obj_float(c, d_res)
    ch.getri_float(d_res)

    d_id = d_mat.gemm(d_res)
    res_id = np.array(d_id)

    err = np.amax(np.abs(res_id - np.identity(m)))

    npt.assert_almost_equal(err, 0., decimal=dec)


def test_potri_gpu():

    d_mat = ch.obj_float(c, np.zeros([m, m], dtype=np.int64))
    d_mat.random(int(time.perf_counter() * 1e3))
    a = np.array(d_mat)
    a = a.T
    a.reshape(a.T.shape)

    identity = ch.obj_float(c, np.identity(m, dtype=np.float32))

    d_res = d_mat.gemm(d_mat, op_a='n', op_b='t', beta=1, matC=identity)

    b = np.dot(a, a.T)
    mat = np.array(d_res)

    d_mat = ch.obj_float(c, d_res)
    ch.potri_float(d_res)

    d_id = d_mat.gemm(d_res)
    res_id = np.array(d_id)

    err = np.amax(np.abs(res_id - np.identity(m)))
    print("")
    print(err)

    npt.assert_almost_equal(err, 0, decimal=dec)


def test_syevd():

    d_mat = ch.obj_float(c, np.zeros([m, m], dtype=np.int64))
    d_U = ch.obj_float(c, np.zeros([m, m], dtype=np.int64))
    h_EV = ch.host_obj_float(np.zeros(m, dtype=np.float32))
    h_EV2 = ch.host_obj_float(np.zeros(m, dtype=np.float32))

    d_res = d_mat.gemm(d_mat, op_a='n', op_b='t')

    ch.syevd_float(d_res, h_EV, d_U)

    U = np.array(d_U).T
    Mat = np.array(d_mat).T
    EV = np.diag(h_EV)

    err = np.amax(np.abs(Mat - np.dot(np.dot(U, EV), U.T)))

    print("")

    print("out of place, compute U")
    print(err)

    npt.assert_almost_equal(err, 0., decimal=dec)

    d_res = d_mat.gemm(d_mat, op_a='n', op_b='t')
    ch.syevd_float(d_res, h_EV2, computeU=False)

    err = np.amax(np.abs(h_EV - np.array(h_EV2)))

    print("in place, U not computed")
    print(err)
    npt.assert_array_equal(h_EV, h_EV2)
