import numpy as np
import carmaWrap as ch
import numpy.testing as npt
import time

m = 1024
n = 1024
min_mn = min(m, n)

dec = 3
prec = 10**(-dec)

c = ch.context.get_instance()

print("precision: ", prec)

# def test_float_svd():

#     mat = np.random.rand(m, n).astype(np.float32)

#     h_mat = ch.host_obj_float(mat, ch.MA_PAGELOCK)
#     h_eig = ch.host_obj_float(np.random.randn([min_mn]), ch.MA_PAGELOCK)
#     h_U = ch.host_obj_float(np.random.randn((m, m)), ch.MA_PAGELOCK)
#     h_VT = ch.host_obj_float(np.random.randn((n, n)), ch.MA_PAGELOCK)

#     npt.assert_array_equal(mat, np.array(h_mat))

#     ch.svd_cpu_float(h_mat, h_eig, h_U, h_VT)

#     # expected: U.S.V=mat
#     #     U = np.array(h_U)
#     #     S = np.random.randn((m, n))
#     #     S = np.diag(h_eig)
#     #     VT = np.array(h_VT)

#     #     res = np.dot(U, S)
#     #     res = np.dot(res, VT.T)

#     #     iErr = np.argmax(np.abs(res - mat))
#     #     err = np.abs(np.abs(res.item(iErr) - mat.item(iErr)))
#     #     if (mat.item(iErr) != 0):
#     #         err = err / mat.item(iErr)

#     #     npt.assert_almost_equal(err, 0., decimal=dec)
#     #     print("")
#     #     print(err)

#     p_U, p_S, p_V = np.linalg.svd(mat)
#     npt.assert_array_almost_equal(np.array(h_eig), p_S, decimal=dec)

# def test_float_getri_cpu():

#     a = np.random.rand(m, m).astype(np.float32)
#     a = np.dot(a, a.T)

#     a += np.identity(m)

#     h_mat = ch.host_obj_float(a, ch.MA_PAGELOCK)
#     h_mat2 = ch.host_obj_float(a, ch.MA_PAGELOCK)

#     ch.getri_cpu_float(h_mat2)

#     a_1 = np.array(h_mat2)
#     #expected: a.a_1=Id
#     res = np.dot(a, a_1)

#     err = np.amax(np.abs(res - np.identity(m)))

#     npt.assert_almost_equal(err, 0., decimal=dec)

# def test_float_potri_cpu():

#     a = np.random.rand(m, m).astype(np.float32)
#     a = np.dot(a, a.T)
#     a += np.identity(m)

#     h_mat = ch.host_obj_float(a, ch.MA_PAGELOCK)
#     h_mat2 = ch.host_obj_float(a, ch.MA_PAGELOCK)

#     ch.potri_cpu_float(h_mat2)

#     a_1 = np.array(h_mat2)
#     #expected: a.a_1=Id
#     res = np.dot(a, a_1)

#     err = np.amax(np.abs(res - np.identity(m)))
#     print("")
#     print(err)

#     npt.assert_almost_equal(err, 0., decimal=dec)

# def test_float_getri_gpu():

#     d_mat = ch.obj_float(c, np.random.randn([m, m]))
#     d_mat.random(np.int32(time.perf_counter() * 1e3))
#     a = np.array(d_mat)
#     a = a.T
#     a.reshape(a.T.shape)

#     identity = ch.obj_float(c, np.identity(m))

#     d_res = d_mat.gemm(d_mat, 'n', 't', 1, identity, 1)

#     b = np.dot(a, a.T)
#     mat = np.array(d_res)

#     d_mat = ch.obj_float(c, d_res)
#     ch.getri_float(d_res)

#     d_id = d_mat.gemm(d_res)
#     res_id = np.array(d_id)

#     err = np.amax(np.abs(res_id - np.identity(m)))

#     npt.assert_almost_equal(err, 0., decimal=dec)

# def test_float_potri_gpu():

#     d_mat = ch.obj_float(c, np.random.randn([m, m))
#     d_mat.random(np.int32(time.perf_counter() * 1e3))
#     a = np.array(d_mat)
#     a = a.T
#     a.reshape(a.T.shape)

#     identity = ch.obj_float(c, np.identity(m))

#     d_res = d_mat.gemm(d_mat, op_a='n', op_b='t', beta=1, matC=identity)

#     b = np.dot(a, a.T)
#     mat = np.array(d_res)

#     d_mat = ch.obj_float(c, d_res)
#     ch.potri_float(d_res)

#     d_id = d_mat.gemm(d_res)
#     res_id = np.array(d_id)

#     err = np.amax(np.abs(res_id - np.identity(m)))
#     print("")
#     print(err)

#     npt.assert_almost_equal(err, 0, decimal=dec)


def test_float_syevd():

    d_mat = ch.obj_float(c, np.random.randn(m, m))
    d_U = ch.obj_float(c, np.random.randn(m, m))
    d_EV = ch.obj_float(c, np.random.randn(m))
    d_EV2 = ch.obj_float(c, np.random.randn(m))

    d_res = d_mat.gemm(d_mat, op_a='n', op_b='t')

    ch.syevd_float(d_res, d_EV, d_U)

    U = np.array(d_U)
    Mat = np.array(d_res)
    EV = np.diag(d_EV)

    npt.assert_almost_equal(np.dot(np.dot(U, EV), U.T), Mat, decimal=dec - 1)

    err = np.amax(np.abs(Mat - np.dot(np.dot(U, EV), U.T)))

    print("")

    print("out of place, compute U")
    print(err)

    npt.assert_almost_equal(err, 0., decimal=dec)

    d_res = d_mat.gemm(d_mat, op_a='n', op_b='t')
    ch.syevd_float(d_res, d_EV2, computeU=False)

    err = np.amax(np.abs(np.array(d_EV) - np.array(d_EV2)))

    print("in place, U not computed")
    print(err)
    npt.assert_almost_equal(np.array(d_EV), np.array(d_EV2), decimal=dec - 2)


def test_double_syevd():

    d_mat = ch.obj_double(c, np.random.randn(m, m))
    d_U = ch.obj_double(c, np.random.randn(m, m))
    d_EV = ch.obj_double(c, np.random.randn(m))
    d_EV2 = ch.obj_double(c, np.random.randn(m))

    d_res = d_mat.gemm(d_mat, op_a='n', op_b='t')

    ch.syevd_double(d_res, d_EV, d_U)

    U = np.array(d_U)
    Mat = np.array(d_res)
    EV = np.diag(d_EV)

    npt.assert_almost_equal(np.dot(np.dot(U, EV), U.T), Mat, decimal=2 * dec - 1)

    err = np.amax(np.abs(Mat - np.dot(np.dot(U, EV), U.T)))

    print("")

    print("out of place, compute U")
    print(err)

    npt.assert_almost_equal(err, 0., decimal=dec)

    d_res = d_mat.gemm(d_mat, op_a='n', op_b='t')
    ch.syevd_double(d_res, d_EV2, computeU=False)

    err = np.amax(np.abs(np.array(d_EV) - np.array(d_EV2)))

    print("in place, U not computed")
    print(err)
    npt.assert_almost_equal(np.array(d_EV), np.array(d_EV2), decimal=2 * dec - 2)
