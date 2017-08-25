import numpy as np


def create_interp_mat(dimx: int, dimy: int):
    """TODO doc

    :parameters:
        dimx: (int) :

        dimy: (int) :
    """
    n = max(dimx, dimy)
    tmp2, tmp1 = np.indices((n, n))
    tmp1 = tmp1[:dimx, :dimy] - (dimx // 2)
    tmp2 = tmp2[:dimx, :dimy] - (dimy // 2)

    tmp = np.zeros((tmp1.size, 6), np.int32)
    tmp[:, 0] = (tmp1**2).flatten()
    tmp[:, 1] = (tmp2**2).flatten()
    tmp[:, 2] = (tmp1 * tmp2).flatten()
    tmp[:, 3] = tmp1.flatten()
    tmp[:, 4] = tmp2.flatten()
    tmp[:, 5] = 1

    return np.dot(np.linalg.inv(np.dot(tmp.T, tmp)), tmp.T).T


def centroid_gain(E, F):
    """ Returns the mean centroid gain
    :parameters:
        E : (np.array(ndim=2,dtype=np.float32)) : measurements from WFS
        F : (np.array(ndim=2,dtype=np.float32)) : geometric measurements
    :return:
        cgain : (float) : mean centroid gain between the sets of WFS measurements and geometric ones
    """
    cgains = np.zeros(E.shape[1])
    for k in range(E.shape[1]):
        cgains[k] = np.polyfit(E[:, k], F[:, k], 1)[0]

    return np.mean(cgains)
