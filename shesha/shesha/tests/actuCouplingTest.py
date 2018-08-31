# Initialize widget and load MICADO ELT pupil

p_geom = wao.config.p_geom

import shesha.util.make_pupil as mkP
import shesha.util.utilities as util
import scipy.ndimage

import numpy as np
import matplotlib.pyplot as plt
plt.ion()

cent = p_geom.pupdiam / 2. + 0.5

p_tel = wao.config.p_tel

p_tel.t_spiders = 0.51
spup = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent,
                      cent).astype(np.float32)

p_tel.t_spiders = 0.
spup2 = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent,
                       cent).astype(np.float32)

spiders = spup2 - spup

(spidersID, k) = scipy.ndimage.label(spiders)
spidersi = util.pad_array(spidersID, p_geom.ssize).astype(np.float32)
pxListSpider = [np.where(spidersi == i) for i in range(1, k + 1)]

# DM positions in iPupil:
dmposx = wao.config.p_dm0._xpos - 0.5
dmposy = wao.config.p_dm0._ypos - 0.5
dmposMat = np.c_[dmposx, dmposy].T  # one actu per column

pitch = wao.config.p_dm0._pitch

plt.clf()
plt.scatter(*dmposMat)
DISCARD = np.zeros(len(dmposx), dtype=np.bool)
PAIRS = []

# For each of the k pieces of the spider
for k, pxList in enumerate(pxListSpider):
    pts = np.c_[pxList[1], pxList[0]]  # x,y coord of pixels of the spider piece
    # lineEq = [a, b]
    # Which minimizes leqst squares of aa*x + bb*y = 1
    lineEq = np.linalg.pinv(pts).dot(np.ones(pts.shape[0]))
    aa, bb = lineEq[0], lineEq[1]

    # Find any point of the fitted line.
    # For simplicity, the intercept with one of the axes x = 0 / y = 0
    if np.abs(bb) < np.abs(aa):  # near vertical
        onePoint = np.array([1 / aa, 0.])
    else:  # otherwise
        onePoint = np.array([0., 1 / bb])
        x = np.arange(2048)
        y = -aa / bb * x + 1 / bb
        plt.plot(x, y, color='C%u' % (4 + k), label='%u' % k)
        plt.plot()

    # Rotation that aligns the spider piece to the horizontal
    rotation = np.array([[-bb, aa], [-aa, -bb]]) / (aa**2 + bb**2)**.5

    # Rotated the spider mask
    rotatedPx = rotation.dot(pts.T - onePoint[:, None])
    # Min and max coordinates along the spider length - to filter actuators that are on
    # 'This' side of the pupil and not the other side
    minU, maxU = rotatedPx[0].min() - 5. * pitch, rotatedPx[0].max() + 5. * pitch

    # Rotate the actuators
    rotatedActus = rotation.dot(dmposMat - onePoint[:, None])
    selGoodSide = (rotatedActus[0] > minU) & (rotatedActus[0] < maxU)

    # Actuators below this piece of spider
    selDiscard = (np.abs(rotatedActus[1]) < 0.5 * pitch) & selGoodSide
    DISCARD |= selDiscard

    # Actuator 'near' this piece of spider
    selPairable = (np.abs(rotatedActus[1]) > 0.5 * pitch) & \
                    (np.abs(rotatedActus[1]) < 1.5 * pitch) & \
                    selGoodSide

    pairableIdx = np.where(selPairable)[0]  # Indices of these actuators
    uCoord = rotatedActus[0,
                          selPairable]  # Their linear coord along the spider major axis

    order = np.sort(uCoord)  # Sort by linear coordinate
    orderIdx = pairableIdx[np.argsort(uCoord)]  # And keep track of original indexes

    for i in range(0, len(order) - 1):
        # Check if next actu in sorted order is very close
        # Some lonely actuators may be hanging in this list
        if np.abs(order[i] - order[i + 1]) < .2 * pitch:
            PAIRS += [(orderIdx[i], orderIdx[i + 1])]

    plt.scatter(*dmposMat[:, selDiscard], color='C1')
    plt.scatter(*dmposMat[:, selPairable], color='C2')
    plt.legend(loc=0)

for (p, q) in PAIRS:
    plt.scatter(*dmposMat[:, [p, q]])

print('To discard: %u actu' % np.sum(DISCARD))
print('%u pairs to slave' % len(PAIRS))