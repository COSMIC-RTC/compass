import astropy.io.fits as pfits
import time
import tools

def applyVoltGetSlopes(wao, noise=False):
    wao.rtc.applycontrol(0, wao.dms)
    for w in range(len(wao.config.p_wfss)):
        wao.wfs.sensors_trace(w, "dm", wao.tel, wao.atm, wao.dms, rst=1)
        wao.wfs.sensors_compimg(w, noise=noise)
    wao.rtc.docentroids(0)
    return wao.rtc.getcentroids(0)


def measureIMatKL(ampliVec, KL2V, Nslopes):
    iMatKL = np.zeros((KL2V.shape[1], Nslopes))
    wao.rtc.set_openloop(0, 1)  # openLoop
    st = time.time()
    for kl in range(KL2V.shape[1]):
        v = ampliVec[kl] * KL2V[:, kl:kl + 1].T.copy()
        wao.rtc.set_perturbcom(0, v)
        iMatKL[kl, :] = applyVoltGetSlopes(wao)
        print("Doing KL interaction matrix on mode: #%d\r" % kl, end=' ')
    print("Modal interaction matrix done in %3.0f seconds" % (time.time() - st))
    return iMatKL


def normalizeKL2V(KL2V):
    KL2VNorm = KL2V * 0.
    for kl in range(KL2V.shape[1]):
        klmaxVal = np.abs(KL2V[:, kl]).max()
        # Norm du max des modes a 1
        KL2VNorm[:, kl] = KL2V[:, kl] / klmaxVal
    return KL2VNorm


def cropKL2V(KL2V):
    nkl = KL2V.shape[1]
    KL2Vmod = np.zeros((KL2V.shape[0], KL2V.shape[1] - 2))
    KL2Vmod[:, 0:-2] = KL2V[:, 0:-4]
    KL2Vmod[:, -2:] = KL2V[:, -2:]
    return KL2Vmod


def plotVDm(numdm, V, size=16, fignum=False):
    """
    plotVDm(0, KL2V[:,0][:-2], size=25, fignum=15)
    """
    if(wao.config.p_dms[numdm]._j1.shape[0] != V.shape[0]):
        print("Error V should have %d dimension not %d " % (wao.config.p_dms[numdm]._j1.shape[0], V.shape[0]))
    else:
        if(fignum):
            plt.figure(fignum)
        plt.clf()
        plt.scatter(wao.config.p_dms[numdm]._j1, wao.config.p_dms[
                    numdm]._i1, c=V, marker="s", s=size**2)
        plt.colorbar()


def computeKLModesImat(pushDMMic, pushTTArcsec, KL2V, Nslopes):
    modesAmpli = np.ones(KL2V.shape[0])
    # normalisation de la valeur max des KL:
    # Amplitude des modes DM en microns et arc sec pour le TT
    KL2VNorm = normalizeKL2V(KL2V)
    NKL = KL2VNorm.shape[1]
    modesAmpli[0:NKL - 2] = pushDMMic  # DM en microns DM en microns
    modesAmpli[NKL - 2:] = pushTTArcsec  # arc sec pour le TT
    imat = measureIMatKL(modesAmpli, KL2VNorm, Nslopes)
    imat[:-2, :] /= pushDMMic
    imat[-2:, :] /= pushTTArcsec
    return imat, KL2VNorm


def computeCmatKL(DKL, KL2V, nfilt, gains, gainTT):
    nmodes = (DKL.shape[0] - nfilt)
    KL2V2Filt = np.zeros((KL2V.shape[0], KL2V.shape[1] - nfilt))
    DKLfilt = np.zeros((nmodes, DKL.shape[1]))
    # Filter the nfilt modes
    DKLfilt[0:nmodes, :] = DKL[0:nmodes, :]
    DKLfilt[-2:, :] = DKL[-2:, :]  # Concatenating the TT (last 2 columns)
    KL2V2Filt[:, 0:nmodes] = KL2V[:, 0:nmodes]
    KL2V2Filt[:, -2:] = KL2V[:, -2:]
    # Direct inversion
    Dmp = np.linalg.inv(DKLfilt.dot(DKLfilt.T)).dot(DKLfilt)
    for i in range(nmodes-2):
        Dmp[i, :] *= gains[i]
    Dmp[-2:, :] *= gainTT
    # Command matrix
    cmat = KL2V2Filt.dot(Dmp).astype(np.float32)
    return cmat

def testDM(KL2VN, nmode, Nslopes, pushDMMic, pushTTArcsec, wao, disp=False):
    ampliVec = np.ones(KL2VN.shape[0])
    NKL = KL2VN.shape[1]
    ampliVec[0:NKL - 2] = pushDMMic  # DM en microns DM en microns
    ampliVec[NKL - 2:] = pushTTArcsec
    iMatKL = np.zeros((KL2V.shape[1], Nslopes))
    wao.rtc.set_openloop(0, 0)  # openLoop
    v = ampliVec[nmode] * KL2VN[:, nmode:nmode + 1].T.copy()
    wao.rtc.set_perturbcom(0, v)
    s = applyVoltGetSlopes(wao)
    if(disp):
        plt.figure(10)
        plt.clf()
        plotVDm(0, v[0][:-2], fignum=10)
        plt.figure(11)
        plt.clf()
        tools.plpyr(s, wao.config.p_wfs0._isvalid)
    return s, v


#ampli = 1
#Nactu = sum(wao.config.p_rtc.controllers[0].nactu)
#Nslopes = wao.rtc.getCentroids(0).shape[0]
#wao.setIntegratorLaw()
#gain = 0.8
#wao.rtc.set_openloop(0, 1)  # openLoop
#decay = np.ones(Nslopes, dtype=(np.float32))
#wao.rtc.set_decayFactor(0, decay)
#mgain = np.ones(Nactu, dtype=(np.float32))
#wao.rtc.set_mgain(0, mgain)
##mcCompass = pfits.getdata("/home/fvidal/ADOPT/data/cmat39mCompass3500KL.fits").byteswap().newbyteorder()
##wao.rtc.set_cmat(0, mcCompass.copy())

#wao.rtc.set_openloop(0, 0)  # closeLoop
#KL2V = wao.returnkl2V()
#KL2V2 = KL2V.copy()
#pushDMMic = 0.01  # 10nm
#pushTTArcsec = 0.005  # 5 mas
#wao.rtc.do_centroids_ref(0)
#miKL, KL2VN =computeKLModesImat(pushDMMic, pushTTArcsec, KL2V2, Nslopes)
##KL2VN = normalizeKL2V(KL2V2)
#nfilt = 450
#cmat = computeCmatKL(miKL, KL2VN, nfilt, np.linspace(2.2,0.1,Nactu-2-nfilt), 2.); wao.rtc.set_cmat(0, cmat.astype(np.float32).copy())

#gDM = 1; gTT=1.;
#gains = np.ones(Nactu,dtype=np.float32)*gDM;
#gains[-2:]=gTT;
#wao.rtc.set_mgain(0, gains)
#wao.rtc.set_openloop(0, 0)  # openLoop


"""
KL2VN = normalizeKL2V(KL2V)
plt.ion()
plt.show()
s, v = testDM(KL2VN, 4538, Nslopes, 1., 0.1, wao, disp=True)
"""
