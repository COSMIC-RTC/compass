import numpy as np
import h5py
import shesha as ao
import time
import sys, os
sys.path.append(os.environ.get('SHESHA_ROOT')+'/test/roket/tools/')
sys.path.append(os.environ.get('SHESHA_ROOT')+'/test/gamora/')
import gamora
import Dphi
import matplotlib.pyplot as plt
plt.ion()

def compute_Cerr(filename, modal=True, ctype="float"):
    """ Returns the residual error covariance matrix using GROOT from a ROKET file
    :parameter:
        filename : (string) : full path to the ROKET file
        modal : (bool) : if True (default), Cerr is returned in the Btt modal basis,
                         in the actuator basis if False
        ctype : (string) : "float" or "double"
    :return:
        Cerr : (np.ndarray(dim=2, dtype=np.float32)) : residual error covariance matrix
    """

    f = h5py.File(filename,'r')
    Lambda_tar = f.attrs["target.Lambda"][0]
    Lambda_wfs = f.attrs["wfs.Lambda"]
    dt = f.attrs["ittime"]
    gain = f.attrs["gain"]
    wxpos = f.attrs["wfs.xpos"][0]
    wypos = f.attrs["wfs.ypos"][0]
    r0 = f.attrs["r0"] * (Lambda_tar/Lambda_wfs)**(6./5.)
    RASC = 180./np.pi * 3600.
    xpos = f["dm.xpos"][:]
    ypos = f["dm.ypos"][:]
    p2m = f.attrs["tel_diam"] / f.attrs["pupdiam"]
    pupshape = long(2 ** np.ceil(np.log2(f.attrs["pupdiam"]) + 1))
    xactu = (xpos - pupshape/2) * p2m
    yactu = (ypos - pupshape/2) * p2m
    H = f.attrs["atm.alt"]
    L0 = f.attrs["L0"]
    speed = f.attrs["windspeed"]
    theta = f.attrs["winddir"]*np.pi/180.
    frac = f.attrs["frac"]

    Htheta = np.linalg.norm([wxpos,wypos])/RASC*H
    vdt = speed*dt/gain
    angleht = np.arctan2(wypos,wxpos)
    fc = 1/(2*(xactu[1] - xactu[0]))
    scale = (1/r0)**(5/3.) * frac * (Lambda_tar/(2*np.pi))**2
    Nact = f["Nact"][:]
    Nact = np.linalg.inv(Nact)
    P = f["P"][:]
    Btt = f["Btt"][:]
    Tf = Btt[:-2,:-2].dot(P[:-2,:-2])
    IF, T = gamora.get_IF(filename)
    IF = IF.T
    T = T.T
    N = IF.shape[0]
    deltaTT = T.T.dot(T)/N
    deltaF = IF.T.dot(T)/N
    pzt2tt = np.linalg.inv(deltaTT).dot(deltaF.T)

    if(ctype == "float"):
        groot = ao.groot_init(Nact.shape[0], int(f.attrs["nscreens"]), angleht, fc, vdt.astype(np.float32),\
                                Htheta.astype(np.float32), f.attrs["L0"], theta,
                                scale.astype(np.float32), xactu.astype(np.float32),
                                yactu.astype(np.float32), pzt2tt.astype(np.float32),
                                Tf.astype(np.float32), Nact.astype(np.float32))
    elif(ctype == "double"):
        groot = ao.groot_initD(Nact.shape[0], int(f.attrs["nscreens"]), angleht, fc, vdt.astype(np.float64),\
                                Htheta.astype(np.float64), f.attrs["L0"].astype(np.float64), theta.astype(np.float64),
                                scale.astype(np.float64), xactu.astype(np.float64),
                                yactu.astype(np.float64), pzt2tt.astype(np.float64),
                                Tf.astype(np.float64), Nact.astype(np.float64))
    else:
        raise TypeError("Unknown ctype : must be float or double")

    groot.compute_Cerr()
    Cerr = groot.get_Cerr()
    cov_err_groot = np.zeros((Nact.shape[0]+2,Nact.shape[0]+2))
    cov_err_groot[:-2,:-2] = Cerr
    cov_err_groot[-2:,-2:] = groot.get_TTcomp()
    if (modal):
        cov_err_groot = P.dot(cov_err_groot).dot(P.T)

    f.close()
    return cov_err_groot

def compute_Cerr_cpu(filename, modal=True):
    """ Returns the residual error covariance matrix using CPU version of GROOT
    from a ROKET file
    :parameter:
        filename : (string) : full path to the ROKET file
        modal : (bool) : if True (default), Cerr is returned in the Btt modal basis,
                         in the actuator basis if False
    :return:
        Cerr : (np.ndarray(dim=2, dtype=np.float32)) : residual error covariance matrix
    """
    f = h5py.File(filename,'r')

    tabx, taby = Dphi.tabulateIj0()
    Lambda_tar = f.attrs["target.Lambda"][0]
    Lambda_wfs = f.attrs["wfs.Lambda"]
    dt = f.attrs["ittime"]
    gain = f.attrs["gain"]
    wxpos = f.attrs["wfs.xpos"][0]
    wypos = f.attrs["wfs.ypos"][0]
    r0 = f.attrs["r0"] * (Lambda_tar/Lambda_wfs)**(6./5.)
    RASC = 180./np.pi * 3600.
    xpos = f["dm.xpos"][:]
    ypos = f["dm.ypos"][:]
    p2m = f.attrs["tel_diam"] / f.attrs["pupdiam"]
    pupshape = long(2 ** np.ceil(np.log2(f.attrs["pupdiam"]) + 1))
    xactu = (xpos - pupshape/2) * p2m
    yactu = (ypos - pupshape/2) * p2m
    Ccov = np.zeros((xpos.size,xpos.size))
    Caniso = np.zeros((xpos.size,xpos.size))
    Cbp = np.zeros((xpos.size,xpos.size))
    xx = np.tile(xactu,(xactu.shape[0],1))
    yy = np.tile(yactu,(yactu.shape[0],1))
    xij = xx - xx.T
    yij = yy - yy.T

    for l in range(f.attrs["nscreens"]):
        H = f.attrs["atm.alt"][l]
        L0 = f.attrs["L0"][l]
        speed = f.attrs["windspeed"][l]
        theta = f.attrs["winddir"][l]*np.pi/180.
        frac = f.attrs["frac"][l]

        Htheta = np.linalg.norm([wxpos,wypos])/RASC*H
        vdt = speed*dt/gain
        # Covariance matrices models on actuators space
        M = np.zeros((xpos.size,xpos.size))
        Mvdt = M.copy()
        Mht = M.copy()
        Mhvdt = M.copy()
        angleht = np.arctan2(wypos,wxpos)
        fc = xactu[1] - xactu[0]

        M = np.linalg.norm([xij,yij],axis=0)
        Mvdt = np.linalg.norm([xij - vdt*np.cos(theta), yij - vdt*np.sin(theta)],axis=0)
        Mht = np.linalg.norm([xij - Htheta*np.cos(angleht), yij - Htheta*np.sin(angleht)],axis=0)
        Mhvdt = np.linalg.norm([xij - vdt*np.cos(theta) - Htheta*np.cos(angleht), yij - vdt*np.sin(theta) - Htheta*np.sin(angleht)],axis=0)

        Ccov +=  0.5 * (Dphi.dphi_lowpass(Mhvdt,fc,L0,tabx,taby) - Dphi.dphi_lowpass(Mht,fc,L0,tabx,taby) \
                - Dphi.dphi_lowpass(Mvdt,fc,L0,tabx,taby) + Dphi.dphi_lowpass(M,fc,L0,tabx,taby)) * (1./r0)**(5./3.) * frac

        Caniso += 0.5*(Dphi.dphi_lowpass(Mht,fc,L0,tabx,taby) - Dphi.dphi_lowpass(M,fc,L0,tabx,taby)) * (1./r0)**(5./3.) * frac
        Cbp += 0.5*(Dphi.dphi_lowpass(Mvdt,fc,L0,tabx,taby) - Dphi.dphi_lowpass(M,fc,L0,tabx,taby)) * (1./r0)**(5./3.) * frac

    Sp = (Lambda_tar/(2*np.pi))**2
    Ctt = (Caniso+Caniso.T)*Sp
    Ctt += ((Cbp+Cbp.T)*Sp)
    Ctt += (Ccov*Sp)

    P = f["P"][:]
    Btt = f["Btt"][:]
    Tf = Btt[:-2,:-2].dot(P[:-2,:-2])

    IF, T = gamora.get_IF(filename)
    IF = IF.T
    T = T.T
    N = IF.shape[0]
    deltaTT = T.T.dot(T)/N
    deltaF = IF.T.dot(T)/N
    pzt2tt = np.linalg.inv(deltaTT).dot(deltaF.T)

    Nact = f["Nact"][:]
    N1 = np.linalg.inv(Nact)
    Ctt = N1.dot(Ctt).dot(N1)
    ttcomp = pzt2tt.dot(Ctt).dot(pzt2tt.T)
    Ctt = Tf.dot(Ctt).dot(Tf.T)
    cov_err = np.zeros((Ctt.shape[0]+2,Ctt.shape[0]+2))
    cov_err[:-2,:-2] = Ctt
    cov_err[-2:,-2:] = ttcomp
    if(modal):
        cov_err = P.dot(cov_err).dot(P.T)
    f.close()

    return cov_err

def compare_GPU_vs_CPU(filename):
    """ Compare results of GROOT vs its CPU version in terms of execution time
    and precision on the PSF renconstruction
    :parameter:
        filename : (string) : full path to the ROKET file

    """
    tic = time.time()
    cov_err_gpu_s = compute_Cerr(filename)
    tac = time.time()
    gpu_time_s = tac - tic

    tic = time.time()
    cov_err_gpu_d = compute_Cerr(filename, ctype="double")
    tac = time.time()
    gpu_time_d = tac - tic


    tic = time.time()
    cov_err_cpu = compute_Cerr_cpu(filename)
    tac = time.time()
    cpu_time = tac - tic

    otftel, otf2, psf_cpu, gpu = gamora.psf_rec_Vii(filename,fitting=False,\
                                                cov=cov_err_cpu.astype(np.float32))
    otftel, otf2, psf_gpu_s, gpu = gamora.psf_rec_Vii(filename,fitting=False,\
                                                cov=cov_err_gpu_s.astype(np.float32))
    otftel, otf2, psf_gpu_d, gpu = gamora.psf_rec_Vii(filename,fitting=False,\
                                                cov=cov_err_gpu_d.astype(np.float32))

    print "-----------------------------------------"
    print "CPU time : ", cpu_time, " s "
    print "GPU time simple precision : ", gpu_time_s, " s "
    print "GPU time double precision : ", gpu_time_d, " s "
    print "Max absolute difference in PSFs simple precision : ", np.abs(psf_cpu-psf_gpu_s).max()
    print "Max absolute difference in PSFs double precision : ", np.abs(psf_cpu-psf_gpu_d).max()
    gamora.cutsPSF(filename, psf_cpu, psf_gpu_s)
    gamora.cutsPSF(filename, psf_cpu, psf_gpu_d)
