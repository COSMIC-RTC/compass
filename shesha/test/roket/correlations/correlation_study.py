"""
Created on Wed Oct 5 14:28:23 2016

@author: fferreira
"""

import numpy as np
import h5py
import glob
import sys
sys.path.append('/home/fferreira/compass/shesha/test/roket/tools/')
sys.path.append('/home/fferreira/compass/shesha/test/psf_reconstruction/')
import Dphi
import roket_exploitation as rexp
import psf_rec as precs
import matplotlib.pyplot as plt
plt.ion()
import matplotlib
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)


def compute_psf(filename):
    otftel, otf, psf, gpu = precs.psf_rec_Vii(filename)
    return psf


def compute_psf_independence(filename):
    cov_err = rexp.get_coverr_independence(filename)
    otfteli, otf2i, psfi, gpu = precs.psf_rec_Vii(filenames[11],covmodes=cov_err)
    return psfi


def compute_and_compare_PSFs(filename,correction=False,synth=False):
    f = h5py.File(filename,'r')
    psf_compass = np.fft.fftshift(f["psf"][:])
    psf = compute_psf(filename)
    psfi = compute_psf_independence(filename)
    psfc = 0
    psfs = 0
    if(correction or synth):
        Caniso, Cbp, Cmodel = compute_covariance_model(filename)
#    if(correction):
        # Cmodel_filtered = filter_piston_TT(filename,Cmodel)
        # Ctt = add_TT_model(filename,Cmodel)
        # Ctt[:-2,:-2] = Cmodel_filtered
        # cov_err = rexp.get_coverr_independence(filename)
        # P = f["P"][:]
        # cov_err = P.dot(cov_err).dot(P.T) + P.dot(Ctt+Ctt.T).dot(P.T)
        # otftelc, otf2c, psfc, gpu = precs.psf_rec_Vii(filename,cov=cov_err.astype(np.float32))
    if(synth):
        Cmodel_filtered = filter_piston_TT(filename,Cmodel)
        Ctt = add_TT_model(filename,Cmodel)
        Ctt[:-2,:-2] = Cmodel_filtered
        Ctt = Ctt + Ctt.T
        Caniso_filtered = filter_piston_TT(filename,Caniso)
        tmp = add_TT_model(filename,Caniso)
        tmp[:-2,:-2] = Caniso_filtered
        Ctt += tmp
        Cbp_filtered = filter_piston_TT(filename,Cbp)
        tmp = add_TT_model(filename,Cbp)
        tmp[:-2,:-2] = Cbp_filtered
        Ctt += tmp
        contributors = ["noise","aliasing","non linearity","filtered modes"]
        cov_err = rexp.get_coverr_independence_contributors(filename,contributors)
        P = f["P"][:]
        cov_err = P.dot(cov_err).dot(P.T) + P.dot(Ctt).dot(P.T)
        otftels, otf2s, psfs, gpu = precs.psf_rec_Vii(filename,cov=cov_err.astype(np.float32))

    Lambda_tar = f.attrs["target.Lambda"][0]
    RASC = 180/np.pi*3600.
    pixsize = Lambda_tar*1e-6  / (psf.shape[0] * f.attrs["tel_diam"]/f.attrs["pupdiam"]) * RASC
    x = (np.arange(psf.shape[0]) - psf.shape[0]/2) * pixsize / (Lambda_tar*1e-6/f.attrs["tel_diam"] * RASC)
    plt.figure()
    plt.semilogy(x,psf[psf.shape[0]/2,:])
    plt.semilogy(x,psfi[psf.shape[0]/2,:],color="green")
    plt.semilogy(x,psf_compass[psf.shape[0]/2,:],color="red")
    if(correction):
        #plt.semilogy(x,psfc[psf.shape[0]/2,:],color="purple")
        if(synth):
            plt.semilogy(x,psfs[psf.shape[0]/2,:],color="black")
            #plt.legend(["PSF rec","PSF ind. assumption", "PSF COMPASS", "PSF corrected", "PSF synth"])
            plt.legend(["PSF rec","PSF ind. assumption", "PSF COMPASS", "PSF synth"])
        else:
            plt.legend(["PSF rec","PSF ind. assumption", "PSF COMPASS", "PSF corrected"])
    elif(synth):
        plt.semilogy(x,psfs[psf.shape[0]/2,:],color="purple")
        plt.legend(["PSF rec","PSF ind. assumption", "PSF COMPASS", "PSF synth"])

    else:
        plt.legend(["PSF rec","PSF ind. assumption", "PSF COMPASS"])
    plt.xlabel("Angular distance [units of lambda/D]")
    plt.ylabel("Normalized intensity")
    plt.title(filename.split('/')[-1])

    f.close()
    return psf_compass, psf, psfi, psfc, psfs


def filter_piston_TT(filename,C):
    IF, T = rexp.get_IF(filename)
    IF = IF.T
    T=T.T
    N = IF.shape[0]
    n = IF.shape[1]

    delta = IF.T.dot(IF).toarray()/N

    # Tip-tilt + piston
    Tp = np.ones((T.shape[0],T.shape[1]+1))
    Tp[:,:2] = T.copy()#.toarray()
    deltaT = IF.T.dot(Tp)/N
    # Tip tilt projection on the pzt dm
    tau = np.linalg.inv(delta).dot(deltaT)

    # Famille generatrice sans tip tilt
    G = np.identity(n)
    tdt = tau.T.dot(delta).dot(tau)
    subTT = tau.dot(np.linalg.inv(tdt)).dot(tau.T).dot(delta)
    G -= subTT

    return G.T.dot(C).dot(G)


def filter_TT(filename,C):
    IF, T = rexp.get_IF(filename)
    IF = IF.T
    T=T.T
    N = IF.shape[0]
    n = IF.shape[1]

    delta = IF.T.dot(IF).toarray()/N

    deltaT = IF.T.dot(T)/N
    # Tip tilt projection on the pzt dm
    tau = np.linalg.inv(delta).dot(deltaT)

    # Famille generatrice sans tip tilt
    G = np.identity(n)
    tdt = tau.T.dot(delta).dot(tau)
    subTT = tau.dot(np.linalg.inv(tdt)).dot(tau.T).dot(delta)
    G -= subTT

    return G.T.dot(C).dot(G)


def compute_covariance_model(filename):
    f = h5py.File(filename,'r')
    L0 = f.attrs["L0"][0]
    Lambda_tar = f.attrs["target.Lambda"][0]
    Lambda_wfs = f.attrs["wfs.Lambda"]
    speed = f.attrs["windspeed"][0]
    dt = f.attrs["ittime"]
    gain = f.attrs["gain"]
    H = f.attrs["atm.alt"][0]
    theta = f.attrs["winddir"][0]*np.pi/180.
    wxpos = f.attrs["wfs.xpos"][0]
    wypos = f.attrs["wfs.ypos"][0]
    r0 = f.attrs["r0"] * (Lambda_tar/Lambda_wfs)**(6./5.)
    RASC = 180./np.pi * 3600.

    Htheta = np.linalg.norm([wxpos,wypos])/RASC*H
    vdt = speed*dt/gain
    gamma = np.arctan2(wypos,wxpos) - theta
    rho =  np.sqrt(Htheta**2 + (vdt)**2 - 2*Htheta*vdt*np.cos(np.pi-gamma))
    # Covariance matrices models on actuators space
    xpos = f["dm.xpos"][:]
    ypos = f["dm.ypos"][:]
    p2m = f.attrs["tel_diam"] / f.attrs["pupdiam"]
    pupshape = long(2 ** np.ceil(np.log2(f.attrs["pupdiam"]) + 1))
    xactu = (xpos - pupshape/2) * p2m
    yactu = (ypos - pupshape/2) * p2m
    M = np.zeros((xpos.size,xpos.size))
    Mvdt = M.copy()
    Mht = M.copy()
    Mhvdt = M.copy()
    angleht = np.arctan2(wypos,wxpos)
    anglehvdt = gamma/2. - theta
    fc = xactu[1] - xactu[0]
    #fc = 0.05

    for i in range(xpos.size):
        for j in range(xpos.size):
            Mvdt[i,j] = (np.sqrt((xactu[i]-(xactu[j]-vdt*np.cos(theta)))**2 + (yactu[i]-(yactu[j]-vdt*np.sin(theta)))**2))
            M[i,j] = (np.sqrt((xactu[i]-xactu[j])**2 + (yactu[i]-yactu[j])**2))
            Mht[i,j] = (np.sqrt((xactu[i]-(xactu[j]-Htheta*np.cos(angleht)))**2 + (yactu[i]-(yactu[j]-Htheta*np.sin(angleht)))**2))
            #Mhvdt[i,j] = (np.sqrt((xactu[i]-(xactu[j]+rho*np.cos(anglehvdt)))**2 + (yactu[i]-(yactu[j]+rho*np.sin(anglehvdt)))**2))
            Mhvdt[i,j] = (np.sqrt(((xactu[i]+vdt*np.cos(theta))-(xactu[j]-Htheta*np.cos(angleht)))**2 + ((yactu[i]+vdt*np.sin(theta))-(yactu[j]-Htheta*np.sin(angleht)))**2))

    Ccov =  0.5 * (Dphi.dphi_lowpass(Mhvdt,fc,L0,tabx,taby) - Dphi.dphi_lowpass(Mht,fc,L0,tabx,taby) \
            - Dphi.dphi_lowpass(Mvdt,fc,L0,tabx,taby) + Dphi.dphi_lowpass(M,fc,L0,tabx,taby)) * (1./r0)**(5./3.)

    Caniso = 0.5*(Dphi.dphi_lowpass(Mht,fc,L0,tabx,taby) - Dphi.dphi_lowpass(M,fc,L0,tabx,taby)) * (1./r0)**(5./3.)
    Cbp = 0.5*(Dphi.dphi_lowpass(Mvdt,fc,L0,tabx,taby) - Dphi.dphi_lowpass(M,fc,L0,tabx,taby)) * (1./r0)**(5./3.)

    Sp = (f.attrs["tel_diam"]/f.attrs["nxsub"])**2/2.
    f.close()
    return (Caniso+Caniso.T)*Sp, (Cbp+Cbp.T)*Sp, Ccov*Sp


def add_TT_model(filename,Cmodel):
    C = np.zeros((Cmodel.shape[0]+2,Cmodel.shape[0]+2))
    IF, T = rexp.get_IF(filename)
    IF = IF.T
    T = T.T
    N = IF.shape[0]
    deltaTT = T.T.dot(T)/N
    deltaF = IF.T.dot(T)/N
    pzt2tt = np.linalg.inv(deltaTT).dot(deltaF.T)

    CTT = Cmodel - filter_TT(filename,Cmodel)
    C[-2:,-2:] = pzt2tt.dot(CTT).dot(pzt2tt.T)
    C[:-2,:-2] = Cmodel

    return C


def load_datas(files):
    nmodes = (files[0])["P"][:].shape[0]
    P = (files[0])["P"][:]
    xpos = files[0].attrs["wfs.xpos"][0]
    ypos = files[0].attrs["wfs.ypos"][0]
    contributors = ["tomography", "bandwidth"]
    Lambda_tar = files[0].attrs["target.Lambda"][0]
    Lambda_wfs = files[0].attrs["wfs.Lambda"][0]
    L0 = files[0].attrs["L0"][0]
    dt = files[0].attrs["ittime"]
    H = files[0].attrs["atm.alt"][0]
    RASC = 180/np.pi * 3600.
    Htheta = np.linalg.norm([xpos,ypos])/RASC*H# np.sqrt(2)*4/RASC*H # Hardcoded for angular separation of sqrt(2)*4 arcsec
    r0 = files[0].attrs["r0"] * (Lambda_tar/Lambda_wfs)**(6./5.)
    nfiles = len(files)
    vartomo = np.zeros((nfiles,nmodes))
    varbp = np.zeros((nfiles,nmodes))
    vartot = np.zeros((nfiles,nmodes))
    theta = np.zeros(nfiles)
    speeds = np.zeros(nfiles)
    gain = np.zeros(nfiles)
    # data[:,0,i] = var(tomo+bp) for file #i
    # data[:,1,i] = var(tomo) for file #i
    # data[:,2,i] = var(bp) for file #i
    # data[:,3,i] = var(tomo)+var(bp) for file #i
    ind = 0
    print "Loading data..."
    for f in files:
        vartot[ind,:] = rexp.variance(f, contributors) * ((2*np.pi/Lambda_tar)**2)
        vartomo[ind,:] = rexp.variance(f, ["tomography"]) * ((2*np.pi/Lambda_tar)**2)
        varbp[ind,:] = rexp.variance(f, ["bandwidth"]) * ((2*np.pi/Lambda_tar)**2)
        theta[ind] = f.attrs["winddir"][0]
        speeds[ind] = f.attrs["windspeed"][0]
        gain[ind] = float('%.1f' % f.attrs["gain"][0])
        ind += 1
        print ind,"/",len(files)

    covar = (vartot - (vartomo+varbp))/2.

    stot = np.sum(vartot,axis=1)
    sbp = np.sum(varbp,axis=1)
    stomo = np.sum(vartomo,axis=1)
    scov = np.sum(covar,axis=1)

    return stot, sbp, stomo, scov, covar


datapath = '/home/fferreira/Data/correlation/'
filenames = glob.glob(datapath + 'roket_8m_1layer_dir*_cpu.h5')

tabx, taby = Dphi.tabulateIj0()
# Illustration du probleme
#psf_compass, psf, psfi, psfc, psfs = compute_and_compare_PSFs(filenames[11],correction=True, synth=True)



'''
files = []
for f in filenames:
    files.append(h5py.File(f,'r'))


# Datas

nmodes = (files[0])["P"][:].shape[0]
P = (files[0])["P"][:]
xpos = files[0].attrs["wfs.xpos"][0]
ypos = files[0].attrs["wfs.ypos"][0]
contributors = ["tomography", "bandwidth"]
Lambda_tar = files[0].attrs["target.Lambda"][0]
Lambda_wfs = files[0].attrs["wfs.Lambda"][0]
L0 = files[0].attrs["L0"][0]
dt = files[0].attrs["ittime"]
H = files[0].attrs["atm.alt"][0]
RASC = 180/np.pi * 3600.
Htheta = np.linalg.norm([xpos,ypos])/RASC*H# np.sqrt(2)*4/RASC*H # Hardcoded for angular separation of sqrt(2)*4 arcsec
r0 = files[0].attrs["r0"] * (Lambda_tar/Lambda_wfs)**(6./5.)
nfiles = len(files)
vartomo = np.zeros((nfiles,nmodes))
varbp = np.zeros((nfiles,nmodes))
vartot = np.zeros((nfiles,nmodes))
theta = np.zeros(nfiles)
speeds = np.zeros(nfiles)
gain = np.zeros(nfiles)

# Illustration du probleme
otftel, otf2, psf, gpu = precs.psf_rec_Vii(filenames[11])
cov_err = rexp.get_coverr_independence(filenames[11])
otfteli, otf2i, psfi, gpu = precs.psf_rec_Vii(filenames[11],covmodes=cov_err)
psf_compass = np.fft.fftshift(files[11]["psf"][:])
RASC = 180/np.pi*3600.
pixsize = Lambda_tar*1e-6  / (psf.shape[0] * 8./640) * RASC
x = (np.arange(psf.shape[0]) - psf.shape[0]/2) * pixsize / (Lambda_tar*1e-6/8. * RASC)
plt.semilogy(x,psf[psf.shape[0]/2,:])
plt.semilogy(x,psfi[psf.shape[0]/2,:],color="green")
plt.semilogy(x,psf_compass[psf.shape[0]/2,:],color="red")
plt.xlabel("Angular distance [units of lambda/D]")
plt.ylabel("Normalized intensity")
plt.legend(["PSF rec","PSF ind. assumption", "PSF COMPASS"])


# data[:,0,i] = var(tomo+bp) for file #i
# data[:,1,i] = var(tomo) for file #i
# data[:,2,i] = var(bp) for file #i
# data[:,3,i] = var(tomo)+var(bp) for file #i
ind = 0
print "Loading data..."
for f in files:
    #vartot[ind,:] = rexp.variance(f, contributors) * ((2*np.pi/Lambda_tar)**2)
    #vartomo[ind,:] = rexp.variance(f, ["tomography"]) * ((2*np.pi/Lambda_tar)**2)
    #varbp[ind,:] = rexp.variance(f, ["bandwidth"]) * ((2*np.pi/Lambda_tar)**2)
    theta[ind] = f.attrs["winddir"][0]
    speeds[ind] = f.attrs["windspeed"][0]
    gain[ind] = float('%.1f' % f.attrs["gain"][0])
    ind += 1
    print ind,"/",len(files)

covar = (vartot - (vartomo+varbp))/2.

stot = np.sum(vartot,axis=1)
sbp = np.sum(varbp,axis=1)
stomo = np.sum(vartomo,axis=1)
scov = np.sum(covar,axis=1)

# Model

print "Building models..."
vdt = speeds*dt/gain
Htheta = np.ones(nfiles) * Htheta
gamma = np.arctan2(ypos,xpos) - theta*np.pi/180.
rho =  np.sqrt(Htheta**2 + (vdt)**2 - 2*Htheta*vdt*np.cos(np.pi-gamma))
# Covariance matrices models on actuators space
xpos = files[11]["dm.xpos"][:]
ypos = files[11]["dm.ypos"][:]
p2m = files[11].attrs["tel_diam"] / files[11].attrs["pupdiam"]
pupshape = long(2 ** np.ceil(np.log2(files[11].attrs["pupdiam"]) + 1))
xactu = (xpos - pupshape/2) * p2m
yactu = (ypos - pupshape/2) * p2m
M = np.zeros((1304,1304))
Mvdt = M.copy()
Mht = M.copy()
Mhvdt = M.copy()
angleht = np.arctan2(files[11].attrs["wfs.ypos"][0],files[11].attrs["wfs.xpos"][0])
anglehvdt = gamma/2. - theta*np.pi/180.
thetar = theta*np.pi/180.

for i in range(1304):
    for j in range(1304):
        Mvdt[i,j] = (np.sqrt((xactu[i]-(xactu[j]+vdt[11]*np.cos(thetar[11])))**2 + (yactu[i]-(yactu[j]+vdt[11]*np.sin(thetar[11])))**2))
        M[i,j] = (np.sqrt((xactu[i]-xactu[j])**2 + (yactu[i]-yactu[j])**2))
        Mht[i,j] = (np.sqrt((xactu[i]-(xactu[j]+Htheta[11]*np.cos(angleht)))**2 + (yactu[i]-(yactu[j]+Htheta[11]*np.sin(angleht)))**2))
        #Mhvdt[i,j] = (np.sqrt((xactu[i]-(xactu[j]+rho[11]*np.cos(anglehvdt[11])))**2 + (yactu[i]-(yactu[j]+rho[11]*np.sin(anglehvdt[11])))**2))
        Mhvdt[i,j] = (np.sqrt(((xactu[i]-vdt[11]*np.cos(thetar[11]))-(xactu[j]+Htheta[11]*np.cos(angleht)))**2 + ((yactu[i]-vdt[11]*np.sin(thetar[11]))-(yactu[j]+Htheta[11]*np.sin(angleht)))**2))

Cmodel =  (Dphi.dphi_lowpass(Mhvdt,0.2,L0,tabx,taby) - Dphi.dphi_lowpass(Mht,0.2,L0,tabx,taby) \
        - Dphi.dphi_lowpass(Mvdt,0.2,L0,tabx,taby) + Dphi.dphi_lowpass(M,0.2,L0,tabx,taby)) * (1/r0)**(5./3.)
'''
'''
mtomo = Dphi.dphi_lowpass(Htheta,0.2,L0, tabx, taby) * (1/r0)**(5./3.)
mbp = Dphi.dphi_lowpass(vdt ,0.2, L0, tabx, taby) * (1/r0)**(5./3.)
mtot = Dphi.dphi_lowpass(rho,0.2,L0,tabx,taby) * (1/r0)**(5./3.)

# Piston correction
print "Computing piston correction..."
pup = precs.get_pup(filenames[11])
r = np.zeros((8192,8192))
p2m = files[11].attrs["tel_diam"]/pup.shape[0]
Npts = files[11]["indx_pup"].size
for k in range(r.shape[0]):
    for j in range(r.shape[0]):
        r[k,j] = np.sqrt((k-r.shape[0]/2+0.5)**2+(j-r.shape[0]/2+0.5)**2) * p2m

ipup = np.zeros((8192,8192))
ipup[3776:3776+640,3776:3776+640] = pup
dphi_map = Dphi.dphi_lowpass(r,0.2,L0,tabx,taby) * (1/r0)**(5./3.)
fpup = np.fft.fft2(ipup)#,s=[8192,8192])
fdphi = np.fft.fft2(dphi_map)#,s=[8192,8192])
fconv = fpup * fpup * fdphi
dconv = np.fft.ifft2(fconv).real / Npts / Npts
mini = np.where(dconv == dconv.min())
dutil = dconv[mini[0][0],mini[1][0]:]
#Avdt = dconv[mini[0]+(vdt/p2m).astype(np.int16),mini[1]] - dconv[mini]
Pbp = np.interp(vdt/p2m,np.arange(dutil.shape[0]),dutil) - dutil[0]
Ptomo = (np.interp(Htheta/gain/p2m,np.arange(dutil.shape[0]),dutil) - dutil[0])*gain**2
Ptot = np.interp(rho/p2m,np.arange(dutil.shape[0]),dutil) - dutil[0]

mtomo -= Ptomo
mbp -= Pbp
mtot -= Ptot
mcov = (-mtomo - mbp + mtot)*0.5


# Correction on psf
m = (np.arange(nmodes)+1)**(-5/6.)
m /= m.sum()
m = m * (mcov[11] / (2*np.pi/Lambda_tar)**2)
'''
'''
cov_err2 = P.dot(cov_err).dot(P.T) + 2*np.diag(m)
otftelc, otf2c, psfc, gpu = precs.psf_rec_Vii(filenames[11],cov=cov_err2.astype(np.float32))
plt.figure()
plt.semilogy(x,psf_compass[psf.shape[0]/2,:],color="red")
plt.semilogy(x,psfi[psf.shape[0]/2,:],color="green")
plt.semilogy(x,psfc[psf.shape[0]/2,:],color="blue")
plt.xlabel("Angular distance [units of lambda/D]")
plt.ylabel("Normalized intensity")
plt.legend([ "PSF COMPASS","PSF ind. assumption","PSF corrected"])
'''


'''
xpos = files[11]["dm.xpos"][:]
ypos = files[11]["dm.ypos"][:]
dm_dim = files[11]["dm_dim"].value
xpos -= (pupshape-dm_dim)/2
ypos -= (pupshape-dm_dim)/2
influ = np.load("influ.npy")
influ2 = np.zeros((dm_dim,dm_dim))
tmp = influ2.copy()
indx_pup = files[11]["indx_pup"][:]
pup = np.zeros((dm_dim,dm_dim)).flatten()
pup[indx_pup] = 1
pup = pup.reshape((dm_dim,dm_dim))
ind2 = np.where(pup)
influshape = influ.shape[0]
A = np.zeros((xpos.size,xpos.size))
rr = Htheta[11]
xx = np.cos(0*theta[11]*np.pi/180.)*rr/p2m
yy = np.sin(0*theta[11]*np.pi/180.)*rr/p2m
pup2 = pup.copy()*0.
pup2[(ind2[0]+xx).astype(np.int32),(ind2[1]+yy).astype(np.int32)] = 1.

for i in range(xpos.size):
    influ2 *=0
    influ2[xpos[i]-influshape/2+1:xpos[i]+influshape/2+1,ypos[i]-influshape/2+1:ypos[i]+influshape/2+1] = influ
    influ2 *= pup
    for j in range(xpos.size):
        if(tmp[xpos[j]-influshape/2+1+xx:xpos[j]+influshape/2+1+xx,ypos[j]-influshape/2+1+yy:ypos[j]+influshape/2+1+yy].shape == influ.shape):
            tmp *=0
            tmp[xpos[j]-influshape/2+1+xx:xpos[j]+influshape/2+1+xx,ypos[j]-influshape/2+1+yy:ypos[j]+influshape/2+1+yy] = influ
            tmp *= pup2
            A[i,j] = (influ2*tmp).sum()
        else:
            A[i,j] = 0.
    print i
'''
