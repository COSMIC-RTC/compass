"""
Created on Wed Oct 5 14:28:23 2016

@author: fferreira
"""

import numpy as np
import h5py
import glob
import sys
sys.path.append('/home/fferreira/compass/trunk/shesha/test/roket/tools/')
sys.path.append('/home/fferreira/compass/trunk/shesha/test/psf_reconstruction/')
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

datapath = '/home/fferreira/Data/correlation/'
filenames = glob.glob(datapath + 'roket_8m_1layer_dir*_cpu.h5')
files = []
for f in filenames:
    files.append(h5py.File(f,'r'))

Lambda_tar = files[0].attrs["target.Lambda"][0]

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

tabx, taby = Dphi.tabulateIj0()

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

# Model
print "Building models..."
vdt = speeds*dt/gain
Htheta = np.ones(nfiles) * Htheta

mtomo = Dphi.dphi_lowpass(Htheta,0.2,L0, tabx, taby) * (1/r0)**(5./3.)
mbp = Dphi.dphi_lowpass(vdt ,0.2, L0, tabx, taby) * (1/r0)**(5./3.)
gamma = np.arctan2(ypos,xpos) - theta*np.pi/180.
rho =  np.sqrt(Htheta**2 + (vdt)**2 - 2*Htheta*vdt*np.cos(np.pi-gamma))
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

cov_err2 = P.dot(cov_err).dot(P.T) + 2*np.diag(m)
otftelc, otf2c, psfc, gpu = precs.psf_rec_Vii(filenames[11],cov=cov_err2.astype(np.float32))
plt.figure()
plt.semilogy(x,psf_compass[psf.shape[0]/2,:],color="red")
plt.semilogy(x,psfi[psf.shape[0]/2,:],color="green")
plt.semilogy(x,psfc[psf.shape[0]/2,:],color="blue")
plt.xlabel("Angular distance [units of lambda/D]")
plt.ylabel("Normalized intensity")
plt.legend([ "PSF COMPASS","PSF ind. assumption","PSF corrected"])
