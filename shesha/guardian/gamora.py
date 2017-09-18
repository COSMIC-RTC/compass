"""
Created on Thu Jul 21 09:28:23 2016

@author: fferreira
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py
import naga as ch
from Gamora import gamora_init, Gamora
from scipy.sparse import csr_matrix
from sys import stdout
import time

plt.ion()
gpudevices = np.array([6, 7], dtype=np.int32)
c = ch.naga_context(devices=gpudevices)

#filename = "/home/fferreira/Data/breakdown_offaxis-4_2.h5"


def cutsPSF(filename, psf, psfs):
    f = h5py.File(filename, 'r')
    Lambda_tar = f.attrs["target.Lambda"][0]
    RASC = 180 / np.pi * 3600.
    pixsize = Lambda_tar * 1e-6 / \
        (psf.shape[0] * f.attrs["tel_diam"] / f.attrs["pupdiam"]) * RASC
    x = (np.arange(psf.shape[0]) - psf.shape[0] / 2) * \
        pixsize / (Lambda_tar * 1e-6 / f.attrs["tel_diam"] * RASC)
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.semilogy(x, psf[psf.shape[0] / 2, :], color="blue")
    plt.semilogy(x, psfs[psf.shape[0] / 2, :], color="red")
    plt.xlabel("X-axis angular distance [units of lambda/D]")
    plt.ylabel("Normalized intensity")
    plt.legend(["PSF 1", "PSF 2"])
    plt.xlim(-20, 20)
    plt.ylim(1e-5, 1)
    plt.subplot(2, 1, 2)
    plt.semilogy(x, psf[:, psf.shape[0] / 2], color="blue")
    plt.semilogy(x, psfs[:, psf.shape[0] / 2], color="red")
    plt.xlabel("Y-axis angular distance [units of lambda/D]")
    plt.ylabel("Normalized intensity")
    plt.legend(["PSF 1", "PSF 2"])
    plt.xlim(-20, 20)
    plt.ylim(1e-5, 1)
    f.close()


def get_err(filename):
    f = h5py.File(filename, 'r')
    # Get the sum of error contributors
    err = f["noise"][:]
    err += f["aliasing"][:]
    err += f["tomography"][:]
    err += f["filtered modes"][:]
    err += f["non linearity"][:]
    err += f["bandwidth"][:]
    f.close()

    return err


def get_err_contributors(filename, contributors):
    f = h5py.File(filename, 'r')
    # Get the sum of error contributors
    err = f["noise"][:] * 0.
    for c in contributors:
        err += f[c][:]
    f.close()

    return err


def get_pup(filename):
    f = h5py.File(filename, 'r')
    if (list(f.keys()).count("spup")):
        spup = f["spup"][:]
    else:
        indx_pup = f["indx_pup"][:]
        pup = np.zeros((f["dm_dim"].value, f["dm_dim"].value))
        pup_F = pup.flatten()
        pup_F[indx_pup] = 1.
        pup = pup_F.reshape(pup.shape)
        spup = pup[np.where(pup)[0].min():np.where(pup)[0].max() + 1,
                   np.where(pup)[1].min():np.where(pup)[1].max() + 1]

    f.close()
    return spup


def get_IF(filename):
    f = h5py.File(filename, 'r')
    IF = csr_matrix((f["IF.data"][:], f["IF.indices"][:], f["IF.indptr"][:]))
    if (list(f.keys()).count("TT")):
        T = f["TT"][:]
    else:
        T = IF[-2:, :].toarray()
        IF = IF[:-2, :]
    f.close()
    return IF, T.T.astype(np.float32)


def psf_rec_roket_file(filename, err=None):
    f = h5py.File(filename, 'r')
    if (err is None):
        err = get_err(filename)
    spup = get_pup(filename)
    # Sparse IF matrix
    IF, T = get_IF(filename)
    # Scale factor
    scale = float(2 * np.pi / f.attrs["target.Lambda"][0])
    # Init GPU
    gamora = gamora_init(b"roket", err.shape[0], err.shape[1],
                         IF.data.astype(np.float32), IF.indices, IF.indptr, T,
                         spup.astype(np.float32), scale)
    # Launch computation
    gamora.psf_rec_roket(err)
    # Get psf
    psf = gamora.get_psf()
    f.close()
    return psf, gamora


def psf_rec_roket_file_cpu(filename):
    f = h5py.File(filename, 'r')
    # Get the sum of error contributors
    err = get_err(filename)

    # Retrieving spupil (for file where spupil was not saved)
    indx_pup = f["indx_pup"][:]
    pup = np.zeros((f["dm_dim"].value, f["dm_dim"].value))
    pup_F = pup.flatten()
    pup_F[indx_pup] = 1.
    pup = pup_F.reshape(pup.shape)
    spup = pup[np.where(pup)[0].min():np.where(pup)[0].max() + 1,
               np.where(pup)[1].min():np.where(pup)[1].max() + 1]
    phase = spup.copy()
    mradix = 2
    fft_size = mradix**int((np.log(2 * spup.shape[0]) / np.log(mradix)) + 1)
    amplipup = np.zeros((fft_size, fft_size), dtype=np.complex)
    psf = amplipup.copy()
    psf = psf.astype(np.float32)

    # Sparse IF matrix
    IF, T = get_IF(filename)
    # Scale factor
    scale = float(2 * np.pi / f.attrs["target.Lambda"][0])

    for k in range(err.shape[1]):
        amplipup = np.zeros((fft_size, fft_size), dtype=np.complex)
        phase[np.where(spup)] = IF.T.dot(err[:-2, k])
        phase[np.where(spup)] += T.dot(err[-2:, k])
        amplipup[:phase.shape[0], :phase.shape[1]] = np.exp(-1j * phase * scale)
        amplipup = np.fft.fft2(amplipup)
        psf += np.fft.fftshift(np.abs(amplipup)**2) / \
            IF.shape[1] / IF.shape[1] / err.shape[1]
        print(" Computing and stacking PSF : %d/%d\r" % (i, n), end=' ')
    print("PSF computed and stacked")
    f.close()
    return psf


def psf_rec_Vii(filename, err=None, fitting=True, covmodes=None, cov=None,
                onlypsf=False):
    f = h5py.File(filename, 'r')
    spup = get_pup(filename)
    # Sparse IF matrix
    IF, T = get_IF(filename)
    # Covariance matrix
    P = f["P"][:]
    print("Projecting error buffer into modal space...")
    if ((err is None) and (cov is None)):
        err = get_err(filename)
        err = P.dot(err)
    print("Computing covariance matrix...")
    if (cov is None):
        if (covmodes is None):
            covmodes = err.dot(err.T) / err.shape[1]
        else:
            covmodes = (P.dot(covmodes)).dot(P.T)
    else:
        covmodes = cov
    print("Done")
    Btt = f["Btt"][:]

    # Scale factor
    scale = float(2 * np.pi / f.attrs["target.Lambda"][0])
    # Init GPU
    gpu = gamora_init(b"Vii", Btt.shape[0], f["noise"][:].shape[1],
                      IF.data.astype(np.float32), IF.indices, IF.indptr, T,
                      spup.astype(np.float32), scale, covmodes.shape[0], Btt, covmodes)
    # Launch computation
    # gamora.set_eigenvals(e.astype(np.float32))
    # gamora.set_covmodes(V.astype(np.float32))
    tic = time.time()
    gpu.psf_rec_Vii()

    otftel = gpu.get_otftel()
    otf2 = gpu.get_otfVii()

    otftel /= otftel.max()
    if (list(f.keys()).count("psfortho") and fitting):
        print("\nAdding fitting to PSF...")
        psfortho = f["psfortho"][:]
        otffit = np.real(np.fft.fft2(psfortho))
        otffit /= otffit.max()
        psf = np.fft.fftshift(np.real(np.fft.ifft2(otffit * otf2)))
    else:
        psf = np.fft.fftshift(np.real(np.fft.ifft2(otftel * otf2)))

    psf *= (psf.shape[0] * psf.shape[0] / float(np.where(spup)[0].shape[0]))
    f.close()
    tac = time.time()
    print(" ")
    print("PSF renconstruction took ", tac - tic, " seconds")
    if onlypsf:
        return psf
    else:
        return otftel, otf2, psf, gpu


def psf_rec_vii_cpu(filename):
    f = h5py.File(filename, 'r')
    IF, T = get_IF(filename)
    ratio_lambda = 2 * np.pi / f.attrs["target.Lambda"][0]
    # Telescope OTF
    print("Computing telescope OTF...")
    spup = get_pup(filename)
    mradix = 2
    fft_size = mradix**int((np.log(2 * spup.shape[0]) / np.log(mradix)) + 1)
    pup = np.zeros((fft_size, fft_size))
    pup[:spup.shape[0], :spup.shape[0]] = spup
    pupfft = np.fft.fft2(pup)
    conjpupfft = np.conjugate(pupfft)
    otftel = np.real(np.fft.ifft2(pupfft * conjpupfft))
    den = 1. / otftel
    den[np.where(np.isinf(den))] = 0
    mask = np.ones((fft_size, fft_size))
    mask[np.where(otftel < 1e-5)] = 0
    otftel = otftel / otftel.max()
    print("Done")
    # Covariance matrix
    print("Computing covariance matrix...")
    err = get_err(filename)
    P = f["P"][:]
    err = P.dot(err)
    Btt = f["Btt"][:]
    #modes = IF.T.dot(Btt)
    covmodes = err.dot(err.T) / err.shape[1]
    print("Done")
    # Vii algorithm
    print("Diagonalizing cov matrix...")
    e, V = np.linalg.eig(covmodes)
    print("Done")
    tmp = np.zeros((fft_size, fft_size))
    newmodek = tmp.copy()
    ind = np.where(pup)
    for k in range(err.shape[0]):
        #newmodek[ind] = IF.T.dot(V[:,k])
        #newmodek[ind] = modes.dot(V[:,k])
        tmp2 = Btt.dot(V[:, k])
        newmodek[ind] = IF.T.dot(tmp2[:-2])
        newmodek[ind] += T.T.dot(tmp2[-2:])
        term1 = np.real(np.fft.fft2(newmodek**2) * conjpupfft)
        term2 = np.abs(np.fft.fft2(newmodek))**2
        tmp += ((term1 - term2) * e[k])
        print(" Computing Vii : %d/%d\r" % (k, covmodes.shape[0]), end=' ')
    print("Vii computed")

    dphi = np.real(np.fft.ifft2(2 * tmp)) * den * mask * ratio_lambda**2
    otf2 = np.exp(-0.5 * dphi) * mask
    otf2 = otf2 / otf2.max()

    psf = np.fft.fftshift(np.real(np.fft.ifft2(otftel * otf2)))
    psf *= (fft_size * fft_size / float(np.where(pup)[0].shape[0]))

    f.close()
    return otftel, otf2, psf


def test_Vii(filename):
    a = time.time()
    otftel_cpu, otf2_cpu, psf_cpu = psf_rec_vii_cpu(filename)
    b = time.time()
    otftel_gpu, otf2_gpu, psf_gpu, gamora = psf_rec_Vii(filename)
    c = time.time()
    cputime = b - a
    gputime = c - b
    print("CPU exec time : ", cputime, " s")
    print("GPU exec time : ", gputime, " s")
    print("Speed up : x", cputime / gputime)
    print("---------------------------------")
    print("precision on psf : ", np.abs(psf_cpu - psf_gpu).max() / psf_cpu.max())


def add_fitting_to_psf(otf2, psf_fit):
    print("\nAdding fitting to PSF...")
    otffit = np.real(np.fft.fft2(psf_fit))
    otffit /= otffit.max()
    psf = np.fft.fftshift(np.real(np.fft.ifft2(otffit * otf2)))
    psf *= (psf.shape[0] * psf.shape[0] / float(np.where(spup)[0].shape[0]))

    return psf
