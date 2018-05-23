"""
GAMORA (Gpu Accelerated Module fOr psf Reconstruction Algorithms)

Python module for GPU accelerated PSF reconstruction using Vii functions and ROKET file

Note: GPU devices used are hardcoded here. Change gpudevices if needed
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py
from naga.context import context
from shesha.sutra_bind.Gamora import gamora_init, Gamora
from scipy.sparse import csr_matrix
from sys import stdout
import time
from guardian.tools import roket_exploitation as rexp

plt.ion()

gpudevices = np.array([0], dtype=np.int32)
c = context(devices=gpudevices)


def psf_rec_roket_file(filename, err=None):
    """
    PSF reconstruction using ROKET file. SE PSF is reconstructed
    for each frame and stacked to obtain the LE PSF.
    Used for ROKET debug only.

    :parameters:
        filename: (str): path to the ROKET file
        err: (np.ndarray[ndim=2, dtype=np.float32]): (optionnal) Buffers of command error
    :return:
        psf: (np.ndarray[ndim=2, dtype=np.float32]): LE PSF
        gpu: (Gamora): Gamora GPU object (for manipulation and debug)
    """
    f = h5py.File(filename, 'r')
    if (err is None):
        err = rexp.get_err(filename)
    spup = rexp.get_pup(filename)
    # Sparse IF matrix
    IF, T = rexp.get_IF(filename)
    # Scale factor
    scale = float(2 * np.pi / f.attrs["_Param_target__Lambda"][0])
    # Init GPU
    gpu = gamora_init(b"roket", err.shape[0], err.shape[1],
                      IF.data.astype(np.float32), IF.indices, IF.indptr, T,
                      spup.astype(np.float32), scale)
    # Launch computation
    gpu.psf_rec_roket(err)
    # Get psf
    psf = gpu.get_psf()
    f.close()
    return psf, gpu


def psf_rec_roket_file_cpu(filename):
    """
    PSF reconstruction using ROKET file (CPU version). SE PSF is reconstructed
    for each frame and stacked to obtain the LE PSF.
    Used for ROKET debug only.

    :parameters:
        filename: (str): path to the ROKET file
    :return:
        psf: (np.ndarray[ndim=2, dtype=np.float32]): LE PSF
    """

    f = h5py.File(filename, 'r')
    # Get the sum of error contributors
    err = rexp.get_err(filename)

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
    IF, T = rexp.get_IF(filename)
    # Scale factor
    scale = float(2 * np.pi / f.attrs["_Param_target__Lambda"][0])

    for k in range(err.shape[1]):
        amplipup = np.zeros((fft_size, fft_size), dtype=np.complex)
        phase[np.where(spup)] = IF.T.dot(err[:-2, k])
        phase[np.where(spup)] += T.dot(err[-2:, k])
        amplipup[:phase.shape[0], :phase.shape[1]] = np.exp(-1j * phase * scale)
        amplipup = np.fft.fft2(amplipup)
        psf += np.fft.fftshift(np.abs(amplipup)**2) / \
            IF.shape[1] / IF.shape[1] / err.shape[1]
        print(" Computing and stacking PSF : %d/%d\r" % (k, err.shape[1]), end=' ')
    print("PSF computed and stacked")
    f.close()
    return psf


def psf_rec_Vii(filename, err=None, fitting=True, covmodes=None, cov=None):
    """
    PSF reconstruction using Vii functions with GPU acceleration.

    :parameters:
        filename: (str): path to the ROKET file
        err: (np.ndarray[ndim=2, dtype=np.float32]): (optionnal) Buffers of command error
        fitting: (bool): (optional) Add the fitting error to the PSF or not (True by default)
        covmodes: (np.ndarray[ndim=2, dtype=np.float32]): (optionnal) Error covariance matrix in the modal space
        cov: (np.ndarray[ndim=2, dtype=np.float32]): (optionnal) Error covariance matrix in the DM space

    :return:
        otftel: (np.ndarray[ndim=2, dtype=np.float32]): OTF of the perfect telescope
        otf2: (np.ndarray[ndim=2, dtype=np.float32]): OTF due to residual phase error
        psf: (np.ndarray[ndim=2, dtype=np.float32]): LE PSF
        gpu: (Gamora): Gamora GPU object (for manipulation and debug)
    """

    f = h5py.File(filename, 'r')
    spup = rexp.get_pup(filename)
    # Sparse IF matrix
    IF, T = rexp.get_IF(filename)
    # Covariance matrix
    P = f["P"][:]
    print("Projecting error buffer into modal space...")
    if ((err is None) and (cov is None)):
        err = rexp.get_err(filename)
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
    scale = float(2 * np.pi / f.attrs["_Param_target__Lambda"][0])
    # Init GPU
    gpu = gamora_init(b"Vii", Btt.shape[0], f["noise"][:].shape[1],
                      IF.data.astype(np.float32), IF.indices, IF.indptr, T,
                      spup.astype(np.float32), scale, covmodes.shape[0], Btt,
                      covmodes.astype(np.float32))
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

    return otftel, otf2, psf, gpu


def psf_rec_vii_cpu(filename):
    """
    PSF reconstruction using Vii functions (CPU version)

    :parameters:
        filename: (str): path to the ROKET file

    :return:
        otftel: (np.ndarray[ndim=2, dtype=np.float32]): OTF of the perfect telescope
        otf2: (np.ndarray[ndim=2, dtype=np.float32]): OTF due to residual phase error
        psf: (np.ndarray[ndim=2, dtype=np.float32]): LE PSF
    """

    f = h5py.File(filename, 'r')
    IF, T = rexp.get_IF(filename)
    ratio_lambda = 2 * np.pi / f.attrs["_Param_target__Lambda"][0]
    # Telescope OTF
    print("Computing telescope OTF...")
    spup = rexp.get_pup(filename)
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
    err = rexp.get_err(filename)
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
    """
    Test function comparing results and performance of GPU version
    versus CPU version of Vii PSF reconstruction

    :parameters:
        filename: (str): path to the ROKET file
    """
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


def add_fitting_to_psf(filename, otf, otffit):
    """
    Compute the PSF including the fitting OTF

    :parameters:
        otf: (np.ndarray[ndim=2, dtype=np.float32]): OTF
        otffit: (np.ndarray[ndim=2, dtype=np.float32]): Fitting error OTF
    :return:
        psf: (np.ndarray[ndim=2, dtype=np.float32]): PSF

    """
    print("\nAdding fitting to PSF...")
    spup = rexp.get_pup(filename)
    psf = np.fft.fftshift(np.real(np.fft.ifft2(otffit * otf)))
    psf *= (psf.shape[0] * psf.shape[0] / float(np.where(spup)[0].shape[0]))

    return psf
