## @package   shesha.ao.cmats
## @brief     Computation implementations of command matrix
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @version   5.5.0
## @date      2022/01/24
## @copyright 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>
#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>

# COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or any later version.

# COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License along with COMPASS.
# If not, see <https://www.gnu.org/licenses/>

# Copyright (C) 2011-2024 COSMIC Team <https//://github.com/COSMIC-RTC/compass>


import numpy as np
import time

import shesha.config as conf
import shesha.constants as scons

from shesha.sutra_wrap import Rtc_FFF as Rtc

from shesha.ao.wfs import noise_cov

from typing import List


def generic_imat_inversion(
    M2V: np.ndarray,
    modalIMat: np.ndarray,
    modeSelect: np.ndarray = None,
    modeGains: np.ndarray = None,
) -> np.ndarray:
    """Generic numpy modal interaction matrix inversion function

    Args:
        M2V: (nActu x nModes) : modal basis matrix
        modalIMat: (nSlopes x nModes) : modal interaction matrix
        modeSelect: (nModes, dtype=bool): (Optional): mode selection, mode at False is filtered
        modeGains: (nModes, dtype=bool): (Optional):
        modal gains to apply. These are gain in the reconstruction sens, ie
        they are applied multiplicatively on the command matrix
    """
    if modeSelect is None:
        modeSelect = np.ones(modalIMat.shape[1], dtype=bool)
    if modeGains is None:
        modeGains = np.ones(modalIMat.shape[1], dtype=np.float32)

    return M2V.dot(
        modeGains[:, None]
        * np.linalg.inv(modalIMat[:, modeSelect].T.dot(modalIMat[:, modeSelect])).dot(
            modalIMat[:, modeSelect].T
        )
    )


def cmat_init(
    ncontrol: int,
    rtc: Rtc,
    p_controller: conf.ParamController,
    p_wfss: List[conf.ParamWfs],
    p_atmos: conf.ParamAtmos,
    p_tel: conf.ParamTel,
    p_dms: List[conf.ParamDm],
    nmodes: int = 0,
) -> None:
    """Compute the command matrix on the GPU

    Args:
        ncontrol: (int) : controller index
        rtc: (Rtc) : rtc object
        p_controller: (ParamController) : controller settings
        p_wfss: (list of ParamWfs) : wfs settings
        p_atmos: (ParamAtmos) : atmos settings
        p_tel : (ParamTel) : telescope settings
        p_dms: (list of ParamDm) : dms settings
        M2V : (np.ndarray[ndim=2, dtype=np.float32]): (optional) KL to volts matrix (for KL cmat)
        nmodes: (int) : (optional) number of kl modes
    """
    if p_controller.type == scons.ControllerType.LS:
        cmat_ls_init(ncontrol, rtc, p_controller)
    elif p_controller.type == scons.ControllerType.MV:
        cmat_mv_init(ncontrol, rtc, p_controller, p_wfss, p_atmos, p_tel, p_dms)
    p_controller.set_cmat(np.array(rtc.d_control[ncontrol].d_cmat))


def cmat_mv_init(ncontrol, rtc, p_controller, p_wfss, p_atmos, p_tel, p_dms):
    """Compute the command matrix for a minimum variance controller

    Args:
            ncontrol: (int) :
            rtc: (Rtc) :
            p_controller: (ParamController) : controller settings
            p_wfss: (list of ParamWfs) : wfs settings
            p_atmos: (ParamAtmos) : atmos settings
            p_tel : (ParamTel) : telescope settings
            p_dms: (list of ParamDm) : dms settings
    """
    Cn = np.zeros(p_controller._imat.shape[0], dtype=np.float32)
    ind = 0
    for k in p_controller.nwfs:
        Cn[ind : ind + 2 * p_wfss[k]._nvalid] = noise_cov(k, p_wfss[k], p_atmos, p_tel)
        ind += 2 * p_wfss[k]._nvalid

    rtc.d_control[ncontrol].load_noisemat(Cn)
    print("Building cmat...")
    rtc.d_control[ncontrol].build_cmat(p_controller.maxcond)

    if p_controller.TTcond is None:
        p_controller.set_TTcond(p_controller.maxcond)

    if "tt" in [dm.type for dm in p_dms]:
        rtc.d_control[ncontrol].filter_cmat(p_controller.TTcond)
    print("Done")


def cmat_ls_init(ncontrol, rtc, p_controller):
    """Compute the command matrix for a least square controller

    Args:
                ncontrol: (int) : controller index
                rtc: (Rtc) : rtc object
                p_controller: (ParamController) : controller settings
    """
    print("Doing imat svd...")
    t0 = time.time()
    rtc.d_control[ncontrol].svdec_imat()
    print("svd done in %f s" % (time.time() - t0))
    eigenv = np.array(rtc.d_control[ncontrol].d_eigenvals)
    imat = np.array(rtc.d_control[ncontrol].d_imat)
    maxcond = p_controller.maxcond
    if eigenv[0] < eigenv[eigenv.shape[0] - 1]:
        mfilt = np.where((eigenv / eigenv[eigenv.shape[0] - 3]) < 1.0 / maxcond)[0]
    else:
        mfilt = np.where((1.0 / (eigenv / eigenv[2])) > maxcond)[0]
    nfilt = mfilt.shape[0]

    print("Building cmat...")
    t0 = time.time()
    if not p_controller.do_kl_imat:
        print("Filtering ", nfilt, " modes")
        rtc.d_control[ncontrol].build_cmat(nfilt)
    else:
        # filter imat
        D_filt = imat.copy()
        # Direct inversion
        Dp_filt = np.linalg.inv(D_filt.T.dot(D_filt)).dot(D_filt.T)
        if p_controller.klgain is not None:
            Dp_filt *= p_controller.klgain[None, :]
        cmat_filt = p_controller._M2V.dot(Dp_filt)
        rtc.d_control[ncontrol].set_cmat(cmat_filt)

    print("cmat done in %f s" % (time.time() - t0))


def svd_for_cmat(D):
    DtD = D.T.dot(D)
    return np.linalg.svd(DtD)


def Btt_for_cmat(rtc, dms, p_dms, p_geom):
    """Compute a command matrix in Btt modal basis (see error breakdown) and set
    it on the sutra_rtc. It computes by itself the volts to Btt matrix.

    Args:

        rtc: (Rtc) : rtc object
        dms: (Dms): dms object
        p_dms: (list of ParamDm): dms settings
        p_geom: (ParamGeom): geometry settings
    """
    from shesha.ao import basis

    IFs = basis.compute_IFsparse(dms, p_dms, p_geom).T
    n = IFs.shape[1]
    IFtt = IFs[:, -2:].toarray()
    IFpzt = IFs[:, : n - 2]

    Btt, P = basis.compute_btt(IFpzt, IFtt)
    return Btt, P


def get_cmat(D, nfilt, Btt=None, rtc=None, svd=None):
    """Compute a command matrix from an interaction matrix 'D'

    usage:
        get_cmat(D,nfilt)
        get_cmat(D,nfilt,Btt=BTT,rtc=RTC)
        get_cmat(D,nfilt,svd=SVD)

    Args:
        D: (np.ndarray[ndim=2, dtype=np.float32]): interaction matrix
        nfilt: (int): number of element to filter
        Btt: (np.ndarray[ndim=2, dtype=np.float32]): Btt modal basis
        rtc: (Rtc) : rtc object
        svd: (tuple of np.ndarray[ndim=1, dtype=np.float32): svd of D.T*D (obtained from np.linalg.svd)
    """
    nfilt = max(nfilt, 0)  # nfilt is positive
    if Btt is not None:
        if svd is not None:
            raise ValueError("Btt and SVD cannt be used together")
        if rtc is None:
            raise ValueError("Btt cannot be used without rtc")
        n = Btt.shape[1]
        index = np.concatenate((np.arange(n - nfilt - 2), np.array([n - 2, n - 1]))).astype(int)
        Btt_filt = Btt[:, index]
        # Modal interaction basis
        Dm = D.dot(Btt_filt)
        # Direct inversion
        Dmp = np.linalg.inv(Dm.T.dot(Dm)).dot(Dm.T)
        # Command matrix
        cmat = Btt_filt.dot(Dmp)
    else:
        if svd is not None:
            u = svd[0]
            s = svd[1]
            v = svd[2]
        else:
            u, s, v = svd_for_cmat(D)
        s_filt = 1 / s
        if nfilt > 0:
            s_filt[-nfilt:] = 0
        DtDx = v.T.dot(np.diag(s_filt)).dot(u.T)
        cmat = DtDx.dot(D.T)

    return cmat.astype(np.float32)
