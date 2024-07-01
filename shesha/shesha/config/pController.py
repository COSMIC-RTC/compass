## @package   shesha.config.pController
## @brief     Class for the controller parameters
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   5.5.0
## @date      2022/01/24
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
#  All rights reserved.
#  Distributed under GNU - LGPL
#
#  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
#  General Public License as published by the Free Software Foundation, either version 3 of the License,
#  or any later version.
#
#  COMPASS: End-to-end AO simulation tool using GPU acceleration
#  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
#
#  The final product includes a software package for simulating all the critical subcomponents of AO,
#  particularly in the context of the ELT and a real-time core based on several control approaches,
#  with performances consistent with its integration into an instrument. Taking advantage of the specific
#  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
#  conduct large simulation campaigns called to the ELT.
#
#  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
#  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
#  various systems configurations such as multi-conjugate AO.
#
#  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
#  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
#  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.

import numpy as np
import shesha.config.config_setter_utils as csu
import shesha.constants as scons


class ParamController:
    """Class for the controller parameters

    Attributes:
        type (str): type of controller
        command_law (str): type of command law type for generic controller only
        nwfs (np.ndarray[ndim=1, dtype=np.int32]): index of wfss in controller
        nvalid (int): number of valid subaps
        nslope (int): number of slope to handle
        nslope_buffer (int): number of previous slopes to use in control
        ndm (np.ndarray[ndim=1, dtype=np.int32]): index of dms in controller
        nactu (int): number of controled actuator
        imat (np.ndarray[ndim=2,dtype=np.float32_t]): full interaction matrix
        cmat (np.ndarray[ndim=2,dtype=np.float32_t]): full control matrix
        maxcond (float): max condition number
        TTcond (float): tiptilt condition number for cmat filtering with mv controller
        delay (float): loop delay [frames]
        gain (float): loop gain
        nkl (int): number of KL modes used in imat_kl and used for computation of covmat in case of minimum variance controller
        cured_ndivs (int): subdivision levels in cured
        modopti (bool): Flag for modal optimization
        nrec (int): Number of sample of open loop slopes for modal optimization computation
        nmodes (int): Number of modes for M2V matrix (modal optimization)
        nmode_buffer (int): Number of previous modal vectors to use for control
        gmin (float): Minimum gain for modal optimization
        gmax (float): Maximum gain for modal optimization
        ngain (int): Number of tested gains
        do_kl_imat (bool): set imat kl on-off
        klpush (np.ndarray[ndim=1, dtype=np.float32]): Modal pushes values
        klgain (np.ndarray[ndim=1, dtype=np.float32]): Gain applied to modes at cMat inversion
        nstates (int): Number of states for generic linear controller
        nstate_buffer (int): Number of state vectors to use for control
        close_opti (bool): Flag for modal optimization with close
        mgain_init (float): Initial values of the modal gains
        lfdownup (tuple): Modal gain correction learning factor
        close_learning_factor (float): Autocorrelation learning factor
        close_target (float): Update framerate
        close_update_index (int): Target value
        n_iir_in (int): number of input taps to iir filter
        n_iir_out (int): number of output taps to iir filter
        polc (int): flag to do polc in generic linear controller
        modal (int): flag to use a modal control in generic linenar controller
        kernconv4imat (int): Flag to use kernel convolution when computing imat
        calpix_name (str): topic name of calpix stream
        loopdata_name (str): topic name of calpix stream
    """

    def __init__(self):
        self.__type = None  # type of controller
        self.__command_law = None  # type of command law type for generic controller only
        self.__nwfs = None  # index of wfss in controller
        self.__nvalid = 0  # number of valid subaps
        self.__nslope = 0  # number of slope to handle
        self.__nslope_buffer = 1  # number of previous slopes to use in control
        self.__ndm = None  # index of dms in controller
        self.__nactu = 0  # number of controled actuator
        self.__imat = None  # full interaction matrix
        self.__cmat = None  # full control matrix
        self.__maxcond = None  # max condition number
        self.__TTcond = None  # tiptilt condition number for cmat filtering with mv controller
        self.__delay = None  # loop delay [frames]
        self.__gain = None  # loop gain
        self.__nkl = None  # number of KL modes used in imat_kl and used for computation of covmat in case of minimum variance controller
        self.__cured_ndivs = None  # subdivision levels in cured

        """ MODAL OPTIMIZATION (style Gendron Lena 1994)"""
        self.__modopti = False  # Flag for modal optimization
        self.__nrec = (
            2048  # Number of sample of open loop slopes for modal optimization computation
        )
        self.__nmodes = None  # Number of modes for M2V matrix (modal optimization)
        self.__nmode_buffer = 0  # Number of previous modal vectors to use for control
        self.__gmin = 0.0  # Minimum gain for modal optimization
        self.__gmax = 1.0  # Maximum gain for modal optimization
        self.__ngain = 15  # Number of tested gains

        """ KL (actually BTT) BASIS INITIALIZATION """
        self.__do_kl_imat = False  # set imat kl on-off
        self.__klpush = None  # Modal pushes values
        self.__klgain = None  # Gain applied to modes at cMat inversion
        self.__nstates = 0  # Number of states for generic linear controller
        self.__nstate_buffer = 0  # Number of state vectors to use for control

        """ MODAL OPTIMIZATION CLOSE"""
        self.__close_opti = False  # Flag for modal optimization with close
        self.__mgain_init = 1.0  # Initial values of the modal gains
        self.__lfdownup = (0.01, 0.01)  # Modal gain correction learning factor
        self.__close_learning_factor = 0.3  # Autocorrelation learning factor
        self.__close_target = 0.0  # Update framerate
        self.__close_update_index = 1  # Target value
        self.__n_iir_in = 0  # number of input taps to iir filter
        self.__n_iir_out = 0  # number of output taps to iir filter
        self.__polc = 0  # flag to do polc in generic linear controller
        self.__modal = 0  # flag to use a modal control in generic linenar controller
        self.__kernconv4imat = 1  # Flag to use kernel convolution when computing imat
        self.__calpix_name = "compass_calPix"  # topic name of calpix stream
        self.__loopdata_name = "compass_loopData"  # topic name of loopData stream

    def get_calpix_name(self) -> str:
        """Get the topic name of calpix stream.

        Returns:
            str: The topic name of calpix stream.
        """
        return self.__calpix_name

    def set_calpix_name(self, calpix_name: str) -> None:
        """Set the calpix topic name.

        Args:
            calpix_name (str): The topic name of calpix stream.
        """
        self.__calpix_name = calpix_name

    calpix_name: str = property(get_calpix_name, set_calpix_name)

    def get_loopdata_name(self) -> str:
        """Get the topic name of calpix stream

        Returns:
            str: The topic name of calpix stream
        """
        return self.__loopdata_name

    def set_loopdata_name(self, loopdata_name: str) -> None:
        """Set the loop data topic name

        Args:
            loopdata_name (str): The topic name of loopData stream
        """
        self.__loopdata_name = loopdata_name

    loopdata_name: str = property(get_loopdata_name, set_loopdata_name)

    def get_type(self) -> str:
        """Get the controller type

        Returns:
            str: The type of controller
        """
        return self.__type

    def set_type(self, t: str) -> None:
        """Set the controller type

        Args:
            t (str): The type of controller
        """
        self.__type = scons.check_enum(scons.ControllerType, t)

    type: str = property(get_type, set_type)

    def get_command_law(self) -> str:
        """Get the command law type for generic controller only

        Returns:
            str: Command law type
        """
        return self.__command_law

    def set_command_law(self, t: str) -> None:
        """Set the command law type for generic controller only

        Args:
            t (str): Command law type
        """
        self.__command_law = scons.check_enum(scons.CommandLawType, t)

    command_law: str = property(get_command_law, set_command_law)

    def get_do_kl_imat(self) -> bool:
        """Get type imat, for imat on kl set at 1

        Returns:
            bool: imat kl
        """
        return self.__do_kl_imat

    def set_do_kl_imat(self, n: bool) -> None:
        """Set type imat, for imat on kl set at 1

        Args:
            n (bool): imat kl
        """
        self.__do_kl_imat = csu.enforce_or_cast_bool(n)

    do_kl_imat: bool = property(get_do_kl_imat, set_do_kl_imat)

    def get_klpush(self) -> np.ndarray:
        """Get klpush for imatkl size = number of kl mode

        Returns:
            np.ndarray: klpush
                Array of shape (nkl,) containing klpush values.
        """
        return self.__klpush

    def set_klpush(self, klpush: np.ndarray) -> None:
        """Set klpush for imatkl size = number of kl mode

        Args:
            klpush (np.ndarray): klpush
                Array of shape (nkl,) containing klpush values.
        """
        self.__klpush = csu.enforce_array(klpush, len(klpush), dtype=np.float32)

    klpush: np.ndarray = property(get_klpush, set_klpush)

    def get_klgain(self) -> np.ndarray:
        """Get klgain for imatkl size = number of kl mode

        Returns:
            np.ndarray: klgain
                Array of shape (nkl,) containing klgain values.
        """
        return self.__klgain

    def set_klgain(self, klgain: np.ndarray) -> None:
        """Set klgain for imatkl size = number of kl mode

        Args:
            klgain (np.ndarray): klgain
                Array of shape (nkl,) containing klgain values.
        """
        self.__klgain = csu.enforce_array(klgain, len(klgain), dtype=np.float32)

    klgain: np.ndarray = property(get_klgain, set_klgain)

    def get_nkl(self) -> int:
        """Get the number of KL modes used in imat_kl and used for computation of covmat in case of minimum variance controller

        Returns:
            int: The number of KL modes
        """
        return self.__nkl

    def set_nkl(self, n: int) -> None:
        """Set the number of KL modes used in imat_kl and used for computation of covmat in case of minimum variance controller

        Args:
            n (int): The number of KL modes
        """
        self.__nkl = csu.enforce_int(n)

    nkl: int = property(get_nkl, set_nkl)

    def get_nwfs(self) -> np.ndarray:
        """Get the indices of wfs

        Returns:
            np.ndarray[ndim=1, dtype=np.int32]: The indices of wfs
        """
        return self.__nwfs

    def set_nwfs(self, index_wfs: np.ndarray) -> None:
        """Set the indices of wfs

        Args:
            index_wfs (np.ndarray[ndim=1, dtype=np.int32]): The indices of wfs
        """
        self.__nwfs = csu.enforce_array(
            index_wfs, len(index_wfs), dtype=np.int32, scalar_expand=False
        )

    nwfs: np.ndarray = property(get_nwfs, set_nwfs)

    def get_ndm(self) -> np.ndarray:
        """Get the indices of dms

        Returns:
            np.ndarray[ndim=1, dtype=np.int32]: indices of dms
        """
        return self.__ndm

    def set_ndm(self, index_dm: np.ndarray) -> None:
        """Set the indices of dms

        Args:
            index_dm (np.ndarray[ndim=1, dtype=np.int32]): indices of dms
        """
        self.__ndm = csu.enforce_array(index_dm, len(index_dm), dtype=np.int32, scalar_expand=False)

    ndm: np.ndarray = property(get_ndm, set_ndm)

    def get_nactu(self) -> int:
        """Get the number of actuators

        Returns:
            int: number of actuators
        """
        return self.__nactu

    def set_nactu(self, nb_actu: int) -> None:
        """Set the number of actuators

        Args:
            nb_actu (int): number of actuators
        """
        self.__nactu = csu.enforce_int(nb_actu)

    nactu: int = property(get_nactu, set_nactu)

    def get_nslope(self) -> int:
        """Get the number of slopes

        Returns:
            int: The number of slopes
        """
        return self.__nslope

    def set_nslope(self, nb_slopes: int) -> None:
        """Set the number of slopes

        Args:
            nb_slopes (int): The number of slopes
        """
        self.__nslope = csu.enforce_int(nb_slopes)

    nslope: int = property(get_nslope, set_nslope)

    def get_nslope_buffer(self) -> int:
        """Get the number of slope buffers

        Returns:
            int: The number of slope buffers
        """
        return self.__nslope_buffer

    def set_nslope_buffer(self, nb_slopes_buffer: int) -> None:
        """Set the number of slope buffers

        Args:
            nb_slopes_buffer (int): The number of slope buffers
        """
        self.__nslope_buffer = csu.enforce_int(nb_slopes_buffer)

    nslope_buffer: int = property(get_nslope_buffer, set_nslope_buffer)

    def get_nvalid(self) -> int:
        """Get the number of valid subaps

        Returns:
            int: The number of valid subaps
        """
        return self.__nvalid

    def set_nvalid(self, nvalid: int) -> None:
        """Set the number of valid subaps

        Args:
            nvalid (int): The number of valid subaps
        """
        self.__nvalid = csu.enforce_int(nvalid)

    nvalid: int = property(get_nvalid, set_nvalid)

    def get_maxcond(self) -> float:
        """Get the max condition number

        Returns:
            float: The max condition number
        """
        return self.__maxcond

    def set_maxcond(self, maxcond: float) -> None:
        """Set the max condition number

        Args:
            maxcond (float): The max condition number
        """
        self.__maxcond = csu.enforce_float(maxcond)

    maxcond: float = property(get_maxcond, set_maxcond)

    def get_TTcond(self) -> float:
        """Get the tiptilt condition number for cmat filtering with mv controller.

        Returns:
            float: The tiptilt condition number.
        """
        return self.__TTcond

    def set_TTcond(self, ttcond: float) -> None:
        """Set the tiptilt condition number for cmat filtering with mv controller.

        Args:
            ttcond (float): The tiptilt condition number.
        """
        self.__TTcond = csu.enforce_float(ttcond)

    TTcond: float = property(get_TTcond, set_TTcond)

    def get_delay(self) -> float:
        """Get the loop delay expressed in frames.

        Returns:
            float: The loop delay [frames].
        """
        return self.__delay

    def set_delay(self, delay: float) -> None:
        """Set the loop delay expressed in frames.

        Args:
            delay (float): The loop delay [frames].
        """
        self.__delay = csu.enforce_float(delay)

    delay: float = property(get_delay, set_delay)

    def get_gain(self) -> float:
        """Get the loop gain.

        Returns:
            float: The loop gain.
        """
        return self.__gain

    def set_gain(self, gain: float) -> None:
        """Set the loop gain.

        Args:
            gain (float): The loop gain.
        """
        self.__gain = csu.enforce_float(gain)

    gain: float = property(get_gain, set_gain)

    def get_cured_ndivs(self) -> int:
        """Get the subdivision levels in cured.

        Returns:
            int: The subdivision levels in cured.
        """
        return self.__cured_ndivs

    def set_cured_ndivs(self, ndivs: int) -> None:
        """Set the subdivision levels in cured.

        Args:
            ndivs (int): The subdivision levels in cured.
        """
        self.__cured_ndivs = csu.enforce_int(ndivs)

    cured_ndivs: int = property(get_cured_ndivs, set_cured_ndivs)

    def get_modopti(self) -> bool:
        """Get the flag for modal optimization.

        Returns:
            bool: Flag for modal optimization.
        """
        return self.__modopti

    def set_modopti(self, modopti: bool) -> None:
        """Set the flag for modal optimization.

        Args:
            modopti (bool): Flag for modal optimization.
        """
        self.__modopti = csu.enforce_or_cast_bool(modopti)

    modopti: bool = property(get_modopti, set_modopti)

    def get_nrec(self) -> int:
        """Get the number of samples of open loop slopes for modal optimization computation.

        Returns:
            int: Number of samples.
        """
        return self.__nrec

    def set_nrec(self, nrec: int) -> None:
        """Set the number of samples of open loop slopes for modal optimization computation.

        Args:
            nrec (int): Number of samples.
        """
        self.__nrec = csu.enforce_int(nrec)

    nrec: int = property(get_nrec, set_nrec)

    def get_nmodes(self) -> int:
        """Get the number of modes for M2V matrix (modal optimization)

        Returns:
            int: Number of modes
        """
        return self.__nmodes

    def set_nmodes(self, nmodes: int) -> None:
        """Set the number of modes for M2V matrix (modal optimization)

        Args:
            nmodes (int): Number of modes
        """
        self.__nmodes = csu.enforce_int(nmodes)

    nmodes: int = property(get_nmodes, set_nmodes)

    def get_nmode_buffer(self) -> int:
        """Get the number of mode buffers

        Returns:
            int: Number of mode buffers
        """
        return self.__nmode_buffer

    def set_nmode_buffer(self, nmode_buffer: int) -> None:
        """Set the number of mode buffers

        Args:
            nmode_buffer (int): Number of mode buffers
        """
        self.__nmode_buffer = csu.enforce_int(nmode_buffer)

    nmode_buffer: int = property(get_nmode_buffer, set_nmode_buffer)

    def get_gmin(self) -> float:
        """Get the minimum gain for modal optimization.

        Returns:
            float: The minimum gain for modal optimization.
        """
        return self.__gmin

    def set_gmin(self, gmin: float) -> None:
        """Set the minimum gain for modal optimization.

        Args:
            gmin (float): The minimum gain for modal optimization.
        """
        self.__gmin = csu.enforce_float(gmin)

    gmin: float = property(get_gmin, set_gmin)

    def get_gmax(self) -> float:
        """Get the maximum gain for modal optimization.

        Returns:
            float: The maximum gain for modal optimization.
        """
        return self.__gmax

    def set_gmax(self, gmax: float) -> None:
        """Set the maximum gain for modal optimization.

        Args:
            gmax (float): The maximum gain for modal optimization.
        """
        self.__gmax = csu.enforce_float(gmax)

    gmax: float = property(get_gmax, set_gmax)

    def get_ngain(self) -> int:
        """Get the number of tested gains.

        Returns:
            int: The number of tested gains.
        """
        return self.__ngain

    def set_ngain(self, ngain: int) -> None:
        """Set the number of tested gains.

        Args:
            ngain (int): The number of tested gains.
        """
        self.__ngain = csu.enforce_int(ngain)

    ngain: int = property(get_ngain, set_ngain)

    def get_imat(self) -> np.ndarray:
        """Get the full interaction matrix.

        Returns:
            np.ndarray: Full interaction matrix.
        """
        return self.__imat

    def set_imat(self, imat: np.ndarray) -> None:
        """Set the full interaction matrix.

        Args:
            imat (np.ndarray): Full interaction matrix.
        """
        self.__imat = csu.enforce_arrayMultiDim(
            imat,
            (self.nslope, -1),  # Allow nModes or nActu as second dimension
            dtype=np.float32,
        )

    _imat: np.ndarray = property(get_imat, set_imat)

    def get_cmat(self) -> np.ndarray:
        """Get the full control matrix.

        Returns:
            np.ndarray: Full control matrix.
        """
        return self.__cmat

    def set_cmat(self, cmat: np.ndarray) -> None:
        """Set the full control matrix.

        Args:
            cmat (np.ndarray): Full control matrix.
        """
        self.__cmat = csu.enforce_arrayMultiDim(cmat, (self.nactu, self.nslope), dtype=np.float32)

    _cmat: np.ndarray = property(get_cmat, set_cmat)

    def get_nstates(self) -> int:
        """Get the number of states.

        Returns:
            int: Number of states.
        """
        return self.__nstates

    def set_nstates(self, nstates: int) -> None:
        """Set the number of states.

        Args:
            nstates (int): Number of states.
        """
        self.__nstates = csu.enforce_int(nstates)

    nstates: int = property(get_nstates, set_nstates)

    def get_nstate_buffer(self) -> int:
        """Get the number of state buffer.

        Returns:
            int: The number of state buffer.
        """
        return self.__nstate_buffer

    def set_nstate_buffer(self, nstate_buffer: int) -> None:
        """Set the number of state buffer.

        Args:
            nstate_buffer (int): The number of state buffer.
        """
        self.__nstate_buffer = csu.enforce_int(nstate_buffer)

    nstate_buffer: int = property(get_nstate_buffer, set_nstate_buffer)

    def get_close_opti(self) -> bool:
        """Get the flag for CLOSE modal optimization.

        Returns:
            bool: The CLOSE flag.
        """
        return self.__close_opti

    def set_close_opti(self, close_opti: bool) -> None:
        """Set the flag for CLOSE modal optimization.

        Args:
            close_opti (bool): The CLOSE flag.
        """
        self.__close_opti = close_opti

    close_opti: bool = property(get_close_opti, set_close_opti)

    def get_mgain_init(self) -> float:
        """Get the initial value of modal gains.

        Returns:
            float: The initial value for modal gains.
        """
        return self.__mgain_init

    def set_mgain_init(self, mgain_init: float) -> None:
        """Set the initial value of modal gains.

        Args:
            mgain_init (float): The initial value of modal gains.
        """
        self.__mgain_init = csu.enforce_float(mgain_init)

    mgain_init: float = property(get_mgain_init, set_mgain_init)

    def get_lfdownup(self) -> tuple:
        """Get the autocorrelation learning factors.

        Returns:
            tuple: The learning factors for autocorrelation.
        """
        return self.__lfdownup

    def set_lfdownup(self, qminus: float, qplus: float) -> None:
        """Set the autocorrelation learning factors.

        Args:
            qminus (float): The learning factor when higher than the target.
            qplus (float): The learning factor when lower than the target.
        """
        self.__lfdownup = (csu.enforce_float(qminus), csu.enforce_float(qplus))

    lfdownup: tuple[float, float] = property(get_lfdownup, set_lfdownup)

    def get_close_learning_factor(self) -> float:
        """Get the modal gain learning factor.

        Returns:
            float: The learning factor for modal gain.
        """
        return self.__close_learning_factor

    def set_close_learning_factor(self, close_learning_factor: float) -> None:
        """Set the modal gain optimization learning factor.

        Args:
            close_learning_factor (float): The learning factor.
        """
        self.__close_learning_factor = csu.enforce_float(close_learning_factor)

    lf: float = property(get_close_learning_factor, set_close_learning_factor)

    def get_close_target(self) -> float:
        """Get the autocorrelation target.

        Returns:
            float: The CLOSE autocorrelation target.
        """
        return self.__close_target

    def set_close_target(self, close_target: float) -> None:
        """Set the autocorrelation target.

        Args:
            close_target (float): The close target.
        """
        self.__close_target = csu.enforce_float(close_target)

    close_target: float = property(get_close_target, set_close_target)

    def get_close_update_index(self) -> int:
        """Get the modal gains update rate.

        Returns:
            int: The CLOSE update index.
        """
        return self.__close_update_index

    def set_close_update_index(self, idx: int) -> None:
        """Set the modal gains update rate.

        Args:
            idx (int): The CLOSE update index.
        """
        self.__close_update_index = csu.enforce_int(idx)

    close_update_index: int = property(get_close_update_index, set_close_update_index)

    def get_n_iir_in(self) -> int:
        """Get the number of inputs used in the IIR filter.

        Returns:
            int: The number of IIR inputs.
        """
        return self.__n_iir_in

    def set_n_iir_in(self, n_iir_in: int) -> None:
        """Set the number of inputs used in the IIR filter.

        Args:
            n_iir_in (int): The number of IIR inputs.
        """
        self.__n_iir_in = csu.enforce_int(n_iir_in)

    n_iir_in: int = property(get_n_iir_in, set_n_iir_in)

    def get_n_iir_out(self) -> int:
        """Get the number of outputs used in the IIR filter.

        Returns:
            int: The number of IIR outputs.
        """
        return self.__n_iir_out

    def set_n_iir_out(self, n_iir_out: int) -> None:
        """Set the number of outputs used in the IIR filter.

        Args:
            n_iir_out (int): The number of IIR outputs.
        """
        self.__n_iir_out = csu.enforce_int(n_iir_out)

    n_iir_out: int = property(get_n_iir_out, set_n_iir_out)

    def get_polc(self) -> bool:
        """Get POLC flag (True means using POL slopes)

        Returns:
            bool: POLC flag
        """
        return self.__polc

    def set_polc(self, polc: bool) -> None:
        """Set POLC flag (True means using POL slopes)

        Args:
            polc (bool): POLC flag
        """
        self.__polc = csu.enforce_or_cast_bool(polc)

    polc: bool = property(get_polc, set_polc)

    def get_modal(self):
        """ Get flag to use modal control \n(allows MVM from modes to actu)

        :return: (bool) : modal flag
        """
        return self.__modal

    def set_modal(self, modal):
        """ Set flag to use modal control \n(allows MVM from modes to actu)

        :param modal: (bool) : modal flag
        """
        self.__modal = csu.enforce_or_cast_bool(modal)

    modal = property(get_modal, set_modal)

    def get_kernconv4imat(self) -> int:
        """Get the value of kernconv4imat.

        Returns:
            int: The value of kernconv4imat.
        """
        return self.__kernconv4imat

    def set_kernconv4imat(self, kernconv4imat: int) -> None:
        """Set the value of kernconv4imat.

        Args:
            kernconv4imat (int): The value of kernconv4imat.
        """
        self.__kernconv4imat = kernconv4imat

    kernconv4imat: int = property(get_kernconv4imat, set_kernconv4imat)

    def get_kernconv4imat(self) -> int:
        """Get kernconv4imat, i.e. a flag for using kernel convolution to have better
        sensitivity on SH spot movements for imat computation

        Returns:
            int: The value of kernconv4imat
        """
        return self.__kernconv4imat

    def set_kernconv4imat(self, kernconv4imat: int) -> None:
        """Set kernconv4imat, i.e. a flag for using kernel convolution to have better
        sensitivity on SH spot movements for imat computation

        Args:
            kernconv4imat (int): The value of kernconv4imat
        """
        self.__kernconv4imat = csu.enforce_or_cast_bool(kernconv4imat)

    kernconv4imat: int = property(get_kernconv4imat, set_kernconv4imat)
