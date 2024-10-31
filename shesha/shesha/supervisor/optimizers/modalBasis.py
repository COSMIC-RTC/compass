#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
#
# COMPASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COMPASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COMPASS. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2024 COSMIC Team

from shesha.ao import basis
import shesha.util.utilities as util
import shesha.util.tools as tools
import shesha.util.make_pupil as mkP
import shesha.util.modesDM as modes
import shesha.constants as scons
import scipy.ndimage
from scipy.sparse import csr_matrix
import numpy as np


class ModalBasis(object):
    """This optimizer class handles all the modal basis and DM Influence functions
    related operations.

    Attributes:
        _config : (config) : Configuration parameters module

        _dms : (DmCompass) : DmCompass instance

        _target : (TargetCompass) : TargetCompass instance

        slaved_actus : TODO : docstring

        selected_actus : TODO : docstring

        couples_actus : TODO : docstring

        index_under_spiders : TODO : docstring

        modal_basis : (np.ndarray) : Last modal basis computed

        projection_matrix : (np.ndarray) : Last projection_matrix computed
    """

    def __init__(self, config, dms, target):
        """Instantiate a ModalBasis object

        Args:
            config : (config) : Configuration parameters module

            dms : (DmCompass) : DmCompass instance

            target : (TargetCompass) : TargetCompass instance
        """
        self._config = config
        self._dms = dms
        self._target = target
        self.slaved_actus = None
        self.selected_actus = None
        self.couples_actus = None
        self.index_under_spiders = None
        self.modal_basis = None
        self.projection_matrix = None

    def compute_influ_basis(self, dm_index: int) -> csr_matrix:
        """Computes and return the influence function phase basis of the specified DM
        as a sparse matrix

        Args:
            dm_index : (int) : Index of the DM

        Returns:
            influ_sparse : (csr_matrix) : influence function phases
        """
        return basis.compute_dm_basis(
            self._dms._dms.d_dms[dm_index],
            self._config.p_dms[dm_index],
            self._config.p_geom,
        )

    def compute_influ_delta(self, dm_index: int) -> np.ndarray:
        """Computes and return IF delta for the specified DM

        Args:
            dm_index : (int) : Index of the DM

        Return:
            influ_delta : (np.ndarray) : influence function deltas
        """
        ifsparse = basis.compute_dm_basis(
            self._dms._dms.d_dms[dm_index],
            self._config.p_dms[dm_index],
            self._config.p_geom,
        )
        mpup = self._config.p_geom.get_mpupil()
        # npix_in_pup = np.sum(mpup)
        ifdelta = ifsparse.dot(ifsparse.T) / np.sum(mpup)
        return ifdelta.toarray()

    def compute_modes_to_volts_basis(
        self,
        modal_basis_type: str,
        *,
        merged: bool = False,
        nbpairs: int = None,
        return_delta: bool = False,
    ) -> np.ndarray:
        """Computes a given modal basis ("KL2V", "Btt", "Btt_petal") and return the 2 transfer matrices

        Args:
            modal_basis_type : (str) : modal basis to compute ("KL2V", "Btt", "Btt_petal")

        Kwargs:
            merged : (bool) : TODO description

            nbpairs : (int) : TODO description

        Returns:
            modal_basis : (np.ndarray) : modes to volts matrix

            projection_matrix : (np.ndarray) : volts to modes matrix (None if "KL")
        """
        if modal_basis_type == "KL2V":
            print("Computing KL2V basis...")
            self.modal_basis = basis.compute_KL2V(
                self._config.p_controllers[0],
                self._dms._dms,
                self._config.p_dms,
                self._config.p_geom,
                self._config.p_atmos,
                self._config.p_tel,
            )
            fnz = util.first_non_zero(self.modal_basis, axis=0)
            # Computing the sign of the first non zero element
            # sig = np.sign(modal_basis[[fnz, np.arange(modal_basis.shape[1])]])
            sig = np.sign(
                self.modal_basis[tuple([fnz, np.arange(self.modal_basis.shape[1])])]
            )  # pour remove le future warning!
            self.modal_basis *= sig[None, :]
            # projection_matrix = None
        elif modal_basis_type == "Btt":
            print("Computing Btt basis...")
            self.modal_basis, self.projection_matrix = self.compute_btt_basis(
                merged=merged, nbpairs=nbpairs, return_delta=return_delta
            )
            fnz = util.first_non_zero(self.modal_basis, axis=0)
            # Computing the sign of the first non zero element
            # sig = np.sign(modal_basis[[fnz, np.arange(modal_basis.shape[1])]])
            sig = np.sign(
                self.modal_basis[tuple([fnz, np.arange(self.modal_basis.shape[1])])]
            )  # pour remove le future warning!
            self.modal_basis *= sig[None, :]
        elif modal_basis_type == "Btt_petal":
            print("Computing Btt with a Petal basis...")
            self.modal_basis, self.projection_matrix = self.compute_btt_petal()
        else:
            raise RuntimeError("Unsupported modal basis")

        return self.modal_basis, self.projection_matrix

    def compute_btt_basis(
        self,
        *,
        merged: bool = False,
        nbpairs: int = None,
        return_delta: bool = False,
    ) -> np.ndarray:
        """Computes the so-called Btt modal basis. The <merged> flag allows merto merge
        2x2 the actuators influence functions for actuators on each side of the spider (ELT case)

        Kwargs:
            merged : (bool) : If True, merge 2x2 the actuators influence functions for
                                        actuators on each side of the spider (ELT case). Default
                                        is False

            nbpairs : (int) : Default is None. TODO : description

            return_delta : (bool) : If False (default), the function returns
                                              Btt (modes to volts matrix),
                                              and P (volts to mode matrix).
                                              If True, returns delta = IF.T.dot(IF) / N
                                              instead of P

        Returns:
            Btt : (np.ndarray) : Btt modes to volts matrix

            projection_matrix : (np.ndarray) : volts to Btt modes matrix
        """
        from shesha.ao import basis

        dms_basis = basis.compute_IFsparse(self._dms._dms, self._config.p_dms, self._config.p_geom)
        influ_basis = dms_basis[:-2, :]
        tt_basis = dms_basis[-2:, :].toarray()
        if merged:
            couples_actus, index_under_spiders = self.compute_merged_influ(0, nbpairs=nbpairs)
            influ_basis2 = influ_basis.copy()
            index_remove = index_under_spiders.copy()
            index_remove += list(couples_actus[:, 1])
            print("Pairing Actuators...")
            for i in range(couples_actus.shape[0]):
                influ_basis2[couples_actus[i, 0], :] += influ_basis2[couples_actus[i, 1], :]
            print("Pairing Done")
            boolarray = np.zeros(influ_basis2.shape[0], dtype=bool)
            boolarray[index_remove] = True
            self.slaved_actus = boolarray
            self.selected_actus = ~boolarray
            self.couples_actus = couples_actus
            self.index_under_spiders = index_under_spiders
            influ_basis2 = influ_basis2[~boolarray, :]
            influ_basis = influ_basis2

        self.btt, self.projection_matrix = basis.compute_btt(
            influ_basis.T, tt_basis.T, return_delta=return_delta
        )

        if merged:
            btt2 = np.zeros((len(boolarray) + 2, self.btt.shape[1]))
            btt2[np.r_[~boolarray, True, True], :] = self.btt
            btt2[couples_actus[:, 1], :] = btt2[couples_actus[:, 0], :]

            P2 = np.zeros((self.btt.shape[1], len(boolarray) + 2))
            P2[:, np.r_[~boolarray, True, True]] = self.projection_matrix
            P2[:, couples_actus[:, 1]] = P2[:, couples_actus[:, 0]]
            self.btt = btt2
            self.projection_matrix = P2

        return self.btt, self.projection_matrix

    def compute_merged_influ(self, dm_index: int, *, nbpairs: int = None) -> np.ndarray:
        """Used to compute merged IF from each side of the spider
        for an ELT case (Petalling Effect)

        Args:
            dm_index : (int) : DM index

        Kwargs:
            nbpairs : (int) : Default is None. TODO : description

        Returns:
            pairs : (np.ndarray) : TODO description

            discard : (list) : TODO description
        """
        p_geom = self._config.p_geom

        cent = p_geom.pupdiam / 2.0 + 0.5
        p_tel = self._config.p_tel
        p_tel.t_spiders = 0.51
        spup = (
            mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent, cent).astype(np.float32).T
        )

        p_tel.t_spiders = 0.0
        spup2 = (
            mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent, cent).astype(np.float32).T
        )

        spiders = spup2 - spup

        (spidersID, k) = scipy.ndimage.label(spiders)
        spidersi = util.pad_array(spidersID, p_geom.ssize).astype(np.float32)
        px_list_spider = [np.where(spidersi == i) for i in range(1, k + 1)]

        # DM positions in iPupil:
        dm_posx = self._config.p_dms[dm_index]._xpos - 0.5
        dm_posy = self._config.p_dms[dm_index]._ypos - 0.5
        dm_pos_mat = np.c_[dm_posx, dm_posy].T  # one actu per column

        pitch = self._config.p_dms[dm_index]._pitch
        discard = np.zeros(len(dm_posx), dtype=bool)
        pairs = []

        # For each of the k pieces of the spider
        for k, px_list in enumerate(px_list_spider):
            pts = np.c_[px_list[1], px_list[0]]  # x,y coord of pixels of the spider piece
            # line_eq = [a, b]
            # Which minimizes leqst squares of aa*x + bb*y = 1
            line_eq = np.linalg.pinv(pts).dot(np.ones(pts.shape[0]))
            aa, bb = line_eq[0], line_eq[1]

            # Find any point of the fitted line.
            # For simplicity, the intercept with one of the axes x = 0 / y = 0
            if np.abs(bb) < np.abs(aa):  # near vertical
                one_point = np.array([1 / aa, 0.0])
            else:  # otherwise
                one_point = np.array([0.0, 1 / bb])

            # Rotation that aligns the spider piece to the horizontal
            rotation = np.array([[-bb, aa], [-aa, -bb]]) / (aa**2 + bb**2) ** 0.5

            # Rotated the spider mask
            rotated_px = rotation.dot(pts.T - one_point[:, None])
            # Min and max coordinates along the spider length - to filter actuators that are on
            # 'This' side of the pupil and not the other side
            min_u, max_u = (
                rotated_px[0].min() - 5.0 * pitch,
                rotated_px[0].max() + 5.0 * pitch,
            )

            # Rotate the actuators
            rotated_actus = rotation.dot(dm_pos_mat - one_point[:, None])
            sel_good_side = (rotated_actus[0] > min_u) & (rotated_actus[0] < max_u)
            threshold = 0.05
            # Actuators below this piece of spider
            sel_discard = (np.abs(rotated_actus[1]) < threshold * pitch) & sel_good_side
            discard |= sel_discard

            # Actuator 'near' this piece of spider
            sel_pairable = (
                (np.abs(rotated_actus[1]) > threshold * pitch)
                & (np.abs(rotated_actus[1]) < 1.0 * pitch)
                & sel_good_side
            )

            pairable_index = np.where(sel_pairable)[0]  # Indices of these actuators
            u_coord = rotated_actus[
                0, sel_pairable
            ]  # Their linear coord along the spider major axis

            order = np.sort(u_coord)  # Sort by linear coordinate
            order_index = pairable_index[np.argsort(u_coord)]  # And keep track of original indexes

            # i = 0
            # while i < len(order) - 1:
            if nbpairs is None:
                i = 0
                ii = len(order) - 1
            else:
                i = len(order) // 2 - nbpairs
                ii = len(order) // 2 + nbpairs
            while i < ii:
                # Check if next actu in sorted order is very close
                # Some lonely actuators may be hanging in this list
                if np.abs(order[i] - order[i + 1]) < 0.2 * pitch:
                    pairs += [(order_index[i], order_index[i + 1])]
                    i += 2
                else:
                    i += 1
        print("To discard: %u actu" % np.sum(discard))
        print("%u pairs to slave" % len(pairs))
        if np.sum(discard) == 0:
            discard = []
        else:
            list(np.where(discard)[0])
        return np.asarray(pairs), list(np.where(discard)[0])

    def compute_btt_petal(self) -> np.ndarray:
        """Computes a Btt modal basis with Pistons filtered

        Returns:
            Btt : (np.ndarray) : Btt modes to volts matrix

            P : (np.ndarray) : volts to Btt modes matrix
        """
        pzt_index = np.where([d.type is scons.DmType.PZT for d in self._config.p_dms])[0][0]
        influ_pzt = self.compute_influ_basis(pzt_index)
        petal_dm_index = np.where(
            [d.influ_type is scons.InfluType.PETAL for d in self._config.p_dms]
        )[0][0]
        influ_petal = self.compute_influ_basis(petal_dm_index)
        tt_index = np.where([d.type is scons.DmType.TT for d in self._config.p_dms])[0][0]
        influ_tt = self.compute_influ_basis(tt_index).toarray()

        self.modal_basis, self.projection_matrix = basis.compute_btt(
            influ_pzt.T, influ_tt.T, influ_petal=influ_petal
        )
        return self.modal_basis, self.projection_matrix

    def compute_phase_to_modes(self, modal_basis: np.ndarray) -> np.ndarray:
        """Return the phase to modes matrix by using the given modal basis

        Args:
            modal_basis : (np.ndarray) : Modal basis matrix

        Returns:
            phase_to_modes : (np.ndarray) : phase to modes matrix
        """
        nbmode = modal_basis.shape[1]
        phase = self._target.get_tar_phase(0)
        phase_to_modes = np.zeros((nbmode, phase.shape[0], phase.shape[1]))
        S = np.sum(self._config.p_geom._spupil)
        for i in range(nbmode):
            self._dms.set_command((modal_basis[:, i]).copy())
            # self.next(see_atmos=False)
            self._target.raytrace(0, dms=self._dms, ncpa=False, reset=True)
            phase = self._target.get_tar_phase(0, pupil=True)
            # Normalisation pour les unites rms en microns !!!
            norm = np.sqrt(np.sum((phase) ** 2) / S)
            if norm == 0:
                norm = 1
            phase_to_modes[i] = phase / norm
        return phase_to_modes

    def compute_ipos_in_pupil(self, d_obs=11.4, d_pup=37., pixsize=None, xpos=None, ypos=None, n_pix=None):
        """
        Return actuators indexes located inside and outside of a set pupil. 
        The "in pupil" selected actuators are inside a ring of inner diameter d_obs and outer diameter d_pup.
        
        Input:
        pixsize : size of a pixel in the pupil. If None : take the value given by ADOPT
        xpos, ypos : actuators positions in pixels. If None : take the value given by ADOPT
        n_pix : size in pixels of the support where xpos and ypos are expressed
        d_obs : diameter of the inner ring in meters, defined by the pupil obstruction.
        d_pup : diameter of the outer ring in meters, defined by the telescope pupil diameter.

        Output:
        ipos_in : index of actuators inside the defined pupil
        ipos_out : index of actuators outside of the defined pupil (complement of ipos_pup)

        Example :
        ipos_in, ipos_out = basis.compute_ipos_in_pupil(d_obs=.28 * 40, d_pup=38.5)

        """

        if (xpos is None) and (ypos is None) and (n_pix is None):
            xpos = self._config.p_dms[0]._xpos
            ypos = self._config.p_dms[0]._ypos
            n_pix = self._config.p_geom._ipupil.shape[0]
            d_center = n_pix//2 - 0.5
        else:
            xpos = np.array(xpos)
            ypos = np.array(ypos)
            d_center = n_pix//2 - 0.5
        
        if (pixsize == None):
            pixsize = self._config.p_geom._pixsize
        
        dist = np.sqrt((xpos - d_center)**2 + (ypos - d_center)**2)
        dist *= pixsize
        ipos_in = np.where((dist < (d_pup/2)) * (dist > (d_obs/2)))[0]
        ipos_out = np.where(np.isin(np.arange(len(xpos)), ipos_in)==False)[0]
        
        return ipos_in, ipos_out

    def compute_ipos_spider(self, n_seg=6, dm_diam=42*1.2, d_spi=0.54, pixsize=None, xpos=None, ypos=None, n_pix=None):
        """
        find the index of actuators along spiders
        comes without any guarantee regarding bugs!

        Input:
        com : ADOPT command class
        ao : ADOPT ao class
        n_pix : number of pixels of the 
        pixsize : size of a pixel in the pupil. If None : take the value given by ADOPT
        xpos, ypos : actuators positions in pixels. If None : take the value given by ADOPT
        n_seg : number of fragment of the pupil
        dm_diam : diameter of the deformable mirror in meters
        d_spi : width of the spider arms in meters

        Output:
        ipos_spi1 : index of actuators located on one side of the spider arms
        ipos_spi2 : index of actuators located on the other side of the spider arms

        Example :
        # If running a simulation with COMPASS, ADOPT already know your DM so you don't have to provide pixsize, xpos, ypos and n_pix
        spi1, spi2 = compute_ipos_spider(com, ao, n_seg=6, dm_diam=42, d_spi=0.6)

        """
        if (xpos is None) and (ypos is None) and (n_pix is None):
            xpos = self._config.p_dms[0]._xpos
            ypos = self._config.p_dms[0]._ypos
            n_pix = self._config.p_geom._ipupil.shape[0]
        
        if (pixsize is None):
            pixsize = self._config.p_geom._pixsize
        
        print("Create 2 pupils with 1/2 spiders")
        dm_pup1 = self.create_pupil(n_pix, dm_diam / pixsize, n_seg=n_seg, d_spider = 1e-5, d_spider2 = 2*d_spi / pixsize)
        dm_pup2 = self.create_pupil(n_pix, dm_diam / pixsize, n_seg=n_seg, d_spider = 2*d_spi / pixsize, d_spider2 = 1e-5)

        print("Find actuators indexes along the spiders")
        spi1 = []
        spi2_temp = []

        for i in range(len(xpos)):
            if dm_pup1[int(xpos[i]), int(ypos[i])] ==0:
                spi1 += [i]
            elif dm_pup2[int(xpos[i]), int(ypos[i])] ==0:
                spi2_temp += [i]
        spi1 = np.array(spi1, dtype=np.int32)
        spi2_temp = np.array(spi2_temp, dtype=np.int32)

        x = np.linspace(-1, 1, n_pix)
        xx, yy = np.meshgrid(x, x)
        _, theta = tools.cart2polar(xx, yy)

        xpos = np.array(xpos, dtype=np.int32)
        ypos = np.array(ypos, dtype=np.int32)

        print("Sort the indexes to get compatibles pairs of actuators")
        r1 = np.argsort(np.diag(theta.T[xpos[spi1]][:, ypos[spi1]]))
        r2 = np.argsort(np.diag(theta.T[xpos[spi2_temp]][:, ypos[spi2_temp]]))

        spi1 = spi1[r1]
        spi2 = spi2_temp*0
        tt = len(spi2_temp)/n_seg
        if not tt.is_integer():
            print("Warning, number of actuators along spiders is not multiple of n_seg")
        tt = int(tt)
        for i in range(n_seg):
            spi2[i*tt:(i+1)*tt] = spi2_temp[r2[i*tt:(i+1)*tt][::-1]]
        
        return spi1, spi2


    def compute_Bg(self, ipos_in=None, IFdelta=None, tt_mode=False, pixsize=None, xpos=None, ypos=None):
        """
        Calcul de la base de Gendron : vecteurs propres de la matrice de covariance des distances inter actionneurs puissance 5/3
        
        Input:
        <ipos_in> : (1D np.arr) : Optional (default=None), actuators indexes inside the effective pupil. The number of modes will be equal to the number of given actuators.
        <IFdelta> : (2D np.arr) : Optional (default=None), but required for basis normalization, Influence Functions covariance matrix
        <tt_mode> : (bool) : Optional (default=False), if True filter the tip tilt from the M4 space and add these modes to the Tip Tilt mirror
        <pixsize> : (float) : Optional (default=None), size of a pixel in the pupil. If None : take the value given by ADOPT
        <xpos> : (1D np.arr) : Optional (default=None), actuators positions in pixels. If None : take the value given by ADOPT
        <ypos> : (1D np.arr) : Optional (default=None), actuators positions in pixels. If None : take the value given by ADOPT

        Output:
        Bg : (2D np.arr) : Gendon modal Basis with number of degrees of freedom defined over ipos_in
        Bgext : (2D np.arr) : Gendon modal Basis with 2 lines for the tip tilt mirror

        Example :
        # If running a simulation with COMPASS, ADOPT already know your DM so you don't have to provide pixsize, xpos, ypos
        ipos_in, ipos_out = compute_ipos_in_pupil(com, ao, d_obs=.28 * 40, d_pup=38.5)
        IFdelta = IFsp.dot(IFsp.T) / IFsp.shape[0]          # IFdelta tableau Nactu x Nactu
        Bg, Bgext = compute_Bg(ao, ipos_in=ipos_in, IFdelta=IFdelta)
        """

        if (xpos is None) and (ypos is None):
            xpos = self._config.p_dms[0]._xpos
            ypos = self._config.p_dms[0]._ypos
        else:
            xpos = np.array(xpos)
            ypos = np.array(ypos)
        
        if (pixsize is None):
            pixsize = self._config.p_geom._pixsize
        
        dist = np.sqrt((xpos[:, None] - xpos[None, :])**2 + (ypos[:, None] - ypos[None, :])**2)    # [pixels]
        dist *= pixsize    # [meters]

        n_actu = len(xpos)

        # Calcul des modes
        L0 = 1e4   # valeur fausse, mais proche de l'infini, donc pas grave.
        if ipos_in is not None:
            B, l = modes.KLmodes(xpos[ipos_in], ypos[ipos_in], L0, True)
            n_modes = len(ipos_in)
        else:
            B, l = modes.KLmodes(xpos, ypos, L0, True)
            n_modes = len(xpos)

        # Normalisation de la base 
        if IFdelta is not None:
            if ipos_in is not None:
                var = B.T.dot(IFdelta[ipos_in][:,ipos_in].dot(B))
                Bt = B / (np.sqrt(np.diag(var))[None,:])
                Bg = np.zeros((n_actu, n_modes))
                Bg[ipos_in,:] = Bt
            else:
                var = B.T.dot(IFdelta.dot(B))
                Bg = B / (np.sqrt(np.diag(var))[None,:])
        else:
            Bg = np.zeros((n_actu, n_modes))
            Bg[ipos_in,:] = B
        # On ajoute le Tip Tilt au début
        Bgext, itt, _ = modes.moreLines(Bg, 2)

        if tt_mode:
            Bgext[:, 0:2] = 0.0
            Bgext[np.ix_(itt, [0,1])] = np.eye(2)
            # Normalisation du Tip Tilt par une méthode """adéquate"""
            Bgext[itt, 0:2] *= 0.02

        return Bg, Bgext


    def compute_Br(self, Bg, ipos_in, ipos_out, L0=1e4, r0=0.1, alpha=0, IFdelta=None, tt_mode=False, pixsize=None, xpos=None, ypos=None):
        """
        Extend the modal basis to the actuators located outside of the pupil, so that the mode stays "Kolmo - compatible"
        
        Input:
        <Bg> : (2D np.arr) : Required, Gendron modal basis (or any other basis)
        <ipos_in> : (1D np.arr) : Required, actuators indexes inside the effective pupil and used as "master" to extrapolate
        <ipos_out> : (1D np.arr) : Required, actuators indexes to be extended
        <L0> : (float) : Optional (default=1e4), Outer scale in meters, can be set to add more or less high frequencies during extrapolation
        <r0> : (float) : Optional (default=0.1), Fried parameter in meters, take a physical value but it will not change the world order
        <alpha> : (float) : Optional (default=0), used for regularization -> ON GOING, keep alpha=0 !!
        <IFdelta> : (2D np.arr) : Optional (default=None) but required for basis normalization, Influence Functions covariance matrix
        <pixsize> : (float) : Optional (default=None), size of a pixel in the pupil. If None : take the value given by ADOPT
        <xpos> : (1D np.arr) : Optional (default=None), actuators positions in pixels. If None : take the value given by ADOPT
        <ypos> : (1D np.arr) : Optional (default=None), actuators positions in pixels. If None : take the value given by ADOPT

        Output:
        Br : (2D np.arr) : Gendon modal Basis with number of degrees of freedom defined over ipos_in, extended to all actuators, even the outer ring actuators
        Brext : (2D np.arr) : Gendon modal Basis with 2 lines for the tip tilt mirror
        
        Example :
        # If running a simulation with COMPASS, ADOPT already know your DM so you don't have to provide pixsize, xpos, ypos
        ipos_in, ipos_out = compute_ipos_in_pupil(com, ao, d_obs=.28 * 40, d_pup=38.5)
        IFdelta = IFsp.dot(IFsp.T) / IFsp.shape[0]          # IFdelta tableau Nactu x Nactu
        Br, Brext = compute_Bg(ao, Bg, ipos_in, ipos_out, IFdelta=IFdelta)
        """

        if (xpos is None) and (ypos is None):
            xpos = self._config.p_dms[0]._xpos
            ypos = self._config.p_dms[0]._ypos
        if (pixsize is None):
            pixsize = self._config.p_geom._pixsize
        
        dist = np.sqrt((xpos[:, None] - xpos[None, :])**2 + (ypos[:, None] - ypos[None, :])**2)    # distances entre actus [mètres]
        dist *= pixsize
        Br = modes.computeMmseMatrix(Bg, dist, ipos_out, ipos_in, L0=L0, r0 = r0, alpha = alpha)
        if IFdelta is not None:
            Br = self.normalize_basis(Br, IFdelta)
        Brext, itt, _ = modes.moreLines(Br, 2) 

        return Br, Brext


    def compute_Bp(ao, Br, spi1, spi2, IFdelta=None):
        """
        Pairing of the actuators at each edge of the spider arms

        Input:
        ao : Required, ADOPT ao class
        <Br> : (2D np.arr) : Required, Gendron modal basis extended to the ring (or any other basis)
        <ipos_spi1> : (1D np.arr) : Required, actuators indexes located on one side of the spider arms
        <ipos_spi2> : (1D np.arr) : Required, actuators indexes located on the other side of the spider arms
        <IFdelta> : (2D np.arr) : Optional (default=None) but required for basis normalization, Influence Functions covariance matrix
        Output:
        Bp : (2D np.arr) : Paired modal Basis
        Bpext : (2D np.arr) : Paired modal Basis with 2 lines for the tip tilt mirror
        """
        Bp = Br.copy()
        Bp[spi1] = Bp[spi2]
        if IFdelta is not None:
            Bp = normalize_basis(Bp, IFdelta)
        Bpext, itt, _ = modes.moreLines(Bp, 2)

        return Bp, Bpext


    def compute_Bc(ao, Br, spi1, L0=1e4, r0=0.1, alpha=0, IFdelta=None, pixsize=None, xpos=None, ypos=None):
        """
        Remove pure piston degrees of freedom to avoid petalling while keeping the modes "Kolmo compatible"

        Input:
        ao : Required, ADOPT ao class
        <Br> : (2D np.arr) : Required, Gendron modal basis extended to the ring (or any other basis)
        <ipos_spi1> : (1D np.arr) : Required, actuators indexes located on one side of the spider arms
        <L0> : (float) : Optional (default=1e4), Outer scale in meters, can be set to add more or less high frequencies during extrapolation
        <r0> : (float) : Optional (default=0.1), Fried parameter in meters, take a physical value but it will not change the world order
        <alpha> : (float) : Optional (default=0), used for regularization -> ON GOING, keep alpha=0 !!
        <IFdelta> : (2D np.arr) : Optional (default=None) but required for basis normalization, Influence Functions covariance matrix
        <pixsize> : (float) : Optional (default=None), size of a pixel in the pupil. If None : take the value given by ADOPT
        <xpos> : (1D np.arr) : Optional (default=None), actuators positions in pixels. If None : take the value given by ADOPT
        <ypos> : (1D np.arr) : Optional (default=None), actuators positions in pixels. If None : take the value given by ADOPT

        Output:
        Bc : (2D np.arr) : Continuous modal basis
        Bcext : (2D np.arr) : Continuous modal basis with 2 lines for the tip tilt mirror
        """
        
        if (xpos is None) and (ypos is None):
            xpos = self._config.p_dms[0]._xpos
            ypos = self._config.p_dms[0]._ypos
            n_pix = self._config.p_geom._ipupil.shape[0]
        if (pixsize is None):
            pixsize = self._config.p_geom._pixsize
        
        dist = np.sqrt((xpos[:, None] - xpos[None, :])**2 + (ypos[:, None] - ypos[None, :])**2)    # distances entre actus [mètres]
        dist *= pixsize
        ipos_all = np.arange(len(xpos))
        ipos_interieur = ipos_all[np.where(np.isin(ipos_all, spi1) == False)]
        Bc = mmse.computeMmseMatrix(Br, dist, spi1, ipos_interieur, L0=L0, r0=r0, alpha=alpha)
        if IFdelta is not None:
            Bc = normalize_basis(Bc, IFdelta)
        Bcext, itt, _ = modes.moreLines(Bc, 2)
        return Bc, Bcext


    def filter_mode(B, fmode, IFdelta):
        """
        Filter a mode from the basis

        Input:
        ao : Required, ADOPT ao class
        <B> : (2D np.arr) : Required, Modal basis
        <fmode> : (1D np.arr) : mode expressed over the actuator space to be filtered
        <IFdelta> : (2D np.arr) : Required, Influence Functions covariance matrix

        Output:
        Bf : (2D np.arr) : Filtered modal basis
        """
        Bf = B.copy()
        dd = np.linalg.inv(fmode.T.dot(IFdelta).dot(fmode))
        Bf -= fmode.dot(dd).dot(fmode.T.dot(IFdelta).dot(B))

        return Bf


    def normalize_basis(self, B, IFdelta):
        """
        Normalization of the modal basis in the phase space

        Input:
        <B> : (2D np.arr) : Required, Modal basis
        <IFdelta> : (2D np.arr) : Required, Influence Functions covariance matrix

        Output:
        Bf : (2D np.arr) : Filtered modal basis
        """
        # normalisation de la base
        print("Normalization ...")
        var = B.T.dot(IFdelta.dot(B))
        Bn = B / (np.sqrt(np.diag(var))[None,:])

        return Bn


    def control_unseen_actu(ao, cmat, pos_actu, L0=1e4, r0=0.2, pixsize=None, xpos=None, ypos=None):
        """
        Remove a given actuator as degree of freedom from the command matrix, in a way that it provides a "Kolmo compatible" command.

        Input:
        ao : Required, ADOPT ao class
        <cmat> : (2D np.arr) : Actuators command matrix
        <pos_actu> : (int, 1D np.arr) : Index or array of indexes of the actuators to be "mmse-ifier"
        <L0> : (float) : Optional (default=1e4), Outer scale in meters, can be set to add more or less high frequencies during extrapolation
        <r0> : (float) : Optional (default=0.1), Fried parameter in meters, take a physical value but it will not change the world order
        <alpha> : (float) : Optional (default=0), used for regularization -> ON GOING, keep alpha=0 !!
        <pixsize> : (float) : Optional (default=None), size of a pixel in the pupil. If None : take the value given by ADOPT
        <xpos> : (1D np.arr) : Optional (default=None), actuators positions in pixels. If None : take the value given by ADOPT
        <ypos> : (1D np.arr) : Optional (default=None), actuators positions in pixels. If None : take the value given by ADOPT

        Output:
        cmat_u : (2D np.arr) : New command matrix
        """
        
        if (xpos is None) and (ypos is None):
            xpos = self._config.p_dms[0]._xpos
            ypos = self._config.p_dms[0]._ypos
            n_pix = self._config.p_geom._ipupil.shape[0]
        if pixsize is None:
            pixsize = self._config.p_geom._pixsize
        dist = np.sqrt((xpos[:, None] - xpos[None, :])**2 + (ypos[:, None] - ypos[None, :])**2)    # distances entre actus [mètres]
        dist *= pixsize

        pos_all_actu = np.arange(len(xpos))  # indice de tous les actionneurs

        comp_pos_actu = pos_all_actu[np.where(np.isin(pos_all_actu, pos_actu) == False)]    # indice des actionneurs à ne pas mmse-er (à garder)
        cmat_u = mmse.computeMmseMatrix(cmat, dist, pos_actu, comp_pos_actu, L0=L0, r0=r0)

        return cmat_u


    def create_pupil(self, n_pix, d_tel, n_seg=0, d_obs=0, d_spider=0, d_spider2 = 0, obs_shape=None, custom_obs=None):
        """
        create a pupil with parameters
        all distances expressed in pixels with respect to n_pix
        """
        x = np.linspace(-1, 1, n_pix)
        xx, yy = np.meshgrid(x, x)
        r = np.sqrt(xx**2 + yy**2)

        from scipy.ndimage import rotate
        my_pup = (r < d_tel/n_pix) * 1

        if d_spider2 != 0.:
            my_pup[int(n_pix/2 - d_spider/2):int(n_pix/2 + d_spider2/2), n_pix//2:] = 0
        else:
            my_pup[int(n_pix/2 - d_spider/2):int(n_pix/2 + d_spider/2), n_pix//2:] = 0
        
        my_pup_rot = my_pup
        for i in range(n_seg):
            my_pup_rot = rotate(my_pup_rot, 360 / n_seg, reshape=False)
            my_pup *= my_pup_rot

        if obs_shape is "ELT":
            my_pup[np.where(r < d_obs / n_pix)] = custom_obs[np.where(r < d_obs / n_pix)]
        
        elif obs_shape is "hexa":
            centers = np.c_[np.cos((2 * np.arange(n_seg) + 1) * np.pi/n_seg), np.sin((2 * np.arange(n_seg) +1) * np.pi/n_seg)]
            h = np.abs(np.min(np.asarray([(c[0]) * xx + (c[1]) * yy for c in centers]), axis=0))

            my_pup[np.where(h < d_obs / n_pix)] = 0
        else:
            my_pup *= (r >= d_obs/n_pix) * 1

        return my_pup


    def compute_petal_basis(ao, n_pix, dm_pix, d_obs, n_seg, d_spider, xpos=None, ypos=None):
        """
        Create petal modes following given deformable mirror parameters. Not normalized !!

        Input:
        ao : Required, ADOPT ao class
        <n_pix> : (int) : Required, Deformable Mirror support size in pixels
        <dm_pix> : (float) : Required, DM diameter in pixels
        <d_obs> : (float) : Required, DM obstruction in pixels (d_obs < dm_pix)
        <n_seg> : (int) : Required, number of fragments of the DM
        <d_spider> : (float) : Required, width of the spider arms
        Output:
        modepetal : (2D np.arr) : Petal modal basis
        a : (2D np.arr) : Pupil morphology with indexed fragments
        """
        dm_pup = create_pupil(n_pix, dm_pix, n_seg=n_seg, d_spider = d_spider, d_obs = d_obs)
        a, ns = label(dm_pup)
        if ns != n_seg:
            print(ns)
            print("bug in the number of petals")

        if (xpos == None) and (ypos == None):
            xpos = np.array(ao.dm0.CsX, dtype=np.int32)
            ypos = np.array(ao.dm0.CsY, dtype=np.int32)

        mode_petal = np.zeros((ao.Nactu, n_seg))
        for i in range(n_seg):
            for k in range(ao.dm0.Nactu):
                if a[xpos[k], ypos[k]] == i+1:
                    mode_petal[k, i] = 1
        
        return mode_petal, a
