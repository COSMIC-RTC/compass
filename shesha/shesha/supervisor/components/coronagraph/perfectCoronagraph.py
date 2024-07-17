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

import numpy as np
import shesha.config as conf
import shesha.constants as scons
from shesha.supervisor.components.coronagraph.genericCoronagraph import (
    GenericCoronagraph,
)
from shesha.init.coronagraph_init import init_coronagraph, init_mft
from shesha.supervisor.components.targetCompass import TargetCompass
from sutra import PerfectCoronagraph
from carma import context


class PerfectCoronagraphCompass(GenericCoronagraph):
    """Class supervising perfect coronagraph component

    Attributes:
        _spupil: (np.ndarray[ndim=2, dtype=np.float32]): Telescope pupil mask

        _pupdiam : (int): Number of pixels along the pupil diameter

        _dim_image :(int): Coronagraphic image dimension

        _p_corono: (ParamCoronagraph): Coronagraph parameters

        _target: (TargetCompass): Compass Target used as input for the coronagraph

        _norm_img : (float): Normalization factor for coronagraphic image

        _norm_psf : (float): Normalization factor for PSF

        _coronagraph: (SutraCoronagraph): Sutra coronagraph instance

        _wav_vec: (np.ndarray[ndim=1, dtype=np.float32]): Vector of wavelength

        _AA: (np.ndarray[ndim=3, dtype=np.complex64]): MFT matrix for image computation

        _BB: (np.ndarray[ndim=3, dtype=np.complex64]): MFT matrix for image computation

        _norm0: (np.ndarray[ndim=3, dtype=np.complex64]): MFT matrix for image computation

        _AA_c: (np.ndarray[ndim=3, dtype=np.complex64]): MFT matrix for image computation

        _BB_c: (np.ndarray[ndim=3, dtype=np.complex64]): MFT matrix for psf computation

        _norm0_c: (np.ndarray[ndim=3, dtype=np.complex64]): MFT matrix for psf computation

        _indices_pup: (tuple): Tuple of ndarray containing X and Y indices of illuminated
                                pixels in the pupil
    """

    def __init__(
        self,
        context: context,
        targetCompass: TargetCompass,
        p_corono: conf.ParamCoronagraph,
        p_geom: conf.ParamGeom,
    ):
        """Initialize a perfect coronagraph instance

        Args:
            context: (carma.context): GPU context

            targetCompass: (TargetCompass): Compass Target used as input for the coronagraph

            p_corono: (ParamCoronagraph): Coronagraph parameters

            p_geom: (ParamGeom): Compass geometry parameters
        """
        init_coronagraph(p_corono, p_geom.pupdiam)
        GenericCoronagraph.__init__(self, p_corono, p_geom, targetCompass)
        self._wav_vec = p_corono._wav_vec

        self._AA, self._BB, self._norm0 = init_mft(p_corono, self._pupdiam, planes="lyot_to_image")
        self._AA_c, self._BB_c, self._norm0_c = init_mft(
            p_corono,
            self._pupdiam,
            planes="lyot_to_image",
            center_on_pixel=True,
        )
        self._indices_pup = np.where(self._spupil > 0.0)

        self._coronagraph = PerfectCoronagraph(
            context,
            self._target.sources[0],
            self._dim_image,
            self._dim_image,
            self._wav_vec,
            self._wav_vec.size,
            0,
        )

        self._coronagraph.set_mft(self._AA, self._BB, self._norm0, scons.MftType.IMG)
        self._coronagraph.set_mft(self._AA_c, self._BB_c, self._norm0_c, scons.MftType.PSF)
        self._compute_normalization()

    def _compute_normalization(self):
        """Computes the normalization factor of coronagraphic images (CPU based)"""
        self._target.reset_tar_phase(0)
        self.compute_psf(accumulate=False)
        self._norm_img = np.max(self.get_psf(expo_type=scons.ExposureType.SE))
        self._norm_psf = self._norm_img
