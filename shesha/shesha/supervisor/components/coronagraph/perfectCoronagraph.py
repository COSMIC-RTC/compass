import numpy as np
import shesha.config as conf
import shesha.constants as scons
from shesha.supervisor.components.coronagraph.genericCoronagraph import GenericCoronagraph
from shesha.init.coronagraph_init import init_coronagraph, init_mft, mft_multiplication
from shesha.supervisor.components.targetCompass import TargetCompass
from sutraWrap import PerfectCoronagraph
from carmaWrap import context


class PerfectCoronagraphCompass(GenericCoronagraph):
    """ Class supervising perfect coronagraph component

    Attributes:
        _image_se: (np.ndarray[ndim=2, dtype=np.float32]): Short exposure coronagraphic image

        _image_le: (np.ndarray[ndim=2, dtype=np.float32]): Long exposure coronagraphic image

        _psf_se: (np.ndarray[ndim=2, dtype=np.float32]): Long exposure PSF

        _psf_le: (np.ndarray[ndim=2, dtype=np.float32]): Long exposure PSF

        _spupil: (np.ndarray[ndim=2, dtype=np.float32]): Telescope pupil mask

        _pupdiam : (int): Number of pixels along the pupil diameter

        _dim_image :(int): Coronagraphic image dimension

        _p_corono: (Param_corono): Coronagraph parameters

        _target: (TargetCompass): Compass Target used as input for the coronagraph

        _norm_img : (float): Normalization factor for coronagraphic image

        _norm_psf : (float): Normalization factor for PSF

        _coronagraph: (SutraCoronagraph): Sutra coronagraph instance

        _wav_vec: (np.ndarray[ndim=1, dtype=np.float32]): Vector of wavelength

        _AA: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for image computation

        _BB: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for image computation

        _norm0: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for image computation

        _AA_c: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for image computation

        _BB_c: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for psf computation
        
        _norm0_c: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for psf computation

        _indices_pup: (tuple): Tuple of ndarray containing X and Y indices of illuminated 
                                pixels in the pupil
    """
    def __init__(self, context: context, targetCompass: TargetCompass, 
                 p_corono: conf.Param_corono, p_geom: conf.Param_geom):
        """ Initialize a perfect coronagraph instance

        Args:
            context: (CarmaWrap.context): GPU context

            targetCompass: (TargetCompass): Compass Target used as input for the coronagraph

            p_geom: (Param_geom): Compass geometry parameters

            p_corono: (Param_corono): Coronagraph parameters
        """
        init_coronagraph(p_corono, p_geom.pupdiam)
        GenericCoronagraph.__init__(self, p_corono, p_geom, targetCompass)
        self._wav_vec = p_corono._wav_vec

        self._AA, self._BB, self._norm0 = init_mft(p_corono,
                                                   self._pupdiam,
                                                   planes='lyot_to_image')
        self._AA_c, self._BB_c, self._norm0_c = init_mft(p_corono,
                                                         self._pupdiam,
                                                         planes='lyot_to_image',
                                                         center_on_pixel=True)
        self._indices_pup = np.where(self._spupil > 0.)

        self._coronagraph = PerfectCoronagraph(context, self._target.sources[0], 
                                               self._dim_image, self._dim_image, 
                                               self._wav_vec, self._wav_vec.size, 0)
        
        AA = np.rollaxis(np.array(self._AA), 0, self._wav_vec.size)
        BB = np.rollaxis(np.array(self._BB), 0, self._wav_vec.size)
        AA_c = np.rollaxis(np.array(self._AA_c), 0, self._wav_vec.size)
        BB_c = np.rollaxis(np.array(self._BB_c), 0, self._wav_vec.size)
        self._coronagraph.set_mft(AA, BB, self._norm0, scons.MftType.IMG)
        self._coronagraph.set_mft(AA_c, BB_c, self._norm0_c, scons.MftType.PSF)
        self._compute_normalization()

    def _compute_electric_field(self, input_opd, wavelength):
        """ Computes the electric field for the given OPD and wavelength (CPU version)

        Args:
            input_opd: (np.array): Input phase OPD in micron

            wavelength: (float): Wavelength in meter
        
        Return:
            electric_field: (np.ndarray[ndim=2, dtype=np.complex64]): Electric field
        """
        phase = input_opd * 2 * np.pi / wavelength
        electric_field = np.exp(1j * phase)
        return electric_field

    def _perfect_coro(self, input_opd, center_on_pixel=False, remove_coro=False):
        """ Applies a perfect coronagraph in the pupil plane
        and computes intensity in the focal plane. (CPU based)

        Args:
            input_opd: (np.array): Input OPD screens

            center_on_pixel: (bool, optional): If True, the image is centered on
                one pixel. Default = False, image centered between four pixels.

            remove_coro: (bool, optional): If True, returns the electric field
            without applying a perfect coronagraph. Default = False.

        Returns:
            iamge_intensity: (np.array): Intensity in focal plane, at all wavelengths
        """
        image_intensity = np.zeros((self._dim_image, self._dim_image))

        for i, wavelength in enumerate(self._wav_vec):
            electric_field = self._compute_electric_field(input_opd, wavelength)
            if remove_coro:
                electric_field_after_coro = electric_field * self._spupil
            else:
                electric_field_after_coro = (electric_field - np.mean(electric_field[self._indices_pup])) * self._spupil
            if center_on_pixel:
                image_electric_field = mft_multiplication(electric_field_after_coro,
                                                          self._AA_c[i],
                                                          self._BB_c[i],
                                                          self._norm0_c[i])
            else:
                image_electric_field = mft_multiplication(electric_field_after_coro,
                                                          self._AA[i],
                                                          self._BB[i],
                                                          self._norm0[i])
            image_intensity += np.abs(image_electric_field)**2

        return image_intensity

    def _compute_normalization_cpu(self):
        """ Computes the normalization factor of coronagraphic images (CPU based)
        """
        input_opd = np.zeros((self._pupdiam, self._pupdiam))
        self._norm_img = np.max(self._perfect_coro(input_opd, center_on_pixel=True, remove_coro=True))
        self._norm_psf = self._norm_img

    def _compute_normalization(self):
        """ Computes the normalization factor of coronagraphic images (CPU based)
        """
        self._target.reset_tar_phase(0)
        self.compute_psf(accumulate=False)
        self._norm_img = np.max(self.get_psf(expo_type="se"))
        self._norm_psf = self._norm_img

    def _compute_image(self, input_opd: np.array, *, accumulate: bool = True):
        """ Computes the coronographic image from input phase given as OPD (CPU based)

        Args:
            input_opd: (np.array): Input phase OPD

            accumulate: (bool, optional): If True (default), the computed image is added to the the long exposure image.
        """
        self._image_se = self._perfect_coro(input_opd) / self._norm_img
        self._psf_se = self._perfect_coro(input_opd, center_on_pixel=True, remove_coro=True) / self._norm_psf

        if(accumulate):
            self._cnt += 1
            self._image_le += self._image_se
            self._psf_le += self._psf_se

    def compute_image(self, comp_psf: bool=True, accumulate: bool = True):
        """ Compute the SE coronagraphic image, and accumulate it in the LE image

        Args:
            comp_psf: (bool, optionnal): If True (default), also compute the PSF SE & LE
            accumulate: (bool, optional): If True (default), the computed SE image is accumulated in 
                                            long exposure
        """
        self._coronagraph.compute_image(accumulate=accumulate)
        if comp_psf:
            self.compute_psf(accumulate=accumulate)

    def compute_psf(self, accumulate: bool = True):
        """ Compute the SE psf, and accumulate it in the LE image

        Args:
            accumulate: (bool, optional): If True (default), the computed SE psf is accumulated in 
                                            long exposure
        """
        self._coronagraph.compute_psf(accumulate=accumulate)