import numpy as np
import shesha.config as conf
import shesha.constants as scons
from shesha.supervisor.components.coronagraph.genericCoronagraph import GenericCoronagraph
from shesha.init.coronagraph_init import init_coronagraph, init_mft, mft_multiplication
from shesha.supervisor.components.targetCompass import TargetCompass
from sutraWrap import StellarCoronagraph
from carmaWrap import context
class StellarCoronagraph(GenericCoronagraph):
    """ Class supervising stellar coronagraph component

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

        _AA_apod_to_fpm: (np.ndarray[ndim=3, dtype=np.complex64]): MFT matrix for focal plane

        _BB_apod_to_fpm: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for focal plane

        _norm0_apod_to_fpm: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for focal plane

        _AA_fpm_to_lyot: (np.ndarray[ndim=3, dtype=np.complex64]): MFT matrix for lyot plane

        _BB_fpm_to_lyot: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for lyot plane

        _norm0_fpm_to_lyot: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for lyot plane

        _AA_lyot_to_image_c: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for image computation (centered on pixel)

        _BB_lyot_to_image_c: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for psf computation (centered on pixel)
        
        _norm0_lyot_to_image_c: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for psf computation (centered on pixel)

        _AA_lyot_to_image: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for image computation

        _BB_lyot_to_image: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for psf computation
        
        _norm0_lyot_to_image: (np.ndarray[ndim=2, dtype=np.complex64]): MFT matrix for psf computation

        _indices_pup: (tuple): Tuple of ndarray containing X and Y indices of illuminated 
                                pixels in the pupil
    """
    def __init__(self, context: context, targetCompass: TargetCompass, p_corono: conf.Param_corono, 
                 p_geom: conf.Param_geom):
        """ Initialize a stellar coronagraph instance

        Args:
            context: (CarmaWrap.context): GPU context

            targetCompass: (TargetCompass): Compass Target used as input for the coronagraph

            p_corono: (Param_corono): Coronagraph parameters

            p_geom: (Param_geom): Compass geometry parameters
        """

        init_coronagraph(p_corono, p_geom.pupdiam)
        GenericCoronagraph.__init__(self, p_corono, p_geom, targetCompass)
        self._wav_vec = p_corono._wav_vec

        self._AA_apod_to_fpm, self._BB_apod_to_fpm, self._norm0_apod_to_fpm = init_mft(self._p_corono,
                                                                                       self._pupdiam,
                                                                                       planes='apod_to_fpm')
        self._AA_fpm_to_lyot, self._BB_fpm_to_lyot, self._norm0_fpm_to_lyot = init_mft(self._p_corono,
                                                                                       self._pupdiam,
                                                                                       planes='fpm_to_lyot')
        self._AA_lyot_to_image, self._BB_lyot_to_image, self._norm0_lyot_to_image = init_mft(self._p_corono,
                                                                                             self._pupdiam,
                                                                                             planes='lyot_to_image')
        self._AA_lyot_to_image_c, self._BB_lyot_to_image_c, self._norm0_lyot_to_image_c = init_mft(self._p_corono,
                                                                                                   self._pupdiam,
                                                                                                   planes='lyot_to_image',
                                                                                                   center_on_pixel=True)
        self._coronagraph = StellarCoronagraph(context, self._target.sources[0], 
                                               self._dim_image, self._dim_image, 
                                               self._p_corono._dim_fpm, self._p_corono._dim_fpm, 
                                               self._wav_vec, self._wav_vec.size, 
                                               self._p_corono._babinet_trick, 0)
        
        AA = np.rollaxis(np.array(self._AA_lyot_to_image), 0, self._wav_vec.size)
        BB = np.rollaxis(np.array(self._BB_lyot_to_image), 0, self._wav_vec.size)
        AA_c = np.rollaxis(np.array(self._AA_lyot_to_image_c), 0, self._wav_vec.size)
        BB_c = np.rollaxis(np.array(self._BB_lyot_to_image_c), 0, self._wav_vec.size)
        AA_fpm = np.rollaxis(np.array(self._AA_apod_to_fpm), 0, self._wav_vec.size)
        BB_fpm = np.rollaxis(np.array(self._BB_apod_to_fpm), 0, self._wav_vec.size)
        AA_lyot = np.rollaxis(np.array(self._AA_fpm_to_lyot), 0, self._wav_vec.size)
        BB_lyot = np.rollaxis(np.array(self._BB_fpm_to_lyot), 0, self._wav_vec.size)
        
        self._coronagraph.set_mft(AA, BB, self._norm0_lyot_to_image, scons.MftType.IMG)
        self._coronagraph.set_mft(AA_c, BB_c, self._norm0_lyot_to_image_c, scons.MftType.PSF)
        self._coronagraph.set_mft(AA_fpm, BB_fpm, self._norm0_apod_to_fpm, scons.MftType.FPM)
        self._coronagraph.set_mft(AA_lyot, BB_lyot, self._norm0_fpm_to_lyot, scons.MftType.LYOT)

        self._coronagraph.set_apodizer(self._p_corono._apodizer)
        self._coronagraph.set_lyot_stop(self._p_corono._lyot_stop)
        fpm = np.rollaxis(np.array(self._p_corono._focal_plane_mask), 0, self._wav_vec.size)
        if self._p_corono._babinet_trick:
            fpm = 1. - fpm
        self._coronagraph.set_focal_plane_mask(fpm)

        self._compute_normalization()

    def _compute_electric_field(self, input_opd, wavelength):
        """
        Args:
            input_opd: (np.array): Input phase OPD in micron

            wavelength: (float): Wavelength in meter
        """
        phase = input_opd * 2 * np.pi / wavelength
        electric_field = np.exp(1j * phase)
        return electric_field

    def _propagate_through_coro(self, input_opd, center_on_pixel=False, no_fpm=False):
        """ Propagate the electric field through the coronagraph
        and compute the intensity in the imaging plane.
        """
        image_intensity = np.zeros((self._dim_image, self._dim_image))

        for i, wavelength in enumerate(self._wav_vec):
            EF_before_apod = self._compute_electric_field(input_opd, wavelength) * self._spupil
            EF_after_apod = EF_before_apod * self._p_corono._apodizer

            EF_before_fpm = mft_multiplication(EF_after_apod,
                                               self._AA_apod_to_fpm[i],
                                               self._BB_apod_to_fpm[i],
                                               self._norm0_apod_to_fpm[i])

            if no_fpm:
                fpm = 1.
            else:
                fpm = self._p_corono._focal_plane_mask[i]

            if self._p_corono._babinet_trick:
                EF_after_fpm_babinet = EF_before_fpm * (1. - fpm)  # Babinet's trick
                EF_before_lyot_babinet = mft_multiplication(EF_after_fpm_babinet,
                                                            self._AA_fpm_to_lyot[i],
                                                            self._BB_fpm_to_lyot[i],
                                                            self._norm0_fpm_to_lyot[i])
                EF_before_lyot = EF_after_apod - EF_before_lyot_babinet
            else:
                EF_after_fpm = EF_before_fpm * fpm
                EF_before_lyot = mft_multiplication(EF_after_fpm,
                                                    self._AA_fpm_to_lyot[i],
                                                    self._BB_fpm_to_lyot[i],
                                                    self._norm0_fpm_to_lyot[i])
            EF_after_lyot = EF_before_lyot * self._p_corono._lyot_stop

            if center_on_pixel:
                image_electric_field = mft_multiplication(EF_after_lyot,
                                                          self._AA_lyot_to_image_c[i],
                                                          self._BB_lyot_to_image_c[i],
                                                          self._norm0_lyot_to_image_c[i])
            else:
                image_electric_field = mft_multiplication(EF_after_lyot,
                                                          self._AA_lyot_to_image[i],
                                                          self._BB_lyot_to_image[i],
                                                          self._norm0_lyot_to_image[i])

            image_intensity += np.abs(image_electric_field)**2
        return image_intensity

    def _compute_psf(self, input_opd, center_on_pixel=True):
        """ Compute |TF(exp(i*phi) * pup * apod * lyot_stop)|**2
        """
        psf_intensity = np.zeros((self._dim_image, self._dim_image))
        for i, wavelength in enumerate(self._wav_vec):
            EF_before_apod = self._compute_electric_field(input_opd, wavelength) * self._spupil
            EF_after_lyot = EF_before_apod * self._p_corono._apodizer * self._p_corono._lyot_stop
            if center_on_pixel:
                psf_electric_field = mft_multiplication(EF_after_lyot,
                                                        self._AA_lyot_to_image_c[i],
                                                        self._BB_lyot_to_image_c[i],
                                                        self._norm0_lyot_to_image_c[i])
            else:
                psf_electric_field = mft_multiplication(EF_after_lyot,
                                                        self._AA_lyot_to_image[i],
                                                        self._BB_lyot_to_image[i],
                                                        self._norm0_lyot_to_image[i])
            psf_intensity += np.abs(psf_electric_field)**2
        return psf_intensity

    def _compute_normalization_cpu(self):
        """ Compute the normalization factor of coronagraphic images
        """
        input_opd = np.zeros((self._pupdiam, self._pupdiam))

        self._norm_img = np.max(self._propagate_through_coro(input_opd,
                                                               no_fpm=True,
                                                               center_on_pixel=True))
        self._norm_psf = np.max(self._compute_psf(input_opd))

    def _compute_normalization(self):
        """ Compute the normalization factor of coronagraphic images
        """
        self._target.reset_tar_phase(0)
        self._coronagraph.compute_image_normalization()
        self._norm_img = np.max(self.get_image(expo_type=scons.ExposureType.SE))
        self.compute_psf(accumulate=False)
        self._norm_psf = np.max(self.get_psf(expo_type=scons.ExposureType.SE))

    def _compute_image(self, *, accumulate: bool = True):
        """ Computes the coronographic image from input phase given as OPD

        Args:
            input_opd: (np.array): Input phase OPD

            accumulate: (bool, optional): If True (default), the computed image is added to the the long exposure image.
        """
        self._image_se = self._propagate_through_coro(self._target.get_tar_phase(0)) / self._norm_img
        self._psf_se = self._compute_psf(self._target.get_tar_phase(0)) / self._norm_psf

        if(accumulate):
            self._cnt += 1
            self._image_le += self._image_se
            self._psf_le += self._psf_se
