import numpy as np
import shesha.config as conf
from shesha.supervisor.components.coronagraph.genericCoronagraph import GenericCoronagraph
from shesha.supervisor.components.coronagraph.coronagraph_init import init_coronagraph, init_mft, mft_multiplication


class ClassicalCoronagraph(GenericCoronagraph):
    """ Class supervising coronagraph component
    """
    def __init__(self, p_corono: conf.Param_corono, p_geom: conf.Param_geom):

        init_coronagraph(p_corono, p_geom.pupdiam)
        GenericCoronagraph.__init__(self, p_corono, p_geom)
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
        self._compute_normalization()

    def _compute_electric_field(self, input_opd, wavelength):
        """
        Args:
            input_opd: (np.array): Input phase OPD in micron

            wavelength: (float): Wavelength in meter
        """
        phase = (input_opd + self._aberrations) * 1e-6 * 2 * np.pi / wavelength
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
                                                    self._norm0_fpm_to_lyot)
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

    def _compute_normalization(self):
        """ Compute the normalization factor of coronagraphic images
        """
        input_opd = np.zeros((self._pupdiam, self._pupdiam))

        self._norm_image = np.max(self._propagate_through_coro(input_opd,
                                                               no_fpm=True,
                                                               center_on_pixel=True))
        self._norm_psf = np.max(self._compute_psf(input_opd))

    def compute_image(self, input_opd: np.array, *, accumulate: bool = True):
        """ Computes the coronographic image from input phase given as OPD

        Args:
            input_opd: (np.array): Input phase OPD

            accumulate: (bool, optional): If True (default), the computed image is added to the the long exposure image.
        """
        self._update_aberrations_buffer()

        self.image_se = self._propagate_through_coro(input_opd) / self._norm_image
        self.psf_se = self._compute_psf(input_opd) / self._norm_psf

        if(accumulate):
            self.cnt += 1
            self.image_le += self.image_se
            self.psf_le += self.psf_se



## TODO / questions:
# introduce phase and amplitude aberrations. NCPA dans compass ?
# Johan : sphere image size ? why 256 pixels ?

## further work :
# padded pupil options ?
# pupangle handling ?
# grey pupil ?
# hardcoding sphere parameters in pupil functions
# rewrite correctly VLT pupil (VLT pupil 'range' issue, maybe use meshgrid instead)