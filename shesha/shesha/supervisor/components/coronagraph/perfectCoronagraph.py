import numpy as np
import shesha.config as conf
from shesha.supervisor.components.coronagraph.genericCoronagraph import GenericCoronagraph
from shesha.supervisor.components.coronagraph.coronagraph_init import init_coronagraph, init_mft, mft_multiplication


class PerfectCoronagraph(GenericCoronagraph):
    """ Class supervising perfect coronagraph component
    """
    def __init__(self, p_corono: conf.Param_corono, p_geom: conf.Param_geom):

        init_coronagraph(p_corono, p_geom.pupdiam)
        GenericCoronagraph.__init__(self, p_corono, p_geom)
        self._wav_vec = p_corono._wav_vec

        self._AA, self._BB, self._norm0 = init_mft(p_corono,
                                                   self._pupdiam,
                                                   planes='lyot_to_image')
        self._AA_c, self._BB_c, self._norm0_c = init_mft(p_corono,
                                                         self._pupdiam,
                                                         planes='lyot_to_image',
                                                         center_on_pixel=True)
        self._indices_pup = np.where(self._spupil > 0.)

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

    def _perfect_coro(self, input_opd, center_on_pixel=False, remove_coro=False):
        """ Applies a perfect coronagraph in the pupil plane
        and computes intensity in the focal plane.

        Args:
            input_opd: (np.array): Input OPD screens

            remove_coro: (bool, optional): If True, returns the electric field
            without applying a perfect coronagraph. Default = False.

        Returns:
            electric_field_after_coro: (np.array): The electric field in a pupil plane
                after applying a perfect coronagraph.
        """
        image_intensity = np.zeros((self._dim_image, self._dim_image))

        for i, wavelength in enumerate(self._wav_vec):
            electric_field = self._compute_electric_field(input_opd, wavelength)
            if remove_coro:
                electric_field_after_coro = electric_field
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

    def _compute_normalization(self):
        """ Computes the normalization factor of coronagraphic images
        """
        input_opd = np.zeros((self._pupdiam, self._pupdiam))
        self._norm_image = np.max(self._perfect_coro(input_opd, center_on_pixel=True))
        self._norm_psf = np.max(self._perfect_coro(input_opd, center_on_pixel=True, remove_coro=True))

    def compute_image(self, input_opd: np.array, *, accumulate: bool = True):
        """ Computes the coronographic image from input phase given as OPD

        Args:
            input_opd: (np.array): Input phase OPD

            accumulate: (bool, optional): If True (default), the computed image is added to the the long exposure image.
        """
        self._update_aberrations_buffer()

        self.image_se = self._perfect_coro(input_opd) / self._norm_image
        self.psf_se = self._perfect_coro(input_opd, center_on_pixel=True, remove_coro=True) / self._norm_psf

        if(accumulate):
            self.cnt += 1
            self.image_le += self.image_se
            self.psf_le += self.psf_se
