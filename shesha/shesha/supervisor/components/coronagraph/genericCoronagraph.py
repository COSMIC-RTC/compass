import numpy as np
import shesha.config as conf
from abc import ABC, abstractmethod
from shesha.supervisor.components.coronagraph.coronagraph_utils import compute_contrast

class GenericCoronagraph(ABC):

    def __init__(self, p_corono: conf.Param_corono, p_geom: conf.Param_geom):

        self._spupil = p_geom.get_spupil()
        self._pupdiam = self._spupil.shape[0]
        self._dim_image = p_corono._dim_image
        self._p_corono = p_corono

        self.image_se = np.zeros((self._dim_image, self._dim_image)) # Short exposure image
        self.image_le = np.zeros((self._dim_image, self._dim_image)) # Long exposure image
        self.psf_se = np.zeros((self._dim_image, self._dim_image))  # Short exposure PSF
        self.psf_le = np.zeros((self._dim_image, self._dim_image))  # Long exposure PSF
        self.cnt = 0

        self._aberrations = 0
        self._ab_buffer = None
        self._ab_cnt = 0

    @abstractmethod
    def compute_image(self):
        pass

    def get_image(self, *, expo_type="le"):
        """ Return the coronographic image

        Args:
            expo_type: (str, optional): If "le" (default), return the long exposure image. If "se", return the short exposure one.
        """
        if(expo_type == "le"):
            if (self.cnt):
                img = self.image_le / self.cnt
            else:
                img = self.image_le
        if(expo_type == "se"):
            img = self.image_se
        return img

    def get_psf(self, *, expo_type="le"):
        """ Return the PSF, i.e. the intensity without focal plane mask

        Args:
            expo_type: (str, optional): If "le" (default), return the long exposure psf. If "se", return the short exposure one.
        """
        if(expo_type == "le"):
            if (self.cnt):
                psf = self.psf_le / self.cnt
            else:
                psf = self.psf_le
        if(expo_type == "se"):
            psf = self.psf_se
        return psf

    def reset(self):  # add a reset_abb option ?
        """ Reset long exposure image and PSF
        """
        self.image_le *= 0
        self.psf_le *= 0
        self.cnt = 0

    def set_aberrations(self, aberrations):
        if len(aberrations.shape) == 2:
            ab_buffer = np.zeros(aberrations.shape + (1,))
            ab_buffer[:, :, 0] = aberrations[:, :]
        self._ab_buffer = ab_buffer
        self._ab_cnt = 0

    def reset_aberrations(self):
        self._ab_buffer = None
        self._aberrations = 0
        self._ab_cnt = 0

    def get_aberrations(self):
        return self._ab_buffer

    def _update_aberrations_buffer(self):
        if self._ab_buffer is not None:
            index = self._ab_cnt % self._ab_buffer.shape[2]
            self._aberrations = self._ab_buffer[:, :, index]
            self._ab_cnt += 1

    def get_contrast(self, *, expo_type='le', d_min=None, d_max=None, width=None, normalized_by_psf=True):
        """ Computes average, standard deviation, minimum and maximum of coronagraphic
        image intensity, over rings at several angular distances from the optical axis.

        A ring includes the pixels between the following angular distances :
        d_min + k * width - width / 2 and d_min + k * width + width / 2 (in lambda/D units)
        with k = 0, 1, 2... until d_min + k * width > d_max (excluded).

        Args:
            expo_type: (str, optional): If "le" (default), computes contrast on the long exposure image.
                                        If "se", it uses the short exposure one.

            d_min: (float, optional): Angular radius of the first ring in lambda/D unit.
                                      Default = width

            d_max: (float, optional): Maximum angular distances in lambda/D unit.
                                      Default includes the whole image.

            width: (float, optional): Width of one ring in lambda/D unit.
                                      Default = 1 [lambda/D]

            normalized_by_psf: (bool, optional): If True (default), the coronagraphic image
                                                 is normalized by the maximum of the PSF

        Returns:
            distances: (1D array): angular distances to the optical axis in lambda/D unit

            mean: (1D array): corresponding average intensities

            std: (1D array): corresponding standard deviations

            mini: (1D array): corresponding minimums

            maxi: (1D array): corresponding maximums
        """
        image_sampling = self._p_corono._image_sampling
        if width == None:
            width = image_sampling
        else:
            width = width * image_sampling
        if d_min == None:
            d_min = width
        else:
            d_min = d_min * image_sampling
        if d_max == None:
            d_max = self._dim_image / 2 - width / 2
        else:
            d_max = d_max * image_sampling

        center = self._dim_image / 2 - (1 / 2)
        image = self.get_image(expo_type=expo_type)
        if normalized_by_psf:
            image = image / np.max(self.get_psf(expo_type=expo_type))

        distances, mean, std, mini, maxi = compute_contrast(image, center, d_min, d_max, width)
        angular_distances = distances / image_sampling
        return angular_distances, mean, std, mini, maxi
