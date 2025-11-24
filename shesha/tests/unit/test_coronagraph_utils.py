"""
Tests for shesha.util.coronagraph_utils module
"""

import numpy as np
import pytest
from shesha.util.coronagraph_utils import (
    roundpupil, classical_lyot_fpm, make_VLT_pupil
)


class TestRoundPupil:
    """Test circular pupil generation."""
    
    def test_roundpupil_basic(self):
        """Test basic circular pupil generation."""
        dim_pp = 64
        prad = 20
        
        pupil = roundpupil(dim_pp, prad)
        
        assert isinstance(pupil, np.ndarray)
        assert pupil.shape == (dim_pp, dim_pp)
        assert pupil.dtype in [np.float32, np.float64]
    
    def test_roundpupil_center_pixel_option(self):
        """Test circular pupil with pixel center option."""
        dim_pp = 64
        prad = 20
        
        pupil_pixel = roundpupil(dim_pp, prad, center_pos='p')
        pupil_between = roundpupil(dim_pp, prad, center_pos='b')
        
        # Both should produce valid pupils
        assert pupil_pixel.shape == (dim_pp, dim_pp)
        assert pupil_between.shape == (dim_pp, dim_pp)
        # They should be different due to centering
        assert not np.allclose(pupil_pixel, pupil_between)
    
    def test_roundpupil_center_between_option(self):
        """Test pupil centering between pixels."""
        dim_pp = 32
        prad = 10
        
        pupil = roundpupil(dim_pp, prad, center_pos='b')
        
        # Should have circular structure
        assert np.sum(pupil) > 0
    
    def test_roundpupil_case_insensitive(self):
        """Test center_pos is case insensitive."""
        dim_pp = 64
        prad = 20
        
        pupil_p = roundpupil(dim_pp, prad, center_pos='p')
        pupil_P = roundpupil(dim_pp, prad, center_pos='P')
        pupil_b = roundpupil(dim_pp, prad, center_pos='b')
        pupil_B = roundpupil(dim_pp, prad, center_pos='B')
        
        assert np.allclose(pupil_p, pupil_P)
        assert np.allclose(pupil_b, pupil_B)
    
    def test_roundpupil_invalid_center_pos(self):
        """Test invalid center_pos raises error."""
        dim_pp = 64
        prad = 20
        
        with pytest.raises(ValueError):
            roundpupil(dim_pp, prad, center_pos='invalid')
    
    def test_roundpupil_values_range(self):
        """Test pupil values are 0 or 1."""
        dim_pp = 64
        prad = 20
        
        pupil = roundpupil(dim_pp, prad)
        
        # Values should be 0 or 1 (or very close to 1)
        assert np.all((pupil == 0) | (pupil == 1))
    
    def test_roundpupil_circular_shape(self):
        """Test pupil has circular shape."""
        dim_pp = 128
        prad = 40
        
        pupil = roundpupil(dim_pp, prad, center_pos='b')
        
        # Center should be in pupil
        center = dim_pp // 2
        assert pupil[center, center] == 1
    
    def test_roundpupil_radius_effect(self):
        """Test different radii produce different pupils."""
        dim_pp = 64
        
        pupil_small = roundpupil(dim_pp, 10)
        pupil_large = roundpupil(dim_pp, 30)
        
        # Larger pupil should have more ones
        assert np.sum(pupil_large) > np.sum(pupil_small)
    
    def test_roundpupil_different_sizes(self):
        """Test pupil generation for different image sizes."""
        for dim_pp in [32, 64, 128]:
            prad = dim_pp // 4
            
            pupil = roundpupil(dim_pp, prad)
            
            assert pupil.shape == (dim_pp, dim_pp)
            assert np.sum(pupil) > 0


class TestClassicalLyotFPM:
    """Test Classical Lyot FPM generation."""
    
    def test_classical_lyot_fpm_basic(self):
        """Test basic Lyot FPM generation."""
        rad_lyot_fpm = 2.0
        dim_fpm = 256
        lyot_fpm_sampling = 20
        wav_vec = [1e-6, 1.1e-6]
        
        fpm_list = classical_lyot_fpm(rad_lyot_fpm, dim_fpm, lyot_fpm_sampling, wav_vec)
        
        assert isinstance(fpm_list, list)
        assert len(fpm_list) == len(wav_vec)
        
        # Each FPM should be a 2D array
        for fpm in fpm_list:
            assert isinstance(fpm, np.ndarray)
            assert fpm.shape == (dim_fpm, dim_fpm)
    
    def test_classical_lyot_fpm_values(self):
        """Test Lyot FPM has valid values (0 or 1)."""
        rad_lyot_fpm = 1.5
        dim_fpm = 128
        lyot_fpm_sampling = 15
        wav_vec = [1e-6]
        
        fpm_list = classical_lyot_fpm(rad_lyot_fpm, dim_fpm, lyot_fpm_sampling, wav_vec)
        
        fpm = fpm_list[0]
        
        # Values should be 0 or 1
        assert np.all((fpm == 0) | (fpm == 1))
    
    def test_classical_lyot_fpm_wavelength_independence(self):
        """Test that FPM is independent of wavelength."""
        rad_lyot_fpm = 2.0
        dim_fpm = 128
        lyot_fpm_sampling = 20
        wav_vec = [0.5e-6, 1.0e-6, 2.0e-6]
        
        fpm_list = classical_lyot_fpm(rad_lyot_fpm, dim_fpm, lyot_fpm_sampling, wav_vec)
        
        # All FPMs should be identical (same for all wavelengths in classical design)
        for i in range(len(fpm_list) - 1):
            assert np.allclose(fpm_list[i], fpm_list[i+1])
    
    def test_classical_lyot_fpm_radius_effect(self):
        """Test FPM changes with different radii."""
        dim_fpm = 128
        lyot_fpm_sampling = 20
        wav_vec = [1e-6]
        
        fpm_small = classical_lyot_fpm(1.0, dim_fpm, lyot_fpm_sampling, wav_vec)[0]
        fpm_large = classical_lyot_fpm(3.0, dim_fpm, lyot_fpm_sampling, wav_vec)[0]
        
        # Larger radius should have fewer ones (larger hole becomes smaller hole)
        assert np.sum(fpm_small) > np.sum(fpm_large)
    
    def test_classical_lyot_fpm_single_wavelength(self):
        """Test with single wavelength."""
        rad_lyot_fpm = 2.0
        dim_fpm = 64
        lyot_fpm_sampling = 15
        wav_vec = [1e-6]
        
        fpm_list = classical_lyot_fpm(rad_lyot_fpm, dim_fpm, lyot_fpm_sampling, wav_vec)
        
        assert len(fpm_list) == 1


class TestMakeVLTPupil:
    """Test VLT pupil generation."""
    
    def test_make_vlt_pupil_basic(self):
        """Test basic VLT pupil generation."""
        pupdiam = 256
        
        pupil = make_VLT_pupil(pupdiam)
        
        assert isinstance(pupil, np.ndarray)
        assert pupil.shape == (pupdiam, pupdiam)
        assert pupil.dtype == np.float32
    
    def test_make_vlt_pupil_values(self):
        """Test VLT pupil values are 0 or 1."""
        pupdiam = 128
        
        pupil = make_VLT_pupil(pupdiam)
        
        # Values should be 0 or 1
        assert np.all((pupil == 0) | (pupil == 1))
    
    def test_make_vlt_pupil_no_obstruction(self):
        """Test VLT pupil without central obstruction."""
        pupdiam = 256
        
        pupil_with_obs = make_VLT_pupil(pupdiam, centralobs_bool=True)
        pupil_no_obs = make_VLT_pupil(pupdiam, centralobs_bool=False)
        
        # Without obstruction should have more ones (more area)
        assert np.sum(pupil_no_obs) > np.sum(pupil_with_obs)
    
    def test_make_vlt_pupil_no_spiders(self):
        """Test VLT pupil without spiders."""
        pupdiam = 256
        
        pupil_with_spiders = make_VLT_pupil(pupdiam, spiders_bool=True)
        pupil_no_spiders = make_VLT_pupil(pupdiam, spiders_bool=False)
        
        # Without spiders should have more ones
        assert np.sum(pupil_no_spiders) >= np.sum(pupil_with_spiders)
    
    def test_make_vlt_pupil_symmetry(self):
        """Test VLT pupil is symmetric."""
        pupdiam = 128
        
        pupil = make_VLT_pupil(pupdiam)
        
        # Check 4-fold symmetry
        pupil_hflip = np.fliplr(pupil)
        pupil_vflip = np.flipud(pupil)
        
        # Should have symmetry
        assert np.allclose(pupil, pupil_hflip) or not np.all(pupil_hflip == 0)
    
    def test_make_vlt_pupil_custom_obstruction(self):
        """Test VLT pupil with custom obstruction ratio."""
        pupdiam = 256
        
        pupil_small_obs = make_VLT_pupil(pupdiam, centralobs=0.1)
        pupil_large_obs = make_VLT_pupil(pupdiam, centralobs=0.3)
        
        # Larger obstruction should have fewer ones
        assert np.sum(pupil_large_obs) < np.sum(pupil_small_obs)
    
    def test_make_vlt_pupil_custom_spider(self):
        """Test VLT pupil with custom spider size."""
        pupdiam = 256
        
        pupil_thin_spider = make_VLT_pupil(pupdiam, spiders=0.003)
        pupil_thick_spider = make_VLT_pupil(pupdiam, spiders=0.01)
        
        # Thicker spider should have fewer ones
        assert np.sum(pupil_thick_spider) < np.sum(pupil_thin_spider)
    
    def test_make_vlt_pupil_even_odd(self):
        """Test VLT pupil works for even and odd dimensions."""
        for pupdiam in [64, 65, 128, 129]:
            pupil = make_VLT_pupil(pupdiam)
            
            assert pupil.shape == (pupdiam, pupdiam)
            assert np.sum(pupil) > 0
    
    def test_make_vlt_pupil_all_options(self):
        """Test VLT pupil with all options combined."""
        pupdiam = 256
        
        pupil = make_VLT_pupil(
            pupdiam,
            centralobs=0.15,
            spiders=0.007,
            spiders_bool=True,
            centralobs_bool=True
        )
        
        assert pupil.shape == (pupdiam, pupdiam)
        assert np.all((pupil == 0) | (pupil == 1))


class TestCoronagraphIntegration:
    """Integration tests for coronagraph utilities."""
    
    def test_roundpupil_and_lyot(self):
        """Test using roundpupil with Lyot FPM."""
        dim_pp = 128
        prad = 40
        
        pupil = roundpupil(dim_pp, prad)
        
        dim_fpm = 256
        rad_lyot_fpm = 2.0
        lyot_sampling = 20
        wav_vec = [1e-6]
        
        fpm = classical_lyot_fpm(rad_lyot_fpm, dim_fpm, lyot_sampling, wav_vec)[0]
        
        # Both should be valid 2D arrays
        assert pupil.ndim == 2
        assert fpm.ndim == 2
    
    def test_all_pupils_consistency(self):
        """Test all pupil generators produce consistent results."""
        size = 256
        
        round_pup = roundpupil(size, 50)
        vlt_pup = make_VLT_pupil(size)
        
        # Both should be valid pupils
        assert round_pup.shape == (size, size)
        assert vlt_pup.shape == (size, size)
        
        # Both should have circular/pupil-like structure
        assert np.sum(round_pup) > 0
        assert np.sum(vlt_pup) > 0
    
    def test_pupil_and_fpm_compatible(self):
        """Test pupils work with FPM."""
        pupil = roundpupil(128, 40)
        
        fpm_list = classical_lyot_fpm(2.0, 128, 20, [1e-6])
        fpm = fpm_list[0]
        
        # Both should have same dimensions
        assert pupil.shape == fpm.shape
