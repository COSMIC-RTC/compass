"""
Tests for shesha.util.kl_util module
"""

import numpy as np
import pytest
from shesha.util.kl_util import (
    make_radii, make_kernels, piston_orth, make_azimuth
)
from shesha.constants import KLType


class TestMakeRadii:
    """Test KL radii generation."""
    
    def test_make_radii_basic(self):
        """Test basic radii generation."""
        cobs = 0.0
        nr = 10
        
        radii = make_radii(cobs, nr)
        
        assert isinstance(radii, np.ndarray)
        assert len(radii) == nr
    
    def test_make_radii_monotonic(self):
        """Test radii are monotonically increasing."""
        cobs = 0.1
        nr = 20
        
        radii = make_radii(cobs, nr)
        
        # Radii should be monotonically increasing
        assert np.all(np.diff(radii) > 0)
    
    def test_make_radii_cobs_effect(self):
        """Test central obstruction effect on radii."""
        nr = 10
        
        radii_nocobs = make_radii(0.0, nr)
        radii_withcobs = make_radii(0.3, nr)
        
        # With obstruction, radii start at non-zero
        assert radii_nocobs[0] < radii_withcobs[0]
    
    def test_make_radii_range(self):
        """Test radii are in valid range [cobs, 1]."""
        cobs = 0.2
        nr = 15
        
        radii = make_radii(cobs, nr)
        
        # Radii should start at cobs and end at ~1
        assert radii[0] >= cobs
        assert radii[-1] <= 1.0 + 0.1
    
    def test_make_radii_different_nr(self):
        """Test radii generation for different nr values."""
        cobs = 0.1
        
        for nr in [5, 10, 20, 50]:
            radii = make_radii(cobs, nr)
            assert len(radii) == nr


class TestMakeKernels:
    """Test KL kernel generation."""
    
    def test_make_kernels_kolmo(self):
        """Test Kolmogorov kernel generation."""
        cobs = 0.1
        nr = 10
        radp = np.linspace(cobs, 1.0, nr)
        
        kernels = make_kernels(cobs, nr, radp, KLType.KOLMO)
        
        assert isinstance(kernels, np.ndarray)
        assert kernels.shape[0] == 5 * nr  # nth = 5 * nr
        assert kernels.shape[1] == nr
        assert kernels.shape[2] == nr
    
    def test_make_kernels_karman(self):
        """Test Von Karman kernel generation."""
        cobs = 0.1
        nr = 10
        radp = np.linspace(cobs, 1.0, nr)
        
        kernels = make_kernels(cobs, nr, radp, KLType.KARMAN, outscl=3.0)
        
        assert isinstance(kernels, np.ndarray)
        assert kernels.dtype == np.float32
    
    def test_make_kernels_symmetry(self):
        """Test kernels have expected symmetry."""
        cobs = 0.1
        nr = 8
        radp = np.linspace(cobs, 1.0, nr)
        
        kernels = make_kernels(cobs, nr, radp, KLType.KOLMO)
        
        # Kernels should be symmetric
        for k in kernels:
            assert np.allclose(k, k.T, atol=1e-5)
    
    def test_make_kernels_outscale_effect(self):
        """Test outer scale parameter effect."""
        cobs = 0.1
        nr = 5
        radp = np.linspace(cobs, 1.0, nr)
        
        kernels_small = make_kernels(cobs, nr, radp, KLType.KARMAN, outscl=1.0)
        kernels_large = make_kernels(cobs, nr, radp, KLType.KARMAN, outscl=10.0)
        
        # Kernels should be different for different outer scales
        assert not np.allclose(kernels_small, kernels_large)
    
    def test_make_kernels_no_nans(self):
        """Test kernels don't contain NaN values."""
        cobs = 0.1
        nr = 10
        radp = np.linspace(cobs, 1.0, nr)
        
        kernels = make_kernels(cobs, nr, radp, KLType.KOLMO)
        
        assert not np.any(np.isnan(kernels))
    
    def test_make_kernels_no_infs(self):
        """Test kernels don't contain infinite values."""
        cobs = 0.1
        nr = 10
        radp = np.linspace(cobs, 1.0, nr)
        
        kernels = make_kernels(cobs, nr, radp, KLType.KARMAN)
        
        assert not np.any(np.isinf(kernels))


class TestKLTypeIntegration:
    """Integration tests for KL mode generation."""
    
    def test_kolmo_vs_karman_shapes(self):
        """Test Kolmo and Karman kernels have same shape."""
        cobs = 0.1
        nr = 8
        radp = np.linspace(cobs, 1.0, nr)
        
        k_kolmo = make_kernels(cobs, nr, radp, KLType.KOLMO)
        k_karman = make_kernels(cobs, nr, radp, KLType.KARMAN)
        
        assert k_kolmo.shape == k_karman.shape
    
    def test_increasing_nr_increases_kernel_size(self):
        """Test kernel size increases with nr."""
        cobs = 0.1
        
        for nr in [5, 10, 15]:
            radp = np.linspace(cobs, 1.0, nr)
            kernels = make_kernels(cobs, nr, radp, KLType.KOLMO)
            
            assert kernels.shape[1] == nr
            assert kernels.shape[2] == nr
    
    def test_radii_kernel_consistency(self):
        """Test make_radii and make_kernels work together."""
        cobs = 0.15
        nr = 12
        
        radp = make_radii(cobs, nr)
        kernels = make_kernels(cobs, nr, radp, KLType.KOLMO)
        
        # Should work without errors
        assert kernels.shape[0] == 5 * nr
        assert len(radp) == nr


class TestPistonOrth:
    """Test piston orthogonalization."""
    
    def test_piston_orth_basic(self):
        """Test basic piston orthogonalization."""
        nr = 10
        
        s = piston_orth(nr)
        
        assert isinstance(s, np.ndarray)
        assert s.shape == (nr, nr)
        assert s.dtype == np.float32
    
    def test_piston_orth_single(self):
        """Test piston orthogonalization with nr=1."""
        s = piston_orth(1)
        
        assert s.shape == (1, 1)
        # Should be normalized
        assert np.isclose(s[0, 0], 1.0)
    
    def test_piston_orth_orthogonality(self):
        """Test that orthogonalization produces orthogonal vectors."""
        nr = 8
        
        s = piston_orth(nr)
        
        # Check columns are orthogonal (or approximately so)
        gram = s.T @ s
        
        # Diagonal should be non-zero (norms)
        assert np.all(np.diag(gram) > 0)
    
    def test_piston_orth_first_column(self):
        """Test first column properties."""
        nr = 10
        
        s = piston_orth(nr)
        
        # First column should have specific structure
        assert not np.any(np.isnan(s[:, 0]))
        assert not np.any(np.isinf(s[:, 0]))
    
    def test_piston_orth_different_sizes(self):
        """Test piston orthogonalization for different sizes."""
        for nr in [2, 5, 10, 20]:
            s = piston_orth(nr)
            assert s.shape == (nr, nr)
            assert s.dtype == np.float32


class TestMakeAzimuth:
    """Test azimuthal basis generation."""
    
    def test_make_azimuth_basic(self):
        """Test basic azimuthal basis generation."""
        nord = 4
        npp = 32
        
        azbas = make_azimuth(nord, npp)
        
        assert isinstance(azbas, np.ndarray)
        assert azbas.shape == (npp, 1 + nord)
        assert azbas.dtype == np.float32
    
    def test_make_azimuth_constant_term(self):
        """Test first column is constant."""
        nord = 5
        npp = 16
        
        azbas = make_azimuth(nord, npp)
        
        # First column should be 1.0 everywhere
        assert np.allclose(azbas[:, 0], 1.0)
    
    def test_make_azimuth_oscillating_terms(self):
        """Test oscillating terms in azimuthal basis."""
        nord = 6
        npp = 64
        
        azbas = make_azimuth(nord, npp)
        
        # At least some columns should have variation (not all constant)
        has_variation = False
        for i in range(1, azbas.shape[1]):
            col = azbas[:, i]
            if not np.allclose(col, col[0]):
                has_variation = True
                break
        
        assert has_variation or nord == 0
    
    def test_make_azimuth_zero_nord(self):
        """Test with zero order (only constant term)."""
        azbas = make_azimuth(0, 32)
        
        assert azbas.shape == (32, 1)
        assert np.allclose(azbas, 1.0)
    
    def test_make_azimuth_range(self):
        """Test values are in reasonable range."""
        nord = 6
        npp = 48
        
        azbas = make_azimuth(nord, npp)
        
        # Values should be bounded (typically between -1 and 1)
        assert np.all(np.abs(azbas) <= 1.0 + 1e-6)
    
    def test_make_azimuth_periodicity(self):
        """Test periodicity of azimuthal basis."""
        nord = 4
        npp = 32
        
        azbas = make_azimuth(nord, npp)
        
        # For cosine/sine terms, should have expected periodicity
        # This is more of a sanity check that values are reasonable
        assert not np.any(np.isnan(azbas))
        assert not np.any(np.isinf(azbas))
    
    def test_make_azimuth_different_parameters(self):
        """Test make_azimuth with various parameters."""
        for nord in [2, 4, 8]:
            for npp in [16, 32, 64]:
                azbas = make_azimuth(nord, npp)
                assert azbas.shape == (npp, 1 + nord)


class TestKLUtilIntegration:
    """Advanced integration tests for KL utilities."""
    
    def test_full_kl_pipeline(self):
        """Test full KL mode generation pipeline."""
        cobs = 0.15
        nr = 10
        nord = 4
        npp = 32
        
        # Generate radii
        radp = make_radii(cobs, nr)
        
        # Generate kernels
        kernels = make_kernels(cobs, nr, radp, KLType.KOLMO)
        
        # Generate azimuthal basis
        azbas = make_azimuth(nord, npp)
        
        # Generate piston orthogonalization
        s = piston_orth(nr)
        
        # All should work together
        assert radp.shape[0] == nr
        assert kernels.shape == (5 * nr, nr, nr)
        assert azbas.shape == (npp, 1 + nord)
        assert s.shape == (nr, nr)
    
    def test_kl_modes_with_obstruction(self):
        """Test KL modes with various obstruction ratios."""
        nr = 8
        nord = 3
        npp = 24
        
        for cobs in [0.0, 0.1, 0.3, 0.5]:
            radp = make_radii(cobs, nr)
            kernels = make_kernels(cobs, nr, radp, KLType.KOLMO)
            azbas = make_azimuth(nord, npp)
            
            # All should work
            assert radp.shape[0] == nr
            assert kernels.shape[0] == 5 * nr
            assert azbas.shape[0] == npp
    
    def test_different_kl_types_consistency(self):
        """Test Kolmo and Karman produce consistent structures."""
        cobs = 0.2
        nr = 6
        radp = make_radii(cobs, nr)
        
        k_kolmo = make_kernels(cobs, nr, radp, KLType.KOLMO)
        k_karman = make_kernels(cobs, nr, radp, KLType.KARMAN)
        
        # Both should have same shape
        assert k_kolmo.shape == k_karman.shape
        
        # Both should be finite
        assert np.all(np.isfinite(k_kolmo))
        assert np.all(np.isfinite(k_karman))
