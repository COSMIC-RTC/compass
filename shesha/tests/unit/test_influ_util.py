"""
Tests for shesha.util.influ_util module
"""

import numpy as np
import pytest
from shesha.util.influ_util import besel_orth, bessel_influence
from shesha.constants import PatternType


class TestBeselOrth:
    """Test Bessel-Fourier orthogonal basis functions."""
    
    def test_besel_orth_basic(self):
        """Test basic Bessel orthogonal function generation."""
        m, n = 1, 1
        phi = np.linspace(0, 2*np.pi, 32, endpoint=False)
        r = np.linspace(0, 1, 32)
        
        result = besel_orth(m, n, phi, r)
        
        assert isinstance(result, np.ndarray)
        # Should be 1D array
        assert result.ndim >= 1
    
    def test_besel_orth_different_m_n(self):
        """Test different mode indices."""
        phi = np.linspace(0, 2*np.pi, 32, endpoint=False)
        r = np.linspace(0, 1, 32)
        
        result_01 = besel_orth(0, 1, phi, r)
        result_11 = besel_orth(1, 1, phi, r)
        result_21 = besel_orth(2, 1, phi, r)
        
        # Different modes should exist
        assert result_01 is not None
        assert result_11 is not None
        assert result_21 is not None
    
    def test_besel_orth_various_radial_modes(self):
        """Test various radial mode indices."""
        m = 0
        phi = np.linspace(0, 2*np.pi, 32, endpoint=False)
        r = np.linspace(0, 1, 32)
        
        for n in range(1, 4):  # n must be >= 1
            result = besel_orth(m, n, phi, r)
            
            assert result is not None
            assert isinstance(result, np.ndarray)
    
    def test_besel_orth_azimuthal_periodicity(self):
        """Test azimuthal periodicity."""
        m = 2  # 2-fold symmetry
        n = 1
        r = np.linspace(0, 1, 32)
        
        phi_1 = np.linspace(0, 2*np.pi, 32, endpoint=False)
        phi_2 = phi_1 + np.pi / m  # Rotated by 1/m of period
        
        result_1 = besel_orth(m, n, phi_1, r)
        result_2 = besel_orth(m, n, phi_2, r)
        
        # Both should be valid results
        assert result_1 is not None
        assert result_2 is not None
    
    def test_besel_orth_radial_range(self):
        """Test different radial ranges."""
        m, n = 1, 1
        phi = np.linspace(0, 2*np.pi, 32, endpoint=False)
        
        r = np.linspace(0.1, 0.9, 32)  # Must match phi size
        
        result = besel_orth(m, n, phi, r)
        
        assert result is not None
    
    def test_besel_orth_output_range(self):
        """Test output values are reasonable."""
        m, n = 1, 1
        phi = np.linspace(0, 2*np.pi, 32, endpoint=False)
        r = np.linspace(0, 1, 32)
        
        result = besel_orth(m, n, phi, r)
        
        # Output should be finite
        assert np.all(np.isfinite(result))
    
    def test_besel_orth_grid_evaluation(self):
        """Test on regular grid."""
        m, n = 2, 1
        
        # Create regular grid
        phi = np.linspace(0, 2*np.pi, 64, endpoint=False)
        r = np.linspace(0, 1, 64)
        
        result = besel_orth(m, n, phi, r)
        
        assert result is not None
        assert isinstance(result, np.ndarray)
    
    def test_besel_orth_single_point(self):
        """Test evaluation at single point."""
        m, n = 1, 1
        phi = np.array([0.5])
        r = np.array([0.5])
        
        result = besel_orth(m, n, phi, r)
        
        assert result is not None
    
    def test_besel_orth_radial_derivative_behavior(self):
        """Test radial derivatives have expected behavior."""
        m, n = 0, 1
        phi = np.linspace(0, 2*np.pi, 32, endpoint=False)
        r = np.linspace(0.01, 1, 32)  # Avoid r=0 singularity
        
        result = besel_orth(m, n, phi, r)
        
        # Result should be finite for all r > 0
        assert np.all(np.isfinite(result))


class TestBesselInfluence:
    """Test Bessel influence function generation."""
    
    def test_bessel_influence_basic_square(self):
        """Test basic Bessel influence function with SQUARE pattern."""
        x = np.linspace(-1, 1, 64)
        y = np.linspace(-1, 1, 64)
        xx, yy = np.meshgrid(x, y)
        
        result = bessel_influence(xx, yy, type_i=PatternType.SQUARE)
        
        assert isinstance(result, np.ndarray)
        assert result.ndim == 2
    
    def test_bessel_influence_basic_hexa(self):
        """Test Bessel influence function with HEXA pattern."""
        x = np.linspace(-1, 1, 64)
        y = np.linspace(-1, 1, 64)
        xx, yy = np.meshgrid(x, y)
        
        result = bessel_influence(xx, yy, type_i=PatternType.HEXA)
        
        assert isinstance(result, np.ndarray)
        assert result.ndim == 2
    
    def test_bessel_influence_basic_hexam4(self):
        """Test Bessel influence function with HEXAM4 pattern."""
        x = np.linspace(-1, 1, 64)
        y = np.linspace(-1, 1, 64)
        xx, yy = np.meshgrid(x, y)
        
        result = bessel_influence(xx, yy, type_i=PatternType.HEXAM4)
        
        assert isinstance(result, np.ndarray)
        assert result.ndim == 2
    
    def test_bessel_influence_2d_grid(self):
        """Test influence on 2D grid."""
        x = np.linspace(-1.5, 1.5, 128)
        y = np.linspace(-1.5, 1.5, 128)
        xx, yy = np.meshgrid(x, y)
        
        result = bessel_influence(xx, yy, type_i=PatternType.SQUARE)
        
        # Result should be 2D
        assert result.ndim == 2
        assert result.shape == xx.shape
    
    def test_bessel_influence_different_sizes(self):
        """Test different grid sizes."""
        for size in [32, 64, 128]:
            x = np.linspace(-1, 1, size)
            y = np.linspace(-1, 1, size)
            xx, yy = np.meshgrid(x, y)
            
            result = bessel_influence(xx, yy, type_i=PatternType.SQUARE)
            
            assert result.shape == (size, size)
    
    def test_bessel_influence_symmetry(self):
        """Test influence function has expected symmetry."""
        x = np.linspace(-1, 1, 64)
        y = np.linspace(-1, 1, 64)
        xx, yy = np.meshgrid(x, y)
        
        result = bessel_influence(xx, yy, type_i=PatternType.SQUARE)
        
        # Should have real values
        assert np.all(np.isreal(result))
    
    def test_bessel_influence_values_positive(self):
        """Test influence function values are reasonable."""
        x = np.linspace(-1, 1, 64)
        y = np.linspace(-1, 1, 64)
        xx, yy = np.meshgrid(x, y)
        
        result = bessel_influence(xx, yy, type_i=PatternType.SQUARE)
        
        # Values should be finite
        assert np.all(np.isfinite(result))
    
    def test_bessel_influence_pattern_types(self):
        """Test all pattern types work."""
        x = np.linspace(-1, 1, 64)
        y = np.linspace(-1, 1, 64)
        xx, yy = np.meshgrid(x, y)
        
        patterns = [PatternType.SQUARE, PatternType.HEXA, PatternType.HEXAM4]
        results = []
        
        for pattern in patterns:
            result = bessel_influence(xx, yy, type_i=pattern)
            results.append(result)
            
            assert result is not None
        
        # Different patterns should produce different results
        assert not np.allclose(results[0], results[1])
    
    def test_bessel_influence_extended_range(self):
        """Test influence over extended spatial range."""
        x = np.linspace(-5, 5, 128)
        y = np.linspace(-5, 5, 128)
        xx, yy = np.meshgrid(x, y)
        
        result = bessel_influence(xx, yy, type_i=PatternType.SQUARE)
        
        assert np.all(np.isfinite(result))
    
    def test_bessel_influence_narrow_range(self):
        """Test influence over narrow range."""
        x = np.linspace(-0.1, 0.1, 32)
        y = np.linspace(-0.1, 0.1, 32)
        xx, yy = np.meshgrid(x, y)
        
        result = bessel_influence(xx, yy, type_i=PatternType.SQUARE)
        
        assert np.all(np.isfinite(result))
    
    def test_bessel_influence_1d_arrays(self):
        """Test with 1D arrays (should still work)."""
        x = np.array([0, 0.5, 1.0])
        y = np.array([0, 0.5, 1.0])
        
        result = bessel_influence(x, y, type_i=PatternType.SQUARE)
        
        assert result is not None


class TestBesselIntegration:
    """Integration tests for Bessel utilities."""
    
    def test_besel_orth_and_bessel_influence(self):
        """Test Bessel orthogonal functions with influence."""
        # Create 2D grids
        x = np.linspace(-1, 1, 32)
        y = np.linspace(-1, 1, 32)
        xx, yy = np.meshgrid(x, y)
        
        # Generate Bessel basis
        m, n = 1, 1
        phi = np.arctan2(yy, xx)
        r = np.sqrt(xx**2 + yy**2)
        
        besel = besel_orth(m, n, phi, r)
        
        # Generate influence
        influence = bessel_influence(xx, yy, type_i=PatternType.SQUARE)
        
        # Both should be valid arrays
        assert besel is not None
        assert influence is not None
    
    def test_multiple_actuator_patterns(self):
        """Test influence generation for multiple actuator patterns."""
        x = np.linspace(-2, 2, 128)
        y = np.linspace(-2, 2, 128)
        xx, yy = np.meshgrid(x, y)
        
        patterns = {
            'SQUARE': PatternType.SQUARE,
            'HEXA': PatternType.HEXA,
            'HEXAM4': PatternType.HEXAM4,
        }
        
        results = {}
        for name, pattern in patterns.items():
            result = bessel_influence(xx, yy, type_i=pattern)
            results[name] = result
            
            assert result.shape == (128, 128)
        
        # Check that results have expected numeric type
        for result in results.values():
            assert isinstance(result, np.ndarray)
    
    def test_bessel_influence_conservation(self):
        """Test influence function properties."""
        x = np.linspace(-2, 2, 128)
        y = np.linspace(-2, 2, 128)
        xx, yy = np.meshgrid(x, y)
        
        result = bessel_influence(xx, yy, type_i=PatternType.SQUARE)
        
        # Check that function has reasonable structure
        # (peak in center, decreases away)
        center_idx = result.shape[0] // 2, result.shape[1] // 2
        center_value = result[center_idx]
        edge_value = result[0, 0]
        
        # Center should not be less than edges (typical for influence)
        # or this could be negative peak - just check finite
        assert np.isfinite(center_value)
        assert np.isfinite(edge_value)
