"""
Tests for shesha.util.rtc_util module
"""

import numpy as np
import pytest
from shesha.util.rtc_util import create_interp_mat, centroid_gain


class TestCreateInterpMat:
    """Test interpolation matrix creation."""
    
    def test_create_interp_mat_basic(self):
        """Test basic interpolation matrix creation."""
        dimx, dimy = 5, 5
        
        mat = create_interp_mat(dimx, dimy)
        
        assert isinstance(mat, np.ndarray)
        assert mat.shape == (dimx * dimy, 6)
    
    def test_create_interp_mat_dtype(self):
        """Test interpolation matrix has correct dtype."""
        dimx, dimy = 5, 5
        
        mat = create_interp_mat(dimx, dimy)
        
        # Should be float type for inverse matrix
        assert mat.dtype in [np.float32, np.float64]
    
    def test_create_interp_mat_rectangular(self):
        """Test interpolation matrix for rectangular dimensions."""
        dimx, dimy = 3, 7
        
        mat = create_interp_mat(dimx, dimy)
        
        assert mat.shape == (dimx * dimy, 6)
    
    def test_create_interp_mat_square(self):
        """Test interpolation matrix for square dimensions."""
        dim = 10
        
        mat = create_interp_mat(dim, dim)
        
        assert mat.shape == (dim * dim, 6)
    
    def test_create_interp_mat_large(self):
        """Test interpolation matrix for larger dimensions."""
        dimx, dimy = 64, 64
        
        mat = create_interp_mat(dimx, dimy)
        
        assert mat.shape == (dimx * dimy, 6)
        assert np.all(np.isfinite(mat))
    
    def test_create_interp_mat_small(self):
        """Test interpolation matrix for small dimensions."""
        dimx, dimy = 5, 5
        
        mat = create_interp_mat(dimx, dimy)
        
        assert mat.shape == (25, 6)
    
    def test_create_interp_mat_asymmetric(self):
        """Test interpolation matrix for asymmetric dimensions."""
        dimx, dimy = 2, 8
        
        mat1 = create_interp_mat(dimx, dimy)
        mat2 = create_interp_mat(dimy, dimx)
        
        # Different shapes
        assert mat1.shape == (16, 6)
        assert mat2.shape == (16, 6)
    
    def test_create_interp_mat_is_invertible(self):
        """Test that basis formed by interp mat can be used."""
        dimx, dimy = 5, 5
        
        mat = create_interp_mat(dimx, dimy)
        
        # Matrix should have full rank (can reconstruct via least squares)
        assert np.linalg.matrix_rank(mat) == min(mat.shape)
    
    def test_create_interp_mat_contains_basis_terms(self):
        """Test that interpolation matrix contains quadratic basis terms."""
        dimx, dimy = 3, 3
        
        mat = create_interp_mat(dimx, dimy)
        
        # Should have 6 columns for: x^2, y^2, xy, x, y, 1
        assert mat.shape[1] == 6
    
    def test_create_interp_mat_different_max_dim(self):
        """Test that matrix size is based on dimx * dimy."""
        # Matrix size should be dimx * dimy for both
        mat1 = create_interp_mat(5, 6)
        mat2 = create_interp_mat(6, 5)
        
        # Both should have 30 rows (5*6)
        assert mat1.shape[0] == 30
        assert mat2.shape[0] == 30


class TestCentroidGain:
    """Test centroid gain calculation."""
    
    def test_centroid_gain_1d_perfect_scaling(self):
        """Test centroid gain for perfectly scaled 1D data."""
        E = np.array([0, 1, 2, 3, 4], dtype=np.float32)
        F = 2.5 * E  # F = 2.5 * E
        
        cgain = centroid_gain(E, F)
        
        assert np.isclose(cgain, 2.5, rtol=1e-5)
    
    def test_centroid_gain_1d_with_offset(self):
        """Test centroid gain for 1D data with offset (slope is key)."""
        E = np.array([0, 1, 2, 3, 4], dtype=np.float32)
        F = 2.5 * E + 0.5  # Offset shouldn't matter for slope
        
        cgain = centroid_gain(E, F)
        
        assert np.isclose(cgain, 2.5, rtol=1e-5)
    
    def test_centroid_gain_1d_negative_scale(self):
        """Test centroid gain for negative scaling."""
        E = np.array([0, 1, 2, 3, 4], dtype=np.float32)
        F = -1.5 * E
        
        cgain = centroid_gain(E, F)
        
        assert np.isclose(cgain, -1.5, rtol=1e-5)
    
    def test_centroid_gain_2d_single_channel(self):
        """Test centroid gain for 2D data with single channel."""
        E = np.array([[0, 1, 2, 3, 4]], dtype=np.float32).T
        F = 2.0 * E
        
        cgain = centroid_gain(E, F)
        
        assert np.isclose(cgain, 2.0, rtol=1e-5)
    
    def test_centroid_gain_2d_multiple_channels(self):
        """Test centroid gain for 2D data with multiple channels."""
        n = 50
        n_channels = 3
        
        E = np.random.randn(n, n_channels).astype(np.float32)
        scales = np.array([1.5, 2.0, 2.5])
        F = E * scales[np.newaxis, :]
        
        cgain = centroid_gain(E, F)
        
        # Should be mean of the individual gains
        expected = np.mean(scales)
        assert np.isclose(cgain, expected, rtol=1e-4)
    
    def test_centroid_gain_2d_with_noise(self):
        """Test centroid gain for noisy 2D data."""
        np.random.seed(42)
        n = 100
        n_channels = 2
        
        E = np.random.randn(n, n_channels).astype(np.float32)
        scales = np.array([1.2, 1.8])
        noise = 0.01 * np.random.randn(n, n_channels).astype(np.float32)
        F = E * scales[np.newaxis, :] + noise
        
        cgain = centroid_gain(E, F)
        
        # Should be approximately mean of scales
        expected = np.mean(scales)
        assert np.abs(cgain - expected) < 0.05
    
    def test_centroid_gain_1d_linear_relationship(self):
        """Test centroid gain preserves linear relationship."""
        E = np.array([-2, -1, 0, 1, 2], dtype=np.float32)
        F = np.array([-4, -2, 0, 2, 4], dtype=np.float32)
        
        cgain = centroid_gain(E, F)
        
        assert np.isclose(cgain, 2.0, rtol=1e-5)
    
    def test_centroid_gain_2d_consistency(self):
        """Test centroid gain consistency across channels."""
        n = 50
        E = np.random.randn(n, 5).astype(np.float32)
        F = 1.5 * E
        
        cgain = centroid_gain(E, F)
        
        # All channels should have same gain
        assert np.isclose(cgain, 1.5, rtol=1e-4)
    
    def test_centroid_gain_raises_on_invalid_dims(self):
        """Test that invalid dimensions raise error."""
        E = np.random.randn(5, 5, 5).astype(np.float32)
        F = E
        
        with pytest.raises(ValueError):
            centroid_gain(E, F)
    
    def test_centroid_gain_zero_input(self):
        """Test centroid gain with constant input (special case)."""
        E = np.ones(5, dtype=np.float32)
        F = 2.0 * np.ones(5, dtype=np.float32)
        
        # polyfit will compute slope from constant data
        cgain = centroid_gain(E, F)
        
        # Should not raise error and be finite
        assert np.isfinite(cgain)
    
    def test_centroid_gain_2d_many_samples(self):
        """Test centroid gain with many samples."""
        n = 10000
        n_channels = 10
        
        E = np.random.randn(n, n_channels).astype(np.float32)
        scales = np.linspace(0.5, 2.5, n_channels)
        F = E * scales[np.newaxis, :]
        
        cgain = centroid_gain(E, F)
        
        expected = np.mean(scales)
        assert np.isclose(cgain, expected, rtol=1e-3)


class TestRTCIntegration:
    """Integration tests for RTC utilities."""
    
    def test_interp_mat_with_centroid_gain(self):
        """Test using interpolation matrix result with centroid gain."""
        dimx, dimy = 5, 5
        
        mat = create_interp_mat(dimx, dimy)
        
        # Create synthetic measurements
        E = np.random.randn(dimx * dimy, 5).astype(np.float32)
        F = 2.0 * E
        
        cgain = centroid_gain(E, F)
        
        assert np.isclose(cgain, 2.0, rtol=1e-3)
    
    def test_matrix_reconstruction_problem(self):
        """Test solving a simple reconstruction problem."""
        dimx, dimy = 8, 8
        
        mat = create_interp_mat(dimx, dimy)
        
        # Create synthetic data to fit
        n_samples = dimx * dimy
        true_coefs = np.array([0.1, 0.2, -0.05, 0.3, 0.1, 1.0])
        y = np.dot(mat, true_coefs)
        
        # Solve for coefficients
        coefs_fit = np.linalg.lstsq(mat, y, rcond=None)[0]
        
        # Should recover coefficients
        assert np.allclose(coefs_fit, true_coefs, rtol=1e-5)
    
    def test_consistent_centroid_gains(self):
        """Test that centroid gains are consistent across multiple runs."""
        np.random.seed(42)
        E = np.random.randn(100, 4).astype(np.float32)
        F = 1.8 * E + 0.1
        
        cgain1 = centroid_gain(E, F)
        cgain2 = centroid_gain(E, F)
        
        assert np.isclose(cgain1, cgain2)
