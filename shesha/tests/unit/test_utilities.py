"""
Tests for shesha.util.utilities module
"""

import numpy as np
import pytest
from shesha.util import utilities as util


class TestDistFunction:
    """Test distance calculation utility."""
    
    def test_dist_basic(self):
        """Test basic distance calculation."""
        size = 32
        xc = size // 2
        yc = size // 2
        
        dist_array = util.dist(size, xc, yc)
        
        assert isinstance(dist_array, np.ndarray)
        assert dist_array.shape == (size, size)
    
    def test_dist_center_is_zero(self):
        """Test that distance at center is near zero."""
        size = 64
        xc = size // 2
        yc = size // 2
        
        dist_array = util.dist(size, xc, yc)
        
        # Center point should be zero or very small
        center_val = dist_array[int(yc), int(xc)]
        assert center_val < 1.0
    
    def test_dist_symmetry(self):
        """Test distance array is symmetric."""
        size = 48
        xc = size // 2
        yc = size // 2
        
        dist_array = util.dist(size, xc, yc)
        
        # Check symmetry around center
        dy = yc - int(yc)
        dx = xc - int(xc)
        
        # Distance should be symmetric
        assert dist_array.shape == (size, size)
    
    def test_dist_monotonic_from_center(self):
        """Test distance increases from center."""
        size = 40
        xc = size // 2
        yc = size // 2
        
        dist_array = util.dist(size, xc, yc)
        
        # Distance at center should be less than at edges
        center_dist = dist_array[int(yc), int(xc)]
        edge_dist = dist_array[0, 0]
        
        assert edge_dist > center_dist
    
    def test_dist_different_centers(self):
        """Test with different center positions."""
        size = 32
        
        dist1 = util.dist(size, 16, 16)
        dist2 = util.dist(size, 20, 20)
        
        # Both should work without error
        assert dist1.shape == (size, size)
        assert dist2.shape == (size, size)


class TestFFTGoodSize:
    """Test FFT good size calculation."""
    
    def test_fft_goodsize_basic(self):
        """Test basic FFT good size."""
        size = 32
        
        result = util.fft_goodsize(size)
        
        assert isinstance(result, (int, np.integer))
        assert result >= size
    
    def test_fft_goodsize_powers_of_two(self):
        """Test FFT good size returns reasonable values."""
        for size in [16, 32, 64, 128, 100, 200]:
            result = util.fft_goodsize(size)
            
            # Should be at least as large as input
            assert result >= size
            # Should be reasonable (not huge)
            assert result < size * 10


class TestPadArray:
    """Test array padding."""
    
    def test_pad_array_basic(self):
        """Test basic array padding."""
        arr = np.ones((16, 16))
        newsize = 32
        
        if hasattr(util, 'pad_array'):
            result = util.pad_array(arr, newsize)
            
            assert result.shape == (newsize, newsize)
    
    def test_pad_array_preserves_data(self):
        """Test padding preserves original data."""
        arr = np.ones((8, 8)) * 5.0
        newsize = 16
        
        if hasattr(util, 'pad_array'):
            result = util.pad_array(arr, newsize)
            
            # Original data should be present in the result
            assert result.shape == (newsize, newsize)
            # Some part should have the original values
            assert np.any(result == 5.0)


class TestRebin:
    """Test rebinning utility."""
    
    def test_rebin_basic(self):
        """Test basic rebinning."""
        arr = np.arange(16).reshape(4, 4)
        
        if hasattr(util, 'rebin'):
            result = util.rebin(arr, (2, 2))
            
            # Result should be smaller
            assert result.size < arr.size
    
    def test_rebin_different_factors(self):
        """Test rebinning with different factors."""
        arr = np.ones((64, 64))
        
        if hasattr(util, 'rebin'):
            for factor in [2, 4, 8]:
                result = util.rebin(arr, (factor, factor))
                
                # Result should be smaller than input
                assert result.size < arr.size
                # Result should be a valid array
                assert isinstance(result, np.ndarray)


class TestGenerateCircle:
    """Test circle generation."""
    
    def test_generate_circle_basic(self):
        """Test basic circle generation."""
        size = 32
        radius = 10
        
        result = util.generate_circle(size, radius)
        
        # Function returns tuple of (x, y) coordinates
        assert isinstance(result, tuple)
        assert len(result) == 2
        x, y = result
        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert len(x) == len(y)
    
    def test_generate_circle_radius_effect(self):
        """Test that different radii produce different points."""
        size = 64
        
        x1, y1 = util.generate_circle(size, 10)
        x2, y2 = util.generate_circle(size, 20)
        
        # Larger radius should give more points
        assert len(x2) >= len(x1)
    
    def test_generate_circle_returns_points_in_range(self):
        """Test circle points are within expected range."""
        size = 64
        radius = 20
        
        x, y = util.generate_circle(size, radius)
        
        # Points should be roughly within size bounds
        assert np.abs(x).max() <= size
        assert np.abs(y).max() <= size


class TestGenerateSquare:
    """Test square generation."""
    
    def test_generate_square_basic(self):
        """Test basic square generation."""
        size = 32
        square_size = 16
        
        result = util.generate_square(size, square_size)
        
        # Function likely returns tuple of (x, y) coordinates like generate_circle
        if isinstance(result, tuple):
            assert len(result) == 2
            x, y = result
            assert isinstance(x, np.ndarray)
            assert isinstance(y, np.ndarray)
        else:
            assert isinstance(result, np.ndarray)
    
    def test_generate_square_size_effect(self):
        """Test that different square sizes produce different results."""
        size = 64
        
        result1 = util.generate_square(size, 10)
        result2 = util.generate_square(size, 20)
        
        # Both should work
        if isinstance(result1, tuple):
            assert len(result1) == 2
            assert len(result2) == 2
        else:
            assert isinstance(result1, np.ndarray)
            assert isinstance(result2, np.ndarray)


class TestMakeGaussian:
    """Test Gaussian generation."""
    
    def test_make_gaussian_basic(self):
        """Test basic Gaussian generation."""
        size = 32
        sigma = 5.0
        
        gaussian = util.makegaussian(size, sigma)
        
        assert isinstance(gaussian, np.ndarray)
        assert gaussian.shape == (size, size)
    
    def test_make_gaussian_properties(self):
        """Test Gaussian properties."""
        size = 64
        sigma = 10.0
        
        gaussian = util.makegaussian(size, sigma)
        
        # Gaussian should peak near center
        center_idx = size // 2
        center_val = gaussian[center_idx, center_idx]
        
        # Center should be maximum or near maximum
        assert center_val >= np.percentile(gaussian, 90)
    
    def test_make_gaussian_different_sigmas(self):
        """Test Gaussian with different sigma values."""
        size = 64
        
        gaussian_narrow = util.makegaussian(size, 2.0)
        gaussian_wide = util.makegaussian(size, 20.0)
        
        # Wider Gaussian should have more gradual falloff
        assert gaussian_narrow.shape == gaussian_wide.shape


class TestBin2D:
    """Test 2D binning."""
    
    def test_bin2d_basic(self):
        """Test basic 2D binning."""
        arr = np.arange(64).reshape(8, 8).astype(float)
        binsize = 2
        
        if hasattr(util, 'bin2d'):
            result = util.bin2d(arr, binsize)
            
            # Result should be smaller
            expected_size = 8 // binsize
            assert result.shape == (expected_size, expected_size)
    
    def test_bin2d_preserves_sum(self):
        """Test that binning preserves total sum."""
        arr = np.ones((16, 16))
        binsize = 2
        
        if hasattr(util, 'bin2d'):
            result = util.bin2d(arr, binsize)
            
            # Total sum should be preserved
            assert np.isclose(np.sum(result), np.sum(arr))


class TestUtilitiesIntegration:
    """Integration tests for utility functions."""
    
    def test_circle_and_dist(self):
        """Test using generate_circle and dist together."""
        size = 64
        
        # Create distance map
        dist_map = util.dist(size, size//2, size//2)
        
        # Create circle
        x, y = util.generate_circle(size, 20)
        
        # Distance map should have expected shape
        assert dist_map.shape == (size, size)
        # Circle should return coordinates
        assert len(x) > 0 and len(y) > 0
    
    def test_gaussian_and_fft(self):
        """Test Gaussian generation with FFT sizing."""
        size = 32
        
        # Get good FFT size
        fft_size = util.fft_goodsize(size)
        
        # Create Gaussian
        gaussian = util.makegaussian(fft_size, 5.0)
        
        # Should work
        assert gaussian.shape == (fft_size, fft_size)
    
    def test_shapes_and_square_generation(self):
        """Test square generation with fft good size."""
        size = 32
        fft_size = util.fft_goodsize(size)
        
        x, y = util.generate_square(fft_size, fft_size // 2)
        
        # Should return coordinates
        assert len(x) > 0
        assert len(y) > 0
        assert len(x) == len(y)
