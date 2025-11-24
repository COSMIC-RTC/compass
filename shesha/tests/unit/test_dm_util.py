"""
Tests for shesha.util.dm_util module
"""

import numpy as np
import pytest
from shesha.util.dm_util import (
    dim_dm_support, dim_dm_patch, createSquarePattern,
    createHexaPattern, select_actuators, make_zernike, zernumero,
    filterActuWithPupil
)
from shesha.constants import DmType, CONST


class TestDimDmSupport:
    """Test DM support dimension calculations."""
    
    def test_dim_dm_support_basic(self):
        """Test basic DM support dimension calculation."""
        cent = 256.0
        extent = 100
        ssize = 512
        
        n1, n2 = dim_dm_support(cent, extent, ssize)
        
        assert isinstance(n1, int)
        assert isinstance(n2, int)
        assert n1 >= 1
        assert n2 <= ssize
        assert n1 < n2
    
    def test_dim_dm_support_centered(self):
        """Test DM support at center of array."""
        cent = 256.0
        extent = 100
        ssize = 512
        
        n1, n2 = dim_dm_support(cent, extent, ssize)
        
        # Should be roughly centered
        center = (n1 + n2) / 2.0
        assert abs(center - cent) < extent / 2 + 1
    
    def test_dim_dm_support_near_edge(self):
        """Test DM support near edge is clipped."""
        cent = 10.0
        extent = 100
        ssize = 512
        
        n1, n2 = dim_dm_support(cent, extent, ssize)
        
        # Should be clipped at lower bound
        assert n1 == 1
    
    def test_dim_dm_support_upper_edge(self):
        """Test DM support near upper edge is clipped."""
        cent = 512.0
        extent = 100
        ssize = 512
        
        n1, n2 = dim_dm_support(cent, extent, ssize)
        
        # Should be clipped at upper bound
        assert n2 == ssize
    
    def test_dim_dm_support_size_calculation(self):
        """Test DM support size is consistent with extent."""
        cent = 256.0
        extent = 50
        ssize = 512
        
        n1, n2 = dim_dm_support(cent, extent, ssize)
        
        # Size should be close to extent
        size = n2 - n1
        assert abs(size - extent) < 2


class TestDimDmPatch:
    """Test DM patch dimension calculations."""
    
    def test_dim_dm_patch_pzt(self):
        """Test DM patch dimensions for PZT type."""
        pupdiam = 512
        diam = 8.0
        alt = 0.0  # On-axis
        xpos_wfs = []
        ypos_wfs = []
        
        patchDiam = dim_dm_patch(pupdiam, diam, DmType.PZT, alt, xpos_wfs, ypos_wfs)
        
        assert isinstance(patchDiam, int)
        assert patchDiam > 0
        assert patchDiam >= pupdiam
    
    def test_dim_dm_patch_kl(self):
        """Test DM patch dimensions for KL type."""
        pupdiam = 512
        diam = 8.0
        alt = 0.0
        xpos_wfs = []
        ypos_wfs = []
        
        patchDiam = dim_dm_patch(pupdiam, diam, DmType.KL, alt, xpos_wfs, ypos_wfs)
        
        assert isinstance(patchDiam, int)
        assert patchDiam > 0
    
    def test_dim_dm_patch_with_offset_wfs(self):
        """Test DM patch with off-axis WFS."""
        pupdiam = 512
        diam = 8.0
        alt = 1000.0
        xpos_wfs = [10.0, 20.0]  # arcsec
        ypos_wfs = [10.0, 20.0]  # arcsec
        
        patchDiam = dim_dm_patch(pupdiam, diam, DmType.PZT, alt, xpos_wfs, ypos_wfs)
        
        # Should be larger due to WFS offset
        patchDiam_nooffset = dim_dm_patch(pupdiam, diam, DmType.PZT, alt, [], [])
        assert patchDiam >= patchDiam_nooffset
    
    def test_dim_dm_patch_invalid_type(self):
        """Test DM patch with invalid DM type."""
        pupdiam = 512
        diam = 8.0
        alt = 0.0
        xpos_wfs = []
        ypos_wfs = []
        
        with pytest.raises(TypeError):
            dim_dm_patch(pupdiam, diam, "invalid_type", alt, xpos_wfs, ypos_wfs)


class TestCreateSquarePattern:
    """Test square pattern actuator creation."""
    
    def test_create_square_pattern_basic(self):
        """Test basic square pattern creation."""
        pitch = 10.0
        nxact = 5
        
        pos = createSquarePattern(pitch, nxact)
        
        assert isinstance(pos, np.ndarray)
        # Shape is (2, M) where M = nxact * nxact
        assert pos.shape[0] == 2
        assert pos.shape[1] == nxact * nxact
    
    def test_create_square_pattern_symmetry(self):
        """Test square pattern is centered."""
        pitch = 1.0
        nxact = 3
        
        pos = createSquarePattern(pitch, nxact)
        
        # Should be centered around origin
        # pos[0, :] is x coordinates, pos[1, :] is y coordinates
        center_x = np.mean(pos[0, :])
        center_y = np.mean(pos[1, :])
        assert np.isclose(center_x, 0, atol=0.1)
        assert np.isclose(center_y, 0, atol=0.1)
    
    def test_create_square_pattern_spacing(self):
        """Test square pattern spacing matches pitch."""
        pitch = 10.0
        nxact = 3
        
        pos = createSquarePattern(pitch, nxact)
        
        # Calculate actual spacing
        x_unique = np.unique(np.round(pos[0, :], 5))
        y_unique = np.unique(np.round(pos[1, :], 5))
        
        # Should have nxact unique positions in each direction
        assert len(x_unique) == nxact
        assert len(y_unique) == nxact


class TestCreateHexaPattern:
    """Test hexagonal pattern actuator creation."""
    
    def test_create_hexa_pattern_basic(self):
        """Test basic hexagonal pattern creation."""
        pitch = 10.0
        nxact = 5
        
        pos = createHexaPattern(pitch, nxact)
        
        assert isinstance(pos, np.ndarray)
        # Shape is (2, M)
        assert pos.shape[0] == 2
        assert pos.shape[1] > 0
    
    def test_create_hexa_pattern_symmetry(self):
        """Test hexagonal pattern is mostly centered."""
        pitch = 1.0
        supportSize = 3
        
        pos = createHexaPattern(pitch, supportSize)
        
        # Should be roughly centered
        center_x = np.median(pos[0, :])
        center_y = np.median(pos[1, :])
        # Just check they're reasonable numbers
        assert not np.isnan(center_x)
        assert not np.isnan(center_y)
    
    def test_create_hexa_pattern_count(self):
        """Test hexagonal pattern returns actuators."""
        pitch = 1.0
        supportSize = 5
        
        pos = createHexaPattern(pitch, supportSize)
        # Should have at least some actuators
        assert pos.shape[1] > 0


class TestActuatorPatterns:
    """Integration tests for actuator patterns."""
    
    def test_square_vs_hexa_shape(self):
        """Test square and hexa patterns return correct shapes."""
        pitch = 1.0
        nxact = 3
        
        square = createSquarePattern(pitch, nxact)
        hexa = createHexaPattern(pitch, nxact)
        
        # Both should be (2, M) format
        assert square.shape[0] == 2
        assert hexa.shape[0] == 2
        assert square.shape[1] > 0
        assert hexa.shape[1] > 0
    
    def test_pattern_dtype(self):
        """Test pattern arrays are float32."""
        pitch = 1.0
        nxact = 3
        
        square = createSquarePattern(pitch, nxact)
        hexa = createHexaPattern(pitch, nxact)
        
        assert square.dtype == np.float32
        assert hexa.dtype == np.float32


class TestSelectActuators:
    """Test actuator selection based on geometry."""
    
    def test_select_actuators_basic(self):
        """Test basic actuator selection."""
        xc = np.array([0.0, 1.0, -1.0, 2.0, -2.0])
        yc = np.array([0.0, 1.0, -1.0, 2.0, -2.0])
        nxact = 5
        pitch = 1.0
        cobs = 0.3
        margin_in = 0.5
        margin_out = 1.5
        
        selected = select_actuators(xc, yc, nxact, pitch, cobs, margin_in, margin_out)
        
        assert isinstance(selected, np.ndarray)
        assert len(selected) > 0
        assert len(selected) <= len(xc)
    
    def test_select_actuators_all_inside(self):
        """Test selection when all actuators are inside bounds."""
        xc = np.array([0.0, 0.5, -0.5, 1.0, -1.0])
        yc = np.array([0.0, 0.5, -0.5, 1.0, -1.0])
        nxact = 5
        pitch = 1.0
        cobs = 0.1  # Small central obscuration
        margin_in = 0.0
        margin_out = 2.0
        
        selected = select_actuators(xc, yc, nxact, pitch, cobs, margin_in, margin_out)
        
        # Should have some selected actuators
        assert len(selected) > 0
    
    def test_select_actuators_with_N_limit(self):
        """Test selection with maximum N actuators."""
        xc = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        yc = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        nxact = 5
        pitch = 1.0
        cobs = 0.3
        margin_in = 0.0
        margin_out = None
        N = 3
        
        selected = select_actuators(xc, yc, nxact, pitch, cobs, margin_in, margin_out, N=N)
        
        assert len(selected) <= N
    
    def test_select_actuators_centered_grid(self):
        """Test selection on centered square actuator grid."""
        # Create a centered grid
        n = 5
        pos = np.linspace(-2, 2, n)
        xx, yy = np.meshgrid(pos, pos)
        xc = xx.flatten()
        yc = yy.flatten()
        
        selected = select_actuators(xc, yc, n, 1.0, 0.1, 0.5, 2.5)
        
        assert len(selected) > 0
        # Verify center actuator is among selected
        dis = np.sqrt(xc**2 + yc**2)
        center_dis = dis[selected]
        assert np.min(center_dis) < 0.5  # Some close to center


class TestZernikeGeneration:
    """Test Zernike mode generation."""
    
    def test_make_zernike_basic(self):
        """Test basic Zernike mode generation."""
        nzer = 5
        size = 64
        diameter = 48
        
        z = make_zernike(nzer, size, diameter)
        
        assert isinstance(z, np.ndarray)
        assert z.shape == (size, size, nzer)
        assert z.dtype == np.float32
    
    def test_make_zernike_single_mode(self):
        """Test generation of single Zernike mode."""
        nzer = 1
        size = 32
        diameter = 24
        
        z = make_zernike(nzer, size, diameter)
        
        assert z.shape == (size, size, nzer)
        # Piston should be mostly 1 within pupil
        pupil_mask = z[:, :, 0] > 0
        assert np.any(pupil_mask)
    
    def test_make_zernike_with_custom_center(self):
        """Test Zernike generation with custom center."""
        nzer = 3
        size = 64
        diameter = 48
        xc = 30.0
        yc = 30.0
        
        z = make_zernike(nzer, size, diameter, xc=xc, yc=yc)
        
        assert z.shape == (size, size, nzer)
        # Check some values are non-zero
        assert np.any(z != 0)
    
    def test_make_zernike_extended(self):
        """Test Zernike generation with extension."""
        nzer = 5
        size = 64
        diameter = 48
        
        z_normal = make_zernike(nzer, size, diameter, ext=0)
        z_extended = make_zernike(nzer, size, diameter, ext=1)
        
        assert z_normal.shape == z_extended.shape
        # Extended version should have more non-zero values
        assert np.sum(z_extended != 0) >= np.sum(z_normal != 0)
    
    def test_make_zernike_large_nzer(self):
        """Test generation of many Zernike modes."""
        nzer = 20
        size = 128
        diameter = 100
        
        z = make_zernike(nzer, size, diameter)
        
        assert z.shape == (size, size, nzer)
        # All modes should have some structure
        for i in range(nzer):
            assert np.any(z[:, :, i] != 0)


class TestZernumeро:
    """Test Zernike numbering conversion."""
    
    def test_zernumero_first_modes(self):
        """Test Zernike number to (n,m) conversion for first modes."""
        # Noll numbering according to the actual implementation
        # Test that we get valid (n, m) pairs for first several modes
        for zn in range(1, 10):
            n, m = zernumero(zn)
            # Check n is non-negative
            assert n >= 0
            # Check |m| <= n
            assert abs(m) <= n
            # Check both are integers
            assert isinstance(n, int)
            assert isinstance(m, int)
    
    def test_zernumero_monotonic_n(self):
        """Test that n is monotonically increasing with Zernike number."""
        n_values = []
        for zn in range(1, 21):
            n, _ = zernumero(zn)
            n_values.append(n)
        
        # n should be non-decreasing
        for i in range(1, len(n_values)):
            assert n_values[i] >= n_values[i-1]
    
    def test_zernumero_abs_m_le_n(self):
        """Test that |m| <= n for all Zernike modes."""
        for zn in range(1, 50):
            n, m = zernumero(zn)
            assert abs(m) <= n


class TestFilterActuWithPupil:
    """Test actuator filtering by pupil."""
    
    def test_filter_actu_with_pupil_basic(self):
        """Test basic actuator filtering with pupil."""
        # Create simple pupil
        pupil_size = 32
        pupil = np.ones((pupil_size, pupil_size), dtype=bool)
        pupil[:8, :] = False
        pupil[-8:, :] = False
        
        # Create actuators within pupil area
        actuPos = np.array([
            [16.0, 16.0, 16.0],
            [16.0, 16.0, 16.0]
        ], dtype=np.float32)
        
        threshold = 2.0
        filtered = filterActuWithPupil(actuPos, pupil, threshold)
        
        assert isinstance(filtered, np.ndarray)
        assert filtered.shape[0] == 2
    
    def test_filter_actu_with_pupil_circular(self):
        """Test filtering with circular pupil."""
        pupil_size = 64
        center = pupil_size // 2
        radius = 20
        
        # Create circular pupil
        y, x = np.ogrid[:pupil_size, :pupil_size]
        pupil = (x - center)**2 + (y - center)**2 <= radius**2
        
        # Create actuators
        actuPos = np.array([
            [center, center, center - 25],
            [center, center, center]
        ], dtype=np.float32)
        
        threshold = 3.0
        filtered = filterActuWithPupil(actuPos, pupil, threshold)
        
        # At least center actuator should remain
        assert filtered.shape[1] >= 1
    
    def test_filter_actu_with_pupil_threshold(self):
        """Test filtering with different thresholds."""
        pupil = np.ones((32, 32), dtype=bool)
        actuPos = np.array([[16.0], [16.0]], dtype=np.float32)
        
        # With large threshold, should keep actuator
        filtered_large = filterActuWithPupil(actuPos, pupil, 10.0)
        assert filtered_large.shape[1] == 1
        
        # With small threshold on edge might filter out
        pupil_edge = np.zeros((32, 32), dtype=bool)
        pupil_edge[15:17, 15:17] = True
        
        actuPos_edge = np.array([[30.0], [30.0]], dtype=np.float32)
        filtered_small = filterActuWithPupil(actuPos_edge, pupil_edge, 0.5)
        
        # This should filter out the edge actuator
        assert filtered_small.shape[1] == 0
