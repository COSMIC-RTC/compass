"""
Tests for shesha.constants module
"""

import numpy as np
import pytest
from shesha.constants import (
    CONST, check_enum, DmType, PatternType, KLType, InfluType,
    ApertureType, SpiderType, CentroiderType, WFSType, ControllerType
)


class TestCONST:
    """Test constants conversion factors."""
    
    def test_rad_to_arcsec(self):
        """Test RAD2ARCSEC conversion factor."""
        assert CONST.RAD2ARCSEC == 3600.0 * 360.0 / (2 * np.pi)
        # One radian should convert properly
        one_rad_arcsec = 1.0 * CONST.RAD2ARCSEC
        assert one_rad_arcsec > 0
    
    def test_arcsec_to_rad(self):
        """Test ARCSEC2RAD conversion factor."""
        assert CONST.ARCSEC2RAD == 2.0 * np.pi / (360.0 * 3600.0)
        # One arcsec should convert properly
        one_arcsec_rad = 1.0 * CONST.ARCSEC2RAD
        assert one_arcsec_rad > 0
    
    def test_rad_to_deg(self):
        """Test RAD2DEG conversion factor."""
        assert CONST.RAD2DEG == 180.0 / np.pi
        # Pi radians should be 180 degrees
        pi_in_deg = np.pi * CONST.RAD2DEG
        assert np.isclose(pi_in_deg, 180.0)
    
    def test_deg_to_rad(self):
        """Test DEG2RAD conversion factor."""
        assert CONST.DEG2RAD == np.pi / 180.0
        # 180 degrees should be pi radians
        deg_to_rad = 180.0 * CONST.DEG2RAD
        assert np.isclose(deg_to_rad, np.pi)
    
    def test_conversion_roundtrip_rad_arcsec(self):
        """Test roundtrip conversion rad -> arcsec -> rad."""
        original = 0.001  # radians
        converted = original * CONST.RAD2ARCSEC
        back = converted * CONST.ARCSEC2RAD
        assert np.isclose(original, back)
    
    def test_conversion_roundtrip_deg_rad(self):
        """Test roundtrip conversion deg -> rad -> deg."""
        original = 45.0  # degrees
        converted = original * CONST.DEG2RAD
        back = converted * CONST.RAD2DEG
        assert np.isclose(original, back)


class TestCheckEnum:
    """Test check_enum validation function."""
    
    def test_check_enum_valid_string(self):
        """Test check_enum with valid enum value."""
        result = check_enum(DmType, "pzt")
        assert result == "pzt"
    
    def test_check_enum_invalid_string(self):
        """Test check_enum with invalid enum value."""
        with pytest.raises(ValueError):
            check_enum(DmType, "invalid_type")
    
    def test_check_enum_non_string_input(self):
        """Test check_enum with non-string input."""
        with pytest.raises(ValueError):
            check_enum(DmType, 123)
    
    def test_check_enum_bytes_input(self):
        """Test check_enum with bytes input."""
        with pytest.raises(ValueError):
            check_enum(DmType, b"pzt")


class TestDmType:
    """Test DM type constants."""
    
    def test_dm_type_pzt(self):
        """Test PZT DM type."""
        assert DmType.PZT == "pzt"
    
    def test_dm_type_tt(self):
        """Test TT DM type."""
        assert DmType.TT == "tt"
    
    def test_dm_type_kl(self):
        """Test KL DM type."""
        assert DmType.KL == "kl"
    
    def test_dm_type_values_unique(self):
        """Test DM types are unique."""
        types = [DmType.PZT, DmType.TT, DmType.KL]
        assert len(types) == len(set(types))


class TestPatternType:
    """Test pattern type constants."""
    
    def test_pattern_type_square(self):
        """Test SQUARE pattern type."""
        assert PatternType.SQUARE == "square"
    
    def test_pattern_type_hexa(self):
        """Test HEXA pattern type."""
        assert PatternType.HEXA == "hexa"
    
    def test_pattern_type_hexam4(self):
        """Test HEXAM4 pattern type."""
        assert PatternType.HEXAM4 == "hexaM4"


class TestKLType:
    """Test KL type constants."""
    
    def test_kl_type_kolmo(self):
        """Test KOLMO KL type."""
        assert KLType.KOLMO == "kolmo"
    
    def test_kl_type_karman(self):
        """Test KARMAN KL type."""
        assert KLType.KARMAN == "karman"


class TestInfluType:
    """Test influence function type constants."""
    
    def test_influ_type_default(self):
        """Test DEFAULT influence type."""
        assert InfluType.DEFAULT == "default"
    
    def test_influ_type_radial_schwartz(self):
        """Test RADIALSCHWARTZ influence type."""
        assert InfluType.RADIALSCHWARTZ == "radialSchwartz"
    
    def test_influ_type_square_schwartz(self):
        """Test SQUARESCHWARTZ influence type."""
        assert InfluType.SQUARESCHWARTZ == "squareSchwartz"
    
    def test_influ_type_blacknutt(self):
        """Test BLACKNUTT influence type."""
        assert InfluType.BLACKNUTT == "blacknutt"
    
    def test_influ_type_gaussian(self):
        """Test GAUSSIAN influence type."""
        assert InfluType.GAUSSIAN == "gaussian"
    
    def test_influ_type_bessel(self):
        """Test BESSEL influence type."""
        assert InfluType.BESSEL == "bessel"
    
    def test_influ_type_petal(self):
        """Test PETAL influence type."""
        assert InfluType.PETAL == "petal"


class TestApertureType:
    """Test aperture type constants."""
    
    def test_aperture_types_exist(self):
        """Test aperture types are defined."""
        assert hasattr(ApertureType, 'EELT')
        assert hasattr(ApertureType, 'EELT_NOMINAL')
    
    def test_aperture_type_values_unique(self):
        """Test aperture types have unique values."""
        types = [ApertureType.EELT, ApertureType.EELT_NOMINAL]
        assert len(types) == len(set(types))


class TestSpiderType:
    """Test spider type constants."""
    
    def test_spider_types_exist(self):
        """Test spider types are defined."""
        # These should be defined in the constants module
        assert hasattr(SpiderType, 'FOUR') or hasattr(SpiderType, 'NO')


class TestCentroiderType:
    """Test centroider type constants."""
    
    def test_centroider_types_exist(self):
        """Test centroider types are defined."""
        assert hasattr(CentroiderType, 'COG')


class TestWFSType:
    """Test WFS type constants."""
    
    def test_wfs_types_exist(self):
        """Test WFS types are defined."""
        assert hasattr(WFSType, 'SH')


class TestControllerType:
    """Test controller type constants."""
    
    def test_controller_types_exist(self):
        """Test controller types are defined."""
        assert hasattr(ControllerType, 'GENERIC')
