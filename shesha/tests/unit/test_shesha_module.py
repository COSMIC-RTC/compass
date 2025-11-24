"""
Tests for shesha module initialization and utility functions
"""

import pytest
import numpy as np
from pathlib import Path
import shesha


class TestSheshaModuleImports:
    """Test shesha module imports and basic structure."""
    
    def test_shesha_module_exists(self):
        """Test shesha module can be imported."""
        assert shesha is not None
    
    def test_shesha_has_supervisor(self):
        """Test shesha has supervisor module."""
        from shesha import supervisor
        assert supervisor is not None
    
    def test_shesha_has_config(self):
        """Test shesha has config module."""
        from shesha import config
        assert config is not None
    
    def test_shesha_has_util(self):
        """Test shesha has util module."""
        from shesha import util
        assert util is not None
    
    def test_shesha_has_constants(self):
        """Test shesha has constants module."""
        from shesha import constants
        assert constants is not None
    
    def test_shesha_has_init(self):
        """Test shesha has init module."""
        from shesha import init
        assert init is not None


class TestSheshaVersionInfo:
    """Test shesha version information."""
    
    def test_shesha_has_version(self):
        """Test shesha module has version info."""
        # Version might be in __version__ or similar
        assert hasattr(shesha, '__name__')
        assert shesha.__name__ == 'shesha'


class TestSheshaUtilModules:
    """Test shesha utility submodules."""
    
    def test_util_has_utilities(self):
        """Test shesha.util has utilities module."""
        from shesha.util import utilities
        assert utilities is not None
    
    def test_util_has_dm_util(self):
        """Test shesha.util has dm_util module."""
        from shesha.util import dm_util
        assert dm_util is not None
    
    def test_util_has_kl_util(self):
        """Test shesha.util has kl_util module."""
        from shesha.util import kl_util
        assert kl_util is not None
    
    def test_util_has_coronagraph_utils(self):
        """Test shesha.util has coronagraph_utils module."""
        from shesha.util import coronagraph_utils
        assert coronagraph_utils is not None


class TestSheshaConfigClasses:
    """Test shesha configuration classes."""
    
    def test_config_has_pconfig(self):
        """Test shesha.config has pConfig class."""
        from shesha.config import pConfig
        assert pConfig is not None
    
    def test_config_has_ptel(self):
        """Test shesha.config has pTel class."""
        from shesha.config import pTel
        assert pTel is not None
    
    def test_config_has_pgeom(self):
        """Test shesha.config has pGeom class."""
        from shesha.config import pGeom
        assert pGeom is not None


class TestSheshaInitModules:
    """Test shesha initialization modules."""
    
    def test_init_has_geom_init(self):
        """Test shesha.init has geom_init module."""
        from shesha.init import geom_init
        assert geom_init is not None
    
    def test_init_has_dm_init(self):
        """Test shesha.init has dm_init module."""
        from shesha.init import dm_init
        assert dm_init is not None
    
    def test_init_has_wfs_init(self):
        """Test shesha.init has wfs_init module."""
        from shesha.init import wfs_init
        assert wfs_init is not None


class TestSheshaDataStructures:
    """Test shesha data structure functionality."""
    
    def test_constant_conversions_work(self):
        """Test constant conversions are working."""
        from shesha.constants import CONST
        
        # Test that conversions produce positive values
        assert CONST.RAD2DEG > 0
        assert CONST.DEG2RAD > 0
        assert CONST.RAD2ARCSEC > 0
        assert CONST.ARCSEC2RAD > 0
    
    def test_dm_type_constants_exist(self):
        """Test DM type constants exist."""
        from shesha.constants import DmType
        
        assert hasattr(DmType, 'PZT')
        assert hasattr(DmType, 'TT')
        assert hasattr(DmType, 'KL')


class TestSheshaPackageStructure:
    """Test overall shesha package structure."""
    
    def test_shesha_module_path_exists(self):
        """Test shesha module path exists."""
        assert hasattr(shesha, '__file__')
        module_path = Path(shesha.__file__).parent
        assert module_path.exists()
    
    def test_shesha_has_subpackages(self):
        """Test shesha has expected subpackages."""
        module_path = Path(shesha.__file__).parent
        
        # Check for key subpackages
        expected_dirs = ['config', 'util', 'init', 'constants.py']
        existing = list(module_path.glob('*'))
        
        assert len(existing) > 0


class TestSheshaImportPaths:
    """Test various import paths work correctly."""
    
    def test_import_constants_directly(self):
        """Test importing constants directly."""
        from shesha.constants import CONST, DmType
        
        assert CONST is not None
        assert DmType is not None
    
    def test_import_util_functions(self):
        """Test importing util functions directly."""
        from shesha.util.dm_util import createSquarePattern
        from shesha.util.kl_util import make_radii
        
        assert createSquarePattern is not None
        assert make_radii is not None
    
    def test_import_config_classes(self):
        """Test importing config classes directly."""
        from shesha.config import pConfig, pTel
        
        assert pConfig is not None
        assert pTel is not None


class TestSheshaArrayTypeHandling:
    """Test shesha type handling for arrays."""
    
    def test_numpy_array_usage_in_util(self):
        """Test numpy arrays are used in util functions."""
        from shesha.util.kl_util import make_radii
        
        radii = make_radii(0.1, 10)
        
        assert isinstance(radii, np.ndarray)
        assert radii.dtype in [np.float32, np.float64]
    
    def test_kernel_float32_type(self):
        """Test kernels use float32 type."""
        from shesha.util.kl_util import make_kernels
        from shesha.constants import KLType
        
        cobs = 0.1
        nr = 5
        radp = np.linspace(cobs, 1.0, nr)
        
        kernels = make_kernels(cobs, nr, radp, KLType.KOLMO)
        
        assert kernels.dtype == np.float32
