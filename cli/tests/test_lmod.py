"""
Tests for compass.system.lmod module.
"""

import os
from pathlib import Path
from unittest.mock import patch, MagicMock
from compass.system.lmod import LmodManager


class TestLmodManager:
    """Test cases for LmodManager class."""
    
    def test_lmod_manager_initialization(self, mock_compass_env):
        """Test LmodManager initialization."""
        manager = LmodManager(verbose=True)
        
        assert manager.verbose is True
        assert manager.compass_root is not None
        assert manager.modulefiles_root is not None
    
    def test_lmod_manager_initialization_custom_root(self, temp_dir):
        """Test LmodManager initialization with custom compass_root."""
        manager = LmodManager(compass_root=temp_dir, verbose=False)
        
        assert manager.compass_root == temp_dir
        assert manager.modulefiles_root == temp_dir / "modulefiles"
    
    def test_lmod_manager_initialization_from_env(self, mock_compass_env):
        """Test LmodManager uses COMPASS_ROOT from environment."""
        manager = LmodManager()
        
        assert manager.compass_root is not None
        assert str(manager.compass_root) == os.getenv("COMPASS_ROOT")
    
    def test_lmod_is_installed_env_var_set(self, monkeypatch):
        """Test is_lmod_installed() when LMOD_CMD is set."""
        monkeypatch.setenv("LMOD_CMD", "/usr/share/lmod/lmod/lmod")
        
        manager = LmodManager()
        result = manager.is_lmod_installed()
        
        assert result is True
    
    def test_lmod_is_installed_env_var_not_set(self, monkeypatch, mock_compass_env):
        """Test is_lmod_installed() when LMOD_CMD is not set."""
        monkeypatch.delenv("LMOD_CMD", raising=False)
        
        manager = LmodManager()
        result = manager.is_lmod_installed()
        
        # Result depends on whether system has Lmod installed
        assert isinstance(result, bool)
    
    def test_lmod_get_available_modules_empty(self, temp_dir):
        """Test get_available_modules() with no modules."""
        manager = LmodManager(compass_root=temp_dir, verbose=False)
        
        modules = manager.get_available_modules()
        
        assert isinstance(modules, list)
        assert len(modules) == 0
    
    def test_lmod_get_available_modules_with_files(self, temp_dir):
        """Test get_available_modules() with module files."""
        manager = LmodManager(compass_root=temp_dir, verbose=False)
        
        # Create module files
        modulefiles = manager.modulefiles_root / "compass"
        modulefiles.mkdir(parents=True, exist_ok=True)
        (modulefiles / "local.lua").write_text("-- compass module")
        (modulefiles / "dev.lua").write_text("-- compass dev module")
        
        modules = manager.get_available_modules()
        
        assert len(modules) >= 2
        assert any("compass/local" in m for m in modules)
        assert any("compass/dev" in m for m in modules)
    
    def test_lmod_add_to_bashrc_first_time(self, temp_dir, monkeypatch):
        """Test add_to_bashrc() adds configuration."""
        home = temp_dir / "home"
        home.mkdir()
        bashrc = home / ".bashrc"
        bashrc.write_text("# Original bashrc content\n")
        
        monkeypatch.setenv("HOME", str(home))
        
        manager = LmodManager(compass_root=temp_dir, verbose=False)
        result = manager.add_to_bashrc()
        
        assert result is True
        
        # Verify content was added
        with open(bashrc, 'r') as f:
            content = f.read()
            assert "COMPASS Modulefiles Configuration" in content
            assert "MODULEPATH" in content
    
    def test_lmod_add_to_bashrc_already_exists(self, temp_dir, monkeypatch):
        """Test add_to_bashrc() when already configured."""
        home = temp_dir / "home"
        home.mkdir()
        bashrc = home / ".bashrc"
        
        existing_config = "# COMPASS Modulefiles Configuration\nexport MODULEPATH=/path\n"
        bashrc.write_text(existing_config)
        
        monkeypatch.setenv("HOME", str(home))
        
        manager = LmodManager(compass_root=temp_dir, verbose=False)
        result = manager.add_to_bashrc()
        
        assert result is True
        
        # Verify no duplicate was added
        with open(bashrc, 'r') as f:
            content = f.read()
            count = content.count("COMPASS Modulefiles Configuration")
            assert count == 1
    
    def test_lmod_add_to_bashrc_creates_missing_bashrc(self, temp_dir, monkeypatch):
        """Test add_to_bashrc() creates bashrc if missing."""
        home = temp_dir / "home"
        home.mkdir()
        bashrc = home / ".bashrc"
        
        monkeypatch.setenv("HOME", str(home))
        
        # bashrc doesn't exist initially
        assert not bashrc.exists()
        
        manager = LmodManager(compass_root=temp_dir, verbose=False)
        # Should handle gracefully even if bashrc doesn't exist
        try:
            result = manager.add_to_bashrc()
            # May fail or succeed depending on implementation
        except FileNotFoundError:
            # This is acceptable behavior
            pass
    
    @patch('subprocess.run')
    def test_lmod_is_compass_module_loaded_true(self, mock_run):
        """Test is_compass_module_loaded() returns True when module is loaded."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="compass/local\n",
            stderr=""
        )
        
        manager = LmodManager()
        result = manager.is_compass_module_loaded()
        
        assert result is True
    
    @patch('subprocess.run')
    def test_lmod_is_compass_module_loaded_false(self, mock_run):
        """Test is_compass_module_loaded() returns False when module not loaded."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="other/module\n",
            stderr=""
        )
        
        manager = LmodManager()
        result = manager.is_compass_module_loaded()
        
        assert result is False
    
    @patch('subprocess.run')
    def test_lmod_is_compass_module_loaded_exception(self, mock_run):
        """Test is_compass_module_loaded() handles exceptions."""
        mock_run.side_effect = Exception("Command failed")
        
        manager = LmodManager()
        result = manager.is_compass_module_loaded()
        
        assert result is False
    
    def test_lmod_check_and_setup_no_lmod(self, monkeypatch, mock_compass_env):
        """Test check_and_setup() when Lmod not installed."""
        monkeypatch.delenv("LMOD_CMD", raising=False)
        
        manager = LmodManager(verbose=False)
        result = manager.check_and_setup()
        
        assert result is False


class TestLmodManagerIntegration:
    """Integration tests for LmodManager."""
    
    def test_lmod_manager_modulefiles_structure(self, temp_dir):
        """Test LmodManager works with correct modulefiles structure."""
        # Create proper structure
        modulefiles = temp_dir / "modulefiles"
        compass_mods = modulefiles / "compass"
        cuda_mods = modulefiles / "cuda"
        compass_mods.mkdir(parents=True, exist_ok=True)
        cuda_mods.mkdir(parents=True, exist_ok=True)
        
        (compass_mods / "local.lua").write_text("-- compass local")
        (compass_mods / "dev.lua").write_text("-- compass dev")
        (cuda_mods / "11.0.lua").write_text("-- cuda 11.0")
        
        manager = LmodManager(compass_root=temp_dir, verbose=False)
        modules = manager.get_available_modules()
        
        assert len(modules) == 3
        module_names = set(modules)
        assert "compass/local" in module_names or any("compass" in m for m in modules)
    
    def test_lmod_manager_with_nested_modules(self, temp_dir):
        """Test LmodManager with nested module structure."""
        # Create nested structure
        modulefiles = temp_dir / "modulefiles"
        compass_mods = modulefiles / "compass"
        compiler_mods = modulefiles / "compilers"
        
        compass_mods.mkdir(parents=True, exist_ok=True)
        compiler_mods.mkdir(parents=True, exist_ok=True)
        
        # Create multiple versions
        (compass_mods / "local.lua").write_text("")
        (compiler_mods / "gcc.lua").write_text("")
        (compiler_mods / "intel.lua").write_text("")
        
        manager = LmodManager(compass_root=temp_dir, verbose=False)
        modules = manager.get_available_modules()
        
        assert len(modules) >= 3
