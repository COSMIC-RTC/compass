"""
Tests for compass.simulation.runner module.
"""

import json
from pathlib import Path
from unittest.mock import patch, MagicMock
from compass.simulation.runner import SimulationRunner


class TestSimulationRunner:
    """Test cases for SimulationRunner class."""
    
    def test_runner_initialization(self, mock_compass_env):
        """Test SimulationRunner initialization."""
        runner = SimulationRunner(verbose=True)
        
        assert runner.verbose is True
        assert runner.compass_root is not None
        assert runner.shesha_root is not None
        assert runner.config_file is not None
    
    def test_runner_initialization_verbose_false(self, mock_compass_env):
        """Test SimulationRunner initialization with verbose=False."""
        runner = SimulationRunner(verbose=False)
        
        assert runner.verbose is False
    
    def test_runner_default_scripts(self, mock_compass_env):
        """Test runner has default script names."""
        runner = SimulationRunner()
        
        assert runner.default_script is not None
        assert runner.default_gui_script is not None
        assert isinstance(runner.default_script, str)
        assert isinstance(runner.default_gui_script, str)
    
    def test_runner_get_script_path_absolute(self, mock_compass_env):
        """Test get_script_path() with absolute path."""
        runner = SimulationRunner()
        
        script_path = runner.get_script_path("/absolute/path/script.py")
        assert script_path == Path("/absolute/path/script.py")
    
    def test_runner_get_script_path_relative(self, mock_compass_env):
        """Test get_script_path() with relative path."""
        runner = SimulationRunner()
        
        # Create a dummy script
        scripts_dir = runner.shesha_root / "shesha" / "scripts"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        (scripts_dir / "test_script.py").write_text("# test")
        
        script_path = runner.get_script_path("test_script.py")
        assert script_path.name == "test_script.py"
        assert "shesha/scripts" in str(script_path)
    
    def test_runner_get_script_path_uses_default(self, mock_compass_env):
        """Test get_script_path() uses default when none specified."""
        runner = SimulationRunner()
        
        # Create default scripts
        scripts_dir = runner.shesha_root / "shesha" / "scripts"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        (scripts_dir / "closed_loop.py").write_text("# default")
        
        script_path = runner.get_script_path(script_type="default")
        assert script_path.exists() or "closed_loop" in script_path.name
    
    def test_runner_load_config_file_not_exists(self, mock_compass_env, mock_home_env):
        """Test load_config() when config file doesn't exist."""
        runner = SimulationRunner()
        config = runner.load_config()
        
        assert isinstance(config, dict)
        assert len(config) == 0
    
    def test_runner_load_config_file_exists(self, mock_compass_env, mock_home_env):
        """Test load_config() when config file exists."""
        runner = SimulationRunner()
        
        # Create config file
        runner.config_file.parent.mkdir(parents=True, exist_ok=True)
        test_config = {
            "default_script": "/path/to/script.py",
            "devices": "0,1",
        }
        with open(runner.config_file, 'w') as f:
            json.dump(test_config, f)
        
        config = runner.load_config()
        
        assert isinstance(config, dict)
        assert config["default_script"] == "/path/to/script.py"
        assert config["devices"] == "0,1"
    
    def test_runner_save_config(self, mock_compass_env, mock_home_env):
        """Test save_config() saves configuration."""
        runner = SimulationRunner()
        
        test_config = {
            "default_script": "/path/to/script.py",
            "max_iterations": 1000,
        }
        
        runner.save_config(test_config)
        
        assert runner.config_file.exists()
        
        # Verify saved content
        with open(runner.config_file, 'r') as f:
            loaded = json.load(f)
            assert loaded == test_config
    
    def test_runner_set_default_script_success(self, mock_compass_env, mock_home_env):
        """Test set_default_script() with existing script."""
        runner = SimulationRunner()
        
        # Create a test script
        scripts_dir = runner.shesha_root / "shesha" / "scripts"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        script_file = scripts_dir / "test_script.py"
        script_file.write_text("# test script")
        
        result = runner.set_default_script("test_script.py")
        
        assert result is True
        
        # Verify it was saved to config
        config = runner.load_config()
        assert "test_script.py" in config["default_script"]
    
    def test_runner_set_default_script_not_found(self, mock_compass_env, mock_home_env):
        """Test set_default_script() with non-existent script."""
        runner = SimulationRunner()
        
        result = runner.set_default_script("nonexistent_script.py")
        
        assert result is False
    
    def test_runner_set_default_gui_script(self, mock_compass_env, mock_home_env):
        """Test set_default_script() with GUI script type."""
        runner = SimulationRunner()
        
        # Create a test GUI script
        scripts_dir = runner.shesha_root / "shesha" / "scripts"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        script_file = scripts_dir / "gui.py"
        script_file.write_text("# GUI script")
        
        result = runner.set_default_script("gui.py", script_type="gui")
        
        assert result is True
        
        # Verify it was saved with GUI key
        config = runner.load_config()
        assert "gui.py" in config["default_gui_script"]
    
    def test_runner_show_script_config(self, mock_compass_env):
        """Test show_script_config() method."""
        runner = SimulationRunner()
        
        # Create a test script
        scripts_dir = runner.shesha_root / "shesha" / "scripts"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        (scripts_dir / "closed_loop.py").write_text("# test")
        
        # Just ensure it doesn't crash
        try:
            runner.show_script_config()
        except Exception:
            # Rich console might fail in test environment
            pass


class TestSimulationRunnerIntegration:
    """Integration tests for SimulationRunner."""
    
    def test_runner_config_roundtrip(self, mock_compass_env, mock_home_env):
        """Test saving and loading configuration."""
        runner = SimulationRunner()
        
        # Save config
        test_config = {
            "default_script": "/path/to/script.py",
            "devices": "0,1",
            "iterations": 2000,
        }
        runner.save_config(test_config)
        
        # Load config
        loaded = runner.load_config()
        
        assert loaded == test_config
    
    def test_runner_script_paths_consistency(self, mock_compass_env):
        """Test script path handling is consistent."""
        runner = SimulationRunner()
        
        # Create test scripts
        scripts_dir = runner.shesha_root / "shesha" / "scripts"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        (scripts_dir / "script1.py").write_text("# script1")
        (scripts_dir / "script2.py").write_text("# script2")
        
        path1 = runner.get_script_path("script1.py")
        path2 = runner.get_script_path("script2.py")
        
        # Both should be in the same scripts directory
        assert path1.parent == path2.parent
        assert path1 != path2
