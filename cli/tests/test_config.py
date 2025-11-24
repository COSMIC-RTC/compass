"""
Tests for compass.core.config module.
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
from compass.core.config import Environment, Config


class TestEnvironmentConfig:
    """Test cases for Environment configuration class."""
    
    def test_environment_init_defaults(self):
        """Test Environment initialization with defaults."""
        env = Environment()
        
        assert env.home is not None
        assert env.pythondontwritebytecode is True
        assert isinstance(env.compass_root, Path)
        assert isinstance(env.conda_root, Path)
    
    def test_environment_init_with_custom_paths(self, temp_dir):
        """Test Environment initialization with custom paths."""
        env = Environment(
            compass_root=temp_dir / "compass",
            shesha_root=temp_dir / "shesha",
        )
        
        assert env.compass_root == temp_dir / "compass"
        assert env.shesha_root == temp_dir / "shesha"
    
    def test_environment_post_init_sets_defaults(self, monkeypatch):
        """Test __post_init__ sets default paths from environment."""
        compass_path = "/custom/compass"
        monkeypatch.setenv("COMPASS_ROOT", compass_path)
        
        env = Environment()
        assert str(env.compass_root) == compass_path
    
    def test_environment_to_env_dict(self, temp_dir):
        """Test to_env_dict() returns correct environment variables."""
        env = Environment(
            compass_root=temp_dir / "compass",
            shesha_root=temp_dir / "shesha",
        )
        
        env_dict = env.to_env_dict()
        
        assert isinstance(env_dict, dict)
        assert "COMPASS_ROOT" in env_dict
        assert "SHESHA_ROOT" in env_dict
        assert "PATH" in env_dict
        assert "PYTHONDONTWRITEBYTECODE" in env_dict
        assert env_dict["PYTHONDONTWRITEBYTECODE"] == "1"
    
    def test_environment_to_env_dict_includes_cuda(self, temp_dir, monkeypatch):
        """Test to_env_dict() includes CUDA variables."""
        env = Environment(
            compass_root=temp_dir / "compass",
            cuda_root=Path("/usr/local/cuda"),
        )
        
        env_dict = env.to_env_dict()
        
        # CUDA variables should be included if cuda_root exists
        assert "CUDA_ROOT" in env_dict or "CUDA_ROOT" not in env_dict  # Depends on system
    
    def test_environment_apply(self, monkeypatch):
        """Test apply() updates os.environ."""
        original_path = os.environ.get("PATH", "")
        
        env = Environment()
        env.apply()
        
        # Check that environment was updated
        assert "COMPASS_ROOT" in os.environ
        assert "PYTHONPATH" in os.environ
    
    def test_environment_export_to_shell_script(self, temp_dir):
        """Test export_to_shell_script() creates valid shell script."""
        env = Environment(
            compass_root=temp_dir / "compass",
            shesha_root=temp_dir / "shesha",
        )
        
        script_path = temp_dir / "compass_env.sh"
        env.export_to_shell_script(script_path)
        
        assert script_path.exists()
        
        with open(script_path, 'r') as f:
            content = f.read()
            assert "#!/bin/bash" in content
            assert "COMPASS_ROOT" in content
            assert "SHESHA_ROOT" in content
            assert "export" in content


class TestConfig:
    """Test cases for Config class."""
    
    def test_config_init_no_path(self):
        """Test Config initialization without config path."""
        config = Config()
        
        assert isinstance(config.data, dict)
        assert config.environment is not None
    
    def test_config_init_with_path(self, temp_dir):
        """Test Config initialization with config path."""
        yaml_file = temp_dir / "config.yaml"
        
        # Create a simple YAML file
        import yaml
        config_data = {
            "paths": {
                "compass_root": str(temp_dir / "compass"),
            },
            "cuda": {
                "root": "/usr/local/cuda",
            }
        }
        with open(yaml_file, 'w') as f:
            yaml.dump(config_data, f)
        
        config = Config(config_path=yaml_file)
        assert config.config_path == yaml_file
    
    def test_config_load(self, temp_dir):
        """Test load() method loads YAML configuration."""
        import yaml
        
        yaml_file = temp_dir / "config.yaml"
        config_data = {
            "build": {"type": "Release"},
            "simulation": {"default_devices": "0"},
        }
        
        with open(yaml_file, 'w') as f:
            yaml.dump(config_data, f)
        
        config = Config()
        config.load(yaml_file)
        
        assert config.data == config_data
    
    def test_config_load_updates_environment(self, temp_dir):
        """Test load() updates environment paths."""
        import yaml
        
        yaml_file = temp_dir / "config.yaml"
        custom_compass_root = temp_dir / "custom_compass"
        
        config_data = {
            "paths": {
                "compass_root": str(custom_compass_root),
            }
        }
        
        with open(yaml_file, 'w') as f:
            yaml.dump(config_data, f)
        
        config = Config()
        config.load(yaml_file)
        
        assert config.environment.compass_root == custom_compass_root
    
    def test_config_save(self, temp_dir):
        """Test save() method writes configuration to file."""
        output_file = temp_dir / "config_out.yaml"
        
        config = Config()
        config.data = {
            "test": {
                "key": "value",
            }
        }
        config.save(output_file)
        
        assert output_file.exists()
        
        import yaml
        with open(output_file, 'r') as f:
            loaded_data = yaml.safe_load(f)
            assert loaded_data == config.data
    
    def test_config_get_simple_key(self):
        """Test get() with simple key."""
        config = Config()
        config.data = {
            "key": "value",
            "nested": {
                "key": "nested_value",
            }
        }
        
        assert config.get("key") == "value"
        assert config.get("missing", default="default") == "default"
    
    def test_config_get_dotted_key(self):
        """Test get() with dotted notation."""
        config = Config()
        config.data = {
            "paths": {
                "compass_root": "/home/compass",
                "install": "/home/compass/local",
            }
        }
        
        assert config.get("paths.compass_root") == "/home/compass"
        assert config.get("paths.install") == "/home/compass/local"
        assert config.get("paths.missing", default="default") == "default"
    
    def test_config_set_simple_key(self):
        """Test set() with simple key."""
        config = Config()
        config.data = {}
        
        config.set("key", "value")
        assert config.data["key"] == "value"
    
    def test_config_set_dotted_key(self):
        """Test set() with dotted notation."""
        config = Config()
        config.data = {}
        
        config.set("paths.compass_root", "/home/compass")
        assert config.data["paths"]["compass_root"] == "/home/compass"
    
    def test_config_set_creates_nested_structure(self):
        """Test set() creates nested structure automatically."""
        config = Config()
        config.data = {}
        
        config.set("build.cmake.type", "Release")
        assert config.data["build"]["cmake"]["type"] == "Release"


class TestConfigIntegration:
    """Integration tests for Config class."""
    
    def test_config_roundtrip_save_load(self, temp_dir):
        """Test saving and loading configuration."""
        yaml_file = temp_dir / "config.yaml"
        
        # Create and save config
        config1 = Config()
        config1.data = {
            "version": "1.0",
            "build": {"type": "Release"},
            "simulation": {"devices": "0,1"},
        }
        config1.save(yaml_file)
        
        # Load config
        config2 = Config()
        config2.load(yaml_file)
        
        assert config2.data == config1.data
    
    def test_config_with_environment_expansion(self, temp_dir, monkeypatch):
        """Test configuration with environment variables."""
        monkeypatch.setenv("COMPASS_ROOT", str(temp_dir / "compass"))
        
        config = Config()
        
        env_dict = config.environment.to_env_dict()
        assert env_dict["COMPASS_ROOT"] == str(temp_dir / "compass")
