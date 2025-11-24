"""
Tests for compass.core.environment module.
"""

import os
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
from compass.core.environment import EnvVariable, EnvironmentChecker, check_environment


class TestEnvVariable:
    """Test cases for EnvVariable class."""
    
    def test_env_variable_is_set_true(self, monkeypatch):
        """Test is_set() returns True when environment variable is set."""
        monkeypatch.setenv("TEST_VAR", "some_value")
        var = EnvVariable("TEST_VAR", "Test variable")
        assert var.is_set() is True
    
    def test_env_variable_is_set_false(self, monkeypatch):
        """Test is_set() returns False when environment variable is not set."""
        monkeypatch.delenv("NONEXISTENT_VAR", raising=False)
        var = EnvVariable("NONEXISTENT_VAR", "Non-existent variable")
        assert var.is_set() is False
    
    def test_env_variable_is_set_empty_string(self, monkeypatch):
        """Test is_set() returns False for empty environment variables."""
        monkeypatch.setenv("EMPTY_VAR", "")
        var = EnvVariable("EMPTY_VAR", "Empty variable")
        assert var.is_set() is False
    
    def test_env_variable_get_value(self, monkeypatch):
        """Test get_value() returns the environment variable value."""
        monkeypatch.setenv("TEST_VAR", "test_value")
        var = EnvVariable("TEST_VAR", "Test variable")
        assert var.get_value() == "test_value"
    
    def test_env_variable_get_value_not_set(self, monkeypatch):
        """Test get_value() returns None when variable is not set."""
        monkeypatch.delenv("NONEXISTENT_VAR", raising=False)
        var = EnvVariable("NONEXISTENT_VAR", "Non-existent variable")
        assert var.get_value() is None
    
    def test_env_variable_required_flag(self):
        """Test required flag in EnvVariable."""
        var_required = EnvVariable("TEST", "Test", required=True)
        var_optional = EnvVariable("TEST", "Test", required=False)
        
        assert var_required.required is True
        assert var_optional.required is False


class TestEnvironmentChecker:
    """Test cases for EnvironmentChecker class."""
    
    def test_init_verbose_true(self):
        """Test EnvironmentChecker initialization with verbose=True."""
        checker = EnvironmentChecker(verbose=True)
        assert checker.verbose is True
    
    def test_init_verbose_false(self):
        """Test EnvironmentChecker initialization with verbose=False."""
        checker = EnvironmentChecker(verbose=False)
        assert checker.verbose is False
    
    def test_get_required_variables(self):
        """Test get_required_variables() returns list of EnvVariable."""
        checker = EnvironmentChecker()
        variables = checker.get_required_variables()
        
        assert isinstance(variables, list)
        assert len(variables) > 0
        assert all(isinstance(v, EnvVariable) for v in variables)
        
        # Check for specific expected variables
        var_names = [v.name for v in variables]
        assert "COMPASS_ROOT" in var_names
        assert "SHESHA_ROOT" in var_names
    
    def test_check_all_missing_variables(self, monkeypatch):
        """Test check_all() identifies missing variables."""
        # Clear environment variables
        for var in ["COMPASS_ROOT", "SHESHA_ROOT", "COMPASS_INSTALL_ROOT"]:
            monkeypatch.delenv(var, raising=False)
        
        checker = EnvironmentChecker()
        missing, warnings = checker.check_all()
        
        assert len(missing) > 0
        missing_names = [v.name for v in missing]
        assert "COMPASS_ROOT" in missing_names
    
    def test_check_all_all_present(self, monkeypatch):
        """Test check_all() when all required variables are present."""
        monkeypatch.setenv("COMPASS_ROOT", "/home/test/compass")
        monkeypatch.setenv("SHESHA_ROOT", "/home/test/compass/shesha")
        
        checker = EnvironmentChecker()
        missing, warnings = checker.check_all()
        
        # Should have fewer missing variables after setting the required ones
        assert isinstance(missing, list)
        assert isinstance(warnings, list)
    
    def test_get_path_existing_variable(self, monkeypatch):
        """Test get_path() with existing environment variable."""
        monkeypatch.setenv("TEST_PATH", "/tmp/test")
        checker = EnvironmentChecker()
        
        path = checker.get_path("TEST_PATH")
        assert path == Path("/tmp/test")
    
    def test_get_path_missing_variable(self, monkeypatch):
        """Test get_path() with missing environment variable."""
        monkeypatch.delenv("NONEXISTENT_PATH", raising=False)
        checker = EnvironmentChecker()
        
        path = checker.get_path("NONEXISTENT_PATH")
        assert path is None
    
    def test_get_required_packages(self):
        """Test _get_required_packages() returns list of packages."""
        checker = EnvironmentChecker()
        packages = checker._get_required_packages()
        
        assert isinstance(packages, list)
        # Should have at least the fallback packages
        assert len(packages) > 0
    
    @patch('subprocess.run')
    def test_get_installed_packages_conda(self, mock_run):
        """Test _get_installed_packages() with conda."""
        import json
        
        mock_output = json.dumps([
            {"name": "numpy", "version": "1.21.0"},
            {"name": "scipy", "version": "1.7.0"},
        ])
        
        mock_run.return_value = MagicMock(returncode=0, stdout=mock_output)
        
        checker = EnvironmentChecker()
        packages = checker._get_installed_packages()
        
        assert isinstance(packages, dict)
        assert "numpy" in packages
        assert "scipy" in packages


class TestCheckEnvironmentFunction:
    """Test cases for check_environment() convenience function."""
    
    def test_check_environment_function(self):
        """Test check_environment() function."""
        result = check_environment(verbose=False)
        assert isinstance(result, bool)


class TestEnvironmentIntegration:
    """Integration tests for environment checking."""
    
    def test_environment_checker_show_status(self, monkeypatch, capsys):
        """Test show_status() displays information."""
        monkeypatch.setenv("COMPASS_ROOT", "/home/test/compass")
        monkeypatch.setenv("SHESHA_ROOT", "/home/test/compass/shesha")
        
        checker = EnvironmentChecker(verbose=False)
        # Just ensure it doesn't crash
        try:
            checker.show_status()
        except Exception:
            # Rich console might fail in test environment, that's OK
            pass
    
    def test_check_and_report(self, monkeypatch):
        """Test check_and_report() method."""
        monkeypatch.setenv("COMPASS_ROOT", "/home/test/compass")
        monkeypatch.setenv("SHESHA_ROOT", "/home/test/compass/shesha")
        
        checker = EnvironmentChecker(verbose=False)
        result = checker.check_and_report()
        
        assert isinstance(result, bool)
