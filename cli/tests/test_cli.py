"""
Tests for compass.cli module - CLI command handlers and interface.
"""

import sys
from pathlib import Path
from unittest.mock import patch, MagicMock
from compass import __version__


class TestCliBasics:
    """Basic test cases for CLI functionality."""
    
    def test_cli_version_available(self):
        """Test that CLI version is available."""
        assert __version__ is not None
        assert isinstance(__version__, str)
        # Version should be in a common format like X.Y.Z
        parts = __version__.split('.')
        assert len(parts) >= 2


class TestCliImports:
    """Test that CLI modules can be imported."""
    
    def test_import_builder(self):
        """Test Builder can be imported from CLI."""
        from compass.deployment.build import Builder
        assert Builder is not None
    
    def test_import_simulation_runner(self):
        """Test SimulationRunner can be imported from CLI."""
        from compass.simulation.runner import SimulationRunner
        assert SimulationRunner is not None
    
    def test_import_environment_checker(self):
        """Test EnvironmentChecker can be imported from CLI."""
        from compass.core.environment import EnvironmentChecker
        assert EnvironmentChecker is not None
    
    def test_import_lmod_manager(self):
        """Test LmodManager can be imported from CLI."""
        from compass.system.lmod import LmodManager
        assert LmodManager is not None
    
    def test_import_config(self):
        """Test Config can be imported from CLI."""
        from compass.core.config import Config, Environment
        assert Config is not None
        assert Environment is not None


class TestCliEnvironmentSetup:
    """Test CLI environment setup functionality."""
    
    def test_environment_has_required_attributes(self, mock_compass_env):
        """Test Environment has required attributes."""
        from compass.core.config import Environment
        
        env = Environment()
        
        assert hasattr(env, 'compass_root')
        assert hasattr(env, 'shesha_root')
        assert hasattr(env, 'to_env_dict')
        assert hasattr(env, 'export_to_shell_script')
    
    def test_environment_checker_has_required_methods(self):
        """Test EnvironmentChecker has required methods."""
        from compass.core.environment import EnvironmentChecker
        
        checker = EnvironmentChecker()
        
        assert hasattr(checker, 'get_required_variables')
        assert hasattr(checker, 'check_all')
        assert hasattr(checker, 'check_and_report')
        assert hasattr(checker, 'show_status')


class TestCliBuildModule:
    """Test CLI build module functionality."""
    
    def test_builder_has_required_methods(self, mock_compass_env):
        """Test Builder has required methods."""
        from compass.deployment.build import Builder
        
        builder = Builder()
        
        assert hasattr(builder, 'get_components')
        assert hasattr(builder, 'clean_build')
        assert hasattr(builder, 'clean_install')
        assert hasattr(builder, 'build_all')
        assert hasattr(builder, 'show_status')
    
    def test_builder_components_structure(self, mock_compass_env):
        """Test Builder components have correct structure."""
        from compass.deployment.build import Builder
        
        builder = Builder()
        components = builder.get_components()
        
        for component in components:
            assert hasattr(component, 'name')
            assert hasattr(component, 'root_path')
            assert hasattr(component, 'build_dir')
            assert hasattr(component, 'dependencies')


class TestCliSimulationModule:
    """Test CLI simulation module functionality."""
    
    def test_simulation_runner_has_required_methods(self, mock_compass_env):
        """Test SimulationRunner has required methods."""
        from compass.simulation.runner import SimulationRunner
        
        runner = SimulationRunner()
        
        assert hasattr(runner, 'get_script_path')
        assert hasattr(runner, 'load_config')
        assert hasattr(runner, 'save_config')
        assert hasattr(runner, 'set_default_script')
        assert hasattr(runner, 'show_script_config')
    
    def test_simulation_runner_initialization(self, mock_compass_env):
        """Test SimulationRunner initializes correctly."""
        from compass.simulation.runner import SimulationRunner
        
        runner = SimulationRunner(verbose=True)
        
        assert runner.verbose is True
        assert runner.compass_root is not None
        assert runner.shesha_root is not None
        assert runner.config_file is not None


class TestCliSystemModule:
    """Test CLI system module functionality."""
    
    def test_lmod_manager_has_required_methods(self):
        """Test LmodManager has required methods."""
        from compass.system.lmod import LmodManager
        
        manager = LmodManager()
        
        assert hasattr(manager, 'is_lmod_installed')
        assert hasattr(manager, 'get_available_modules')
        assert hasattr(manager, 'add_to_bashrc')
        assert hasattr(manager, 'is_compass_module_loaded')
        assert hasattr(manager, 'check_and_setup')
    
    def test_lmod_manager_initialization(self):
        """Test LmodManager initializes correctly."""
        from compass.system.lmod import LmodManager
        
        manager = LmodManager(verbose=True)
        
        assert manager.verbose is True
        assert manager.compass_root is not None
        assert manager.modulefiles_root is not None


class TestCliConfigModule:
    """Test CLI config module functionality."""
    
    def test_config_has_required_methods(self):
        """Test Config has required methods."""
        from compass.core.config import Config
        
        config = Config()
        
        assert hasattr(config, 'load')
        assert hasattr(config, 'save')
        assert hasattr(config, 'get')
        assert hasattr(config, 'set')
        assert hasattr(config, 'data')
        assert hasattr(config, 'environment')
    
    def test_config_initialization(self):
        """Test Config initializes correctly."""
        from compass.core.config import Config
        
        config = Config()
        
        assert isinstance(config.data, dict)
        assert config.environment is not None


class TestCliLoggerModule:
    """Test CLI logger module functionality."""
    
    def test_logger_setup_function_available(self):
        """Test setup_logger function is available."""
        from compass.core.logger import setup_logger
        
        logger = setup_logger(name="test")
        assert logger is not None
        assert logger.name == "test"
    
    def test_logger_mixin_available(self):
        """Test LoggerMixin is available."""
        from compass.core.logger import LoggerMixin
        
        class TestClass(LoggerMixin):
            pass
        
        obj = TestClass()
        assert obj.logger is not None


class TestCliVersionAndMetadata:
    """Test CLI version and metadata."""
    
    def test_version_format(self):
        """Test version has correct format."""
        from compass import __version__
        
        # Version should contain numbers and dots
        assert any(c.isdigit() for c in __version__)
        assert '.' in __version__
    
    def test_package_info_available(self):
        """Test package info is available."""
        import compass
        
        assert hasattr(compass, '__version__')


class TestCliModuleStructure:
    """Test overall CLI module structure."""
    
    def test_core_module_structure(self):
        """Test core module has expected submodules."""
        from compass import core
        
        assert hasattr(core, 'config')
        assert hasattr(core, 'environment')
        assert hasattr(core, 'logger')
    
    def test_deployment_module_structure(self):
        """Test deployment module has expected submodules."""
        from compass import deployment
        
        assert hasattr(deployment, 'build')
    
    def test_simulation_module_structure(self):
        """Test simulation module has expected submodules."""
        from compass import simulation
        
        assert hasattr(simulation, 'runner')
    
    def test_system_module_structure(self):
        """Test system module has expected submodules."""
        from compass import system
        
        assert hasattr(system, 'lmod')
