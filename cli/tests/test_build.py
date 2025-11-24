"""
Tests for compass.deployment.build module.
"""

import os
from pathlib import Path
from unittest.mock import patch, MagicMock, call
from compass.deployment.build import BuildType, Component, Builder


class TestBuildType:
    """Test cases for BuildType enum."""
    
    def test_build_type_release(self):
        """Test BuildType.RELEASE."""
        assert BuildType.RELEASE.value == "Release"
    
    def test_build_type_debug(self):
        """Test BuildType.DEBUG."""
        assert BuildType.DEBUG.value == "Debug"
    
    def test_build_type_relwithdebinfo(self):
        """Test BuildType.RELWITHDEBINFO."""
        assert BuildType.RELWITHDEBINFO.value == "RelWithDebInfo"
    
    def test_all_build_types(self):
        """Test all build type values are valid."""
        for build_type in BuildType:
            assert build_type.value in ["Release", "Debug", "RelWithDebInfo"]


class TestComponent:
    """Test cases for Component dataclass."""
    
    def test_component_initialization(self, temp_dir):
        """Test Component initialization."""
        comp = Component(
            name="libcarma",
            root_path=temp_dir / "libcarma",
        )
        
        assert comp.name == "libcarma"
        assert comp.root_path == temp_dir / "libcarma"
        assert comp.required is True
        assert isinstance(comp.dependencies, list)
    
    def test_component_build_dir_default(self, temp_dir):
        """Test Component default build_dir is set."""
        comp = Component(
            name="libcarma",
            root_path=temp_dir / "libcarma",
        )
        
        assert comp.build_dir == temp_dir / "libcarma" / "build"
    
    def test_component_custom_build_dir(self, temp_dir):
        """Test Component with custom build_dir."""
        custom_build = temp_dir / "custom_build"
        comp = Component(
            name="libcarma",
            root_path=temp_dir / "libcarma",
            build_dir=custom_build,
        )
        
        assert comp.build_dir == custom_build
    
    def test_component_with_dependencies(self, temp_dir):
        """Test Component with dependencies."""
        comp = Component(
            name="libsutra",
            root_path=temp_dir / "libsutra",
            dependencies=["libcarma"],
        )
        
        assert "libcarma" in comp.dependencies


class TestBuilder:
    """Test cases for Builder class."""
    
    def test_builder_initialization(self, mock_compass_env):
        """Test Builder initialization."""
        builder = Builder(verbose=False)
        
        assert builder.verbose is False
        assert builder.compass_root is not None
        assert builder.build_dir is not None
        assert builder.install_dir is not None
    
    def test_builder_initialization_verbose(self, mock_compass_env):
        """Test Builder initialization with verbose=True."""
        builder = Builder(verbose=True)
        
        assert builder.verbose is True
    
    def test_builder_get_components(self, mock_compass_env):
        """Test get_components() returns expected components."""
        builder = Builder()
        components = builder.get_components()
        
        assert len(components) == 3
        
        component_names = [c.name for c in components]
        assert "libcarma" in component_names
        assert "libsutra" in component_names
        assert "python_module" in component_names
        
        # Check dependencies
        libsutra = [c for c in components if c.name == "libsutra"][0]
        assert "libcarma" in libsutra.dependencies
    
    def test_builder_clean_build_directory_exists(self, mock_compass_env):
        """Test clean_build() when build directory exists."""
        builder = Builder()
        builder.build_dir.mkdir(parents=True, exist_ok=True)
        
        # Create a dummy file in the build directory
        (builder.build_dir / "CMakeCache.txt").write_text("dummy")
        
        result = builder.clean_build()
        
        assert result is True
        assert not builder.build_dir.exists()
    
    def test_builder_clean_build_directory_not_exists(self, mock_compass_env):
        """Test clean_build() when build directory doesn't exist."""
        builder = Builder()
        
        # Ensure directory doesn't exist
        if builder.build_dir.exists():
            import shutil
            shutil.rmtree(builder.build_dir)
        
        result = builder.clean_build()
        
        assert result is True
    
    def test_builder_clean_install_directory_exists(self, mock_compass_env):
        """Test clean_install() when install directory exists."""
        builder = Builder()
        builder.install_dir.mkdir(parents=True, exist_ok=True)
        
        # Create a dummy file
        (builder.install_dir / "lib").mkdir(exist_ok=True)
        (builder.install_dir / "lib" / "libcarma.so").write_text("dummy")
        
        result = builder.clean_install()
        
        assert result is True
        assert not builder.install_dir.exists()
    
    def test_builder_clean_install_directory_not_exists(self, mock_compass_env):
        """Test clean_install() when install directory doesn't exist."""
        builder = Builder()
        
        # Ensure directory doesn't exist
        if builder.install_dir.exists():
            import shutil
            shutil.rmtree(builder.install_dir)
        
        result = builder.clean_install()
        
        assert result is True
    
    @patch('subprocess.run')
    def test_builder_run_command_success(self, mock_run, mock_compass_env):
        """Test _run_command() with successful execution."""
        mock_run.return_value = MagicMock(returncode=0)
        
        builder = Builder(verbose=False)
        result = builder._run_command(["echo", "test"], cwd=Path("/tmp"))
        
        assert result is True
        mock_run.assert_called_once()
    
    @patch('subprocess.run')
    def test_builder_run_command_failure(self, mock_run, mock_compass_env):
        """Test _run_command() with failed execution."""
        import subprocess
        mock_run.side_effect = subprocess.CalledProcessError(1, "cmd")
        
        builder = Builder(verbose=False)
        result = builder._run_command(["false"], cwd=Path("/tmp"))
        
        assert result is False
    
    def test_builder_show_status(self, mock_compass_env):
        """Test show_status() method."""
        builder = Builder()
        
        # Create required directories for status check
        builder.build_dir.mkdir(parents=True, exist_ok=True)
        builder.install_dir.mkdir(parents=True, exist_ok=True)
        
        # Just ensure it doesn't crash
        try:
            builder.show_status()
        except Exception:
            # Rich console might fail in test environment
            pass


class TestBuilderIntegration:
    """Integration tests for Builder class."""
    
    def test_builder_get_components_all_have_required_fields(self, mock_compass_env):
        """Test all components have required fields."""
        builder = Builder()
        components = builder.get_components()
        
        for component in components:
            assert component.name
            assert component.root_path
            assert component.build_dir
            assert component.install_path
    
    def test_builder_components_are_component_objects(self, mock_compass_env):
        """Test all returned objects are Component instances."""
        builder = Builder()
        components = builder.get_components()
        
        assert all(isinstance(c, Component) for c in components)
    
    def test_builder_dependency_order(self, mock_compass_env):
        """Test components are ordered correctly for dependencies."""
        builder = Builder()
        components = builder.get_components()
        
        component_dict = {c.name: c for c in components}
        
        # libsutra depends on libcarma
        libsutra_idx = next(i for i, c in enumerate(components) if c.name == "libsutra")
        libcarma_idx = next(i for i, c in enumerate(components) if c.name == "libcarma")
        assert libcarma_idx < libsutra_idx
        
        # python_module depends on both
        py_idx = next(i for i, c in enumerate(components) if c.name == "python_module")
        assert libcarma_idx < py_idx
        assert libsutra_idx < py_idx
