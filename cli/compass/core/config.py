"""
Configuration and Environment Management

This module provides centralized configuration management for COMPASS.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, Optional, Any
from dataclasses import dataclass, field


@dataclass
class Environment:
    """
    Manages environment variables and paths for COMPASS.
    
    This provides a programmatic interface for environment management.
    """
    
    # Base paths
    home: Path = field(default_factory=lambda: Path.home())
    compass_root: Optional[Path] = None
    compass_install_root: Optional[Path] = None
    shesha_root: Optional[Path] = None
    conda_root: Optional[Path] = None
    
    # CUDA paths
    cuda_root: Path = Path("/usr/local/cuda")
    
    # Additional settings
    pythondontwritebytecode: bool = True
    
    def __post_init__(self):
        """Initialize default paths if not provided."""
        if self.compass_root is None:
            self.compass_root = Path(os.getenv("COMPASS_ROOT", self.home / "compass"))
        
        if self.compass_install_root is None:
            self.compass_install_root = Path(
                os.getenv("COMPASS_INSTALL_ROOT", self.compass_root / "local")
            )
        
        if self.shesha_root is None:
            self.shesha_root = Path(
                os.getenv("SHESHA_ROOT", self.compass_root / "shesha")
            )
        
        if self.conda_root is None:
            # Try different conda variants
            conda_default = self.home / "miniconda3"
            if not conda_default.exists():
                conda_default = self.home / "anaconda3"
            if not conda_default.exists():
                conda_default = self.home / "miniforge3"
            self.conda_root = Path(os.getenv("CONDA_ROOT", conda_default))
    
    def to_env_dict(self) -> Dict[str, str]:
        """
        Convert configuration to environment variable dictionary.
        
        Returns:
            Dictionary of environment variable names to values
        """
        env = {}
        
        # Base paths
        env["COMPASS_ROOT"] = str(self.compass_root)
        env["COMPASS_INSTALL_ROOT"] = str(self.compass_install_root)
        env["SHESHA_ROOT"] = str(self.shesha_root)
        env["CONDA_ROOT"] = str(self.conda_root)
        
        # CUDA
        if self.cuda_root.exists():
            env["CUDA_ROOT"] = str(self.cuda_root)
            env["CUDA_INC_PATH"] = str(self.cuda_root / "include")
            env["CUDA_LIB_PATH"] = str(self.cuda_root / "lib")
            env["CUDA_LIB_PATH_64"] = str(self.cuda_root / "lib64")
        
        # Python settings
        if self.pythondontwritebytecode:
            env["PYTHONDONTWRITEBYTECODE"] = "1"
        
        # Build PATH
        path_components = []
        if self.conda_root.exists():
            path_components.append(str(self.conda_root / "bin"))
        if self.cuda_root.exists():
            path_components.append(str(self.cuda_root / "bin"))
        if self.compass_install_root:
            path_components.append(str(self.compass_install_root / "bin"))
        path_components.append(os.environ.get("PATH", "/usr/local/bin:/usr/bin:/bin"))
        env["PATH"] = ":".join(path_components)
        
        # Build LD_LIBRARY_PATH
        ld_lib_components = []
        if self.cuda_root.exists():
            ld_lib_components.extend([
                str(self.cuda_root / "lib64"),
                str(self.cuda_root / "lib"),
            ])
        if self.compass_install_root:
            ld_lib_components.append(str(self.compass_install_root / "lib"))
        ld_lib_components.append(os.environ.get("LD_LIBRARY_PATH", ""))
        env["LD_LIBRARY_PATH"] = ":".join(filter(None, ld_lib_components))
        
        # Build PYTHONPATH
        python_components = []
        if self.shesha_root:
            python_components.append(str(self.shesha_root))
        if self.compass_install_root:
            python_dir = self.compass_install_root / "python"
            if python_dir.exists():
                python_components.append(str(python_dir))
        python_components.append(os.environ.get("PYTHONPATH", ""))
        env["PYTHONPATH"] = ":".join(filter(None, python_components))
        
        # Build PKG_CONFIG_PATH
        pkg_config_components = []
        if self.compass_install_root:
            pkg_config_components.append(str(self.compass_install_root / "lib" / "pkgconfig"))
        pkg_config_components.append(os.environ.get("PKG_CONFIG_PATH", ""))
        env["PKG_CONFIG_PATH"] = ":".join(filter(None, pkg_config_components))
        
        return env
    
    def apply(self):
        """Apply environment variables to current process."""
        env_dict = self.to_env_dict()
        os.environ.update(env_dict)
    
    def export_to_shell_script(self, output_path: Path):
        """
        Export environment to a shell script for sourcing.
        
        Args:
            output_path: Path to write the shell script
        """
        env_dict = self.to_env_dict()
        
        with open(output_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# COMPASS Environment Configuration\n")
            f.write("# Auto-generated by compass CLI\n\n")
            
            for key, value in env_dict.items():
                if key not in ["HOME"]:  # Don't override HOME
                    f.write(f'export {key}="{value}"\n')
            
            f.write("\necho 'COMPASS environment configured:'\n")
            f.write("echo '  COMPASS_ROOT='$COMPASS_ROOT\n")
            f.write("echo '  SHESHA_ROOT='$SHESHA_ROOT\n")
        
        output_path.chmod(0o755)


class Config:
    """
    Centralized configuration management for COMPASS.
    
    Supports loading from YAML files and environment variables.
    """
    
    def __init__(self, config_path: Optional[Path] = None):
        """
        Initialize configuration.
        
        Args:
            config_path: Path to YAML configuration file
        """
        self.config_path = config_path
        self.data: Dict[str, Any] = {}
        self.environment = Environment()
        
        if config_path and config_path.exists():
            self.load(config_path)
    
    def load(self, config_path: Path):
        """Load configuration from YAML file."""
        with open(config_path, 'r') as f:
            self.data = yaml.safe_load(f)
        
        # Update environment from config
        if 'paths' in self.data:
            paths = self.data['paths']
            if 'compass_root' in paths:
                self.environment.compass_root = Path(paths['compass_root']).expanduser()
            if 'shesha_root' in paths:
                self.environment.shesha_root = Path(paths['shesha_root']).expanduser()
            if 'conda_root' in paths:
                self.environment.conda_root = Path(paths['conda_root']).expanduser()
        
        if 'cuda' in self.data and 'root' in self.data['cuda']:
            self.environment.cuda_root = Path(self.data['cuda']['root'])
    
    def save(self, config_path: Optional[Path] = None):
        """Save configuration to YAML file."""
        output_path = config_path or self.config_path
        if not output_path:
            raise ValueError("No config path specified")
        
        with open(output_path, 'w') as f:
            yaml.dump(self.data, f, default_flow_style=False, sort_keys=False)
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value by key (supports dot notation)."""
        keys = key.split('.')
        value = self.data
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value
    
    def set(self, key: str, value: Any):
        """Set configuration value by key (supports dot notation)."""
        keys = key.split('.')
        data = self.data
        
        for k in keys[:-1]:
            if k not in data:
                data[k] = {}
            data = data[k]
        
        data[keys[-1]] = value
