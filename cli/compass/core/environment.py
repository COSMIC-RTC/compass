"""
Environment validation module.

Checks that the COMPASS environment is properly configured.
"""

import os
import json
import subprocess
from typing import List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
from rich.console import Console
from rich.table import Table


@dataclass
class EnvVariable:
    """Represents a required environment variable."""
    
    name: str
    description: str
    required: bool = True
    
    def is_set(self) -> bool:
        """Check if the environment variable is set."""
        value = os.environ.get(self.name)
        return value is not None and value.strip() != ""
    
    def get_value(self) -> Optional[str]:
        """Get the environment variable value."""
        return os.environ.get(self.name)


class EnvironmentChecker:
    """
    Checks that the COMPASS environment is properly configured.
    
    Verifies that required environment variables are set and paths exist.
    """
    
    def __init__(self, verbose: bool = True):
        """
        Initialize the environment checker.
        
        Args:
            verbose: Enable verbose output (default: True)
        """
        self.verbose = verbose
    
    def get_required_variables(self) -> List[EnvVariable]:
        """
        Get list of required environment variables.
        
        Returns:
            List of required environment variables
        """
        return [
            EnvVariable("COMPASS_ROOT", "Root COMPASS directory"),
            EnvVariable("SHESHA_ROOT", "Shesha Python package path"),
            EnvVariable("COMPASS_INSTALL_ROOT", "COMPASS installation path", required=False),
        ]
    
    def show_status(self):
        """
        Print the status of required environment variables.
        """
        console = Console()
        required_vars = self.get_required_variables()
        console.print("[bold cyan]" + "=" * 60 + "[/bold cyan]")
        console.print("[bold cyan]COMPASS Environment Configuration[/bold cyan]")
        console.print("[bold cyan]" + "=" * 60 + "[/bold cyan]\n")
        table = Table(title="COMPASS Environment Configuration")
        table.add_column("Variable", style="cyan")
        table.add_column("Value", style="green")
        table.add_column("Status", style="white")
        
        for var in required_vars:
            value = var.get_value() or "[red]Not set[/red]"
            status = "✓" if var.is_set() else ("✗" if var.required else "⚠")
            table.add_row(var.name, value, status)
        
        console.print(table)
        
        # Check for Python packages
        console.print(f"\n[bold]Python Environment:[/bold]")
        self._check_python_packages()
    
    def _check_python_packages(self):
        """Check if required Python packages are available."""
        console = Console()
        
        # Try to load requirements from compass requirements.txt
        required_packages = self._get_required_packages()
        
        # Get installed packages from conda/pip with versions
        installed_packages = self._get_installed_packages()
        
        for package in required_packages:
            # Normalize package name for comparison (lowercase, handle variants)
            package_lower = package.lower()
            if package_lower in installed_packages:
                version = installed_packages[package_lower]
                console.print(f"  [green]✓ {package}[/green] [dim]({version})[/dim]")
            else:
                console.print(f"  [red]✗ {package} (not installed)[/red]")
    
    def _get_installed_packages(self) -> dict:
        """
        Get dict of installed packages from conda list and/or pip list.
        
        Returns:
            Dict mapping package names (lowercase) to versions
        """
        installed = {}
        
        # Try conda/mamba list first
        for cmd in ['conda', 'mamba']:
            try:
                result = subprocess.run(
                    [cmd, 'list', '--json'],
                    capture_output=True,
                    text=True,
                    timeout=10,
                    check=False
                )
                if result.returncode == 0:
                    packages = json.loads(result.stdout)
                    for pkg in packages:
                        if isinstance(pkg, dict) and 'name' in pkg:
                            name = pkg['name'].lower()
                            version = pkg.get('version', 'unknown')
                            installed[name] = version
                    break  # Use first available (conda or mamba)
            except (subprocess.TimeoutExpired, FileNotFoundError, json.JSONDecodeError):
                continue
        
        # Also check pip list (may have additional packages or override versions)
        try:
            result = subprocess.run(
                ['pip', 'list', '--format=json'],
                capture_output=True,
                text=True,
                timeout=10,
                check=False
            )
            if result.returncode == 0:
                packages = json.loads(result.stdout)
                for pkg in packages:
                    if isinstance(pkg, dict) and 'name' in pkg:
                        name = pkg['name'].lower()
                        version = pkg.get('version', 'unknown')
                        # Only add if not already present (prefer conda version)
                        if name not in installed:
                            installed[name] = version
        except (subprocess.TimeoutExpired, FileNotFoundError, json.JSONDecodeError):
            pass
        
        return installed
    
    def _get_required_packages(self) -> List[str]:
        """
        Get list of required Python packages from requirements.txt or compass-env.yml.
        
        Returns:
            List of package names
        """
        compass_root = os.getenv("COMPASS_ROOT")
        
        # Try requirements.txt first
        if compass_root:
            requirements_file = Path(compass_root) / "cli/requirements.txt"
            if requirements_file.exists():
                try:
                    packages = []
                    with open(requirements_file, 'r') as f:
                        for line in f:
                            line = line.strip()
                            # Skip empty lines and comments
                            if line and not line.startswith('#'):
                                # Extract package name (remove version specifiers)
                                pkg = line.split('>=')[0].split('==')[0].split('<')[0].split('>')[0].strip()
                                if pkg:
                                    packages.append(pkg)
                    if packages:
                        return packages
                except Exception:
                    pass
        
        # Fallback to essential packages
        return ["numpy", "scipy", "astropy", "pyyaml", "matplotlib"]
    
    def check_all(self) -> Tuple[List[EnvVariable], List[str]]:
        """
        Check all required environment variables.
        
        Returns:
            Tuple of (missing variables, warnings)
        """
        required = self.get_required_variables()
        missing = []
        warnings = []
        
        # Check each required variable
        for var in required:
            if not var.is_set():
                if var.required:
                    missing.append(var)
            elif self.verbose:
                print(f"✓ {var.name} = {var.get_value()}")
        
        return missing, warnings
    
    def check_and_report(self) -> bool:
        """
        Check environment and print a report.
        
        Returns:
            True if all required variables are set, False otherwise
        """
        console = Console()
        missing, warnings = self.check_all()
        
        if warnings:
            for warning in warnings:
                console.print(f"[yellow]⚠[/yellow] {warning}")
        
        if missing:
            console.print("\n[bold red]✗ Missing required environment variables:[/bold red]")
            for var in missing:
                console.print(f"  • {var.name}: {var.description}")
            
            console.print("\n[yellow]To fix this:[/yellow]")
            console.print("  1. Set environment variables: [cyan]source compass_env.sh[/cyan]")
            console.print("  2. Or run: [cyan]compass config export[/cyan] to generate the script")
            
            return False
        
        if not warnings:
            console.print("[green]✓ COMPASS environment is properly configured[/green]")
        
        return True
    
    def get_path(self, var_name: str) -> Optional[Path]:
        """
        Get a path from an environment variable.
        
        Args:
            var_name: Name of the environment variable
        
        Returns:
            Path object if variable is set, None otherwise
        """
        value = os.environ.get(var_name)
        if value:
            return Path(value)
        return None


def check_environment(verbose: bool = False) -> bool:
    """
    Convenience function to check the environment.
    
    Args:
        verbose: Enable verbose output
    
    Returns:
        True if environment is properly configured
    """
    checker = EnvironmentChecker(verbose=verbose)
    return checker.check_and_report()
