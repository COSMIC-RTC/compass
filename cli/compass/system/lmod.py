"""
Lmod (Lua-based Module System) management for COMPASS.

Handles Lmod installation checking, modulefile configuration, and module loading.
"""

import os
import logging
import subprocess
from pathlib import Path
from typing import Optional, List
from rich.console import Console

from compass.core.logger import setup_logger, LoggerMixin


class LmodManager(LoggerMixin):
    """
    Manages Lmod module loading and modulefile setup for COMPASS.
    
    Note: Lmod installation is a prerequisite and must be done via system package manager.
    """
    
    def __init__(self, compass_root: Optional[Path] = None, verbose: bool = True):
        """
        Initialize Lmod manager.
        
        Args:
            compass_root: Root COMPASS directory (default: ~/compass or $COMPASS_ROOT)
            verbose: Enable verbose output (default: True)
        """
        self.verbose = verbose
        self.console = Console()
        
        if compass_root is None:
            compass_root = Path(os.getenv("COMPASS_ROOT", Path.home() / "compass"))
        self.compass_root = compass_root
        
        self.modulefiles_root = self.compass_root / "modulefiles"
        
        # Setup logger via mixin
        if verbose:
            self.logger.setLevel(logging.DEBUG)
    
    def is_lmod_installed(self) -> bool:
        """
        Check if Lmod is installed.
        
        Returns:
            True if Lmod is installed
        """
        # Check if LMOD_CMD environment variable is set (indicates Lmod is active)
        if os.environ.get("LMOD_CMD"):
            return True
        
        # Check for system-wide installation paths
        system_paths = [
            Path("/usr/share/lmod/lmod/init/bash"),
            Path("/usr/share/lmod/lmod/init/profile"),
            Path("/etc/profile.d/z00_lmod.sh"),
            Path("/etc/profile.d/lmod.sh"),
        ]
        
        for path in system_paths:
            if path.exists():
                return True
        
        return False
    
    def get_available_modules(self) -> List[str]:
        """
        Get list of available COMPASS modules.
        
        Returns:
            List of module names
        """
        modules = []
        
        if not self.modulefiles_root.exists():
            return modules
        
        # Look for .lua files in modulefiles subdirectories
        for subdir in self.modulefiles_root.iterdir():
            if subdir.is_dir():
                for module_file in subdir.glob("*.lua"):
                    module_name = f"{subdir.name}/{module_file.stem}"
                    modules.append(module_name)
        
        return modules
    
    def add_to_bashrc(self) -> bool:
        """
        Add modulefile configuration to ~/.bashrc.
        
        Returns:
            True if successful
        """
        bashrc_path = Path.home() / ".bashrc"
        marker = "# COMPASS Modulefiles Configuration"
        
        # Check if already added
        if bashrc_path.exists():
            with open(bashrc_path, 'r') as f:
                if marker in f.read():
                    self.logger.info("COMPASS modulefile configuration already in .bashrc")
                    return True
        
        # Add configuration
        with open(bashrc_path, 'a') as f:
            f.write(f"\n{marker}\n")
            f.write(f"export MODULEPATH={self.modulefiles_root}:$MODULEPATH\n")
        
        self.logger.info("✓ COMPASS modulefile configuration added to .bashrc")
        return True
    
    def is_compass_module_loaded(self) -> bool:
        """
        Check if COMPASS module is currently loaded.
        
        Returns:
            True if compass/local module is loaded
        """
        try:
            result = subprocess.run(
                ["bash", "-c", "module list 2>&1"],
                capture_output=True,
                text=True,
                timeout=5
            )
            output = result.stdout + result.stderr
            return "compass/local" in output.lower()
        except Exception:
            return False
    
    def check_and_setup(self) -> bool:
        """
        Check if Lmod is installed and setup modulefiles.
        
        Returns:
            True if Lmod is ready and modulefiles are setup
        """
        self.console.print("[bold blue]Checking Lmod installation...[/bold blue]\n")
        
        if not self.is_lmod_installed():
            self.console.print("[bold red]✗ Lmod is not installed![/bold red]\n")
            self.console.print("Lmod must be installed via your system package manager:")
            self.console.print("  [cyan]# RHEL/Rocky/CentOS:[/cyan]")
            self.console.print("  [cyan]sudo dnf install Lmod[/cyan]\n")
            self.console.print("  [cyan]# Debian/Ubuntu:[/cyan]")
            self.console.print("  [cyan]sudo apt install lmod[/cyan]\n")
            return False
        
        self.console.print("[green]✓ Lmod is installed[/green]")
        
        # Check modulefiles
        if not self.modulefiles_root.exists():
            self.console.print(f"[yellow]⚠ Modulefiles directory not found: {self.modulefiles_root}[/yellow]")
            return False
        
        self.console.print(f"[green]✓ Modulefiles found at {self.modulefiles_root}[/green]")
        
        # Show available modules
        modules = self.get_available_modules()
        if modules:
            self.console.print(f"\n[bold]Available modules:[/bold]")
            for module in modules:
                self.console.print(f"  • {module}")
        
        self.console.print()
        self.console.print("[bold green]✓ Lmod is ready![/bold green]")
        self.console.print("\n[yellow]To use COMPASS modules:[/yellow]")
        
        # Check if already in MODULEPATH
        modulepath = os.environ.get("MODULEPATH", "")
        if str(self.modulefiles_root) not in modulepath:
            self.console.print("  [cyan]1. Add to your shell:[/cyan]")
            self.console.print(f"     [cyan]export MODULEPATH={self.modulefiles_root}:$MODULEPATH[/cyan]")
            self.console.print("  [cyan]2. Or add to ~/.bashrc (automatic):[/cyan]")
            self.console.print(f"     [cyan]compass lmod setup[/cyan]")
        
        self.console.print("  [cyan]3. Load the module:[/cyan]")
        self.console.print("     [cyan]module load compass/local[/cyan]\n")
        
        return True
    
    def setup(self) -> bool:
        """
        Setup Lmod for COMPASS (add MODULEPATH to bashrc).
        
        Returns:
            True if successful
        """
        self.console.print("[bold blue]🔧 Setting up COMPASS Modulefiles[/bold blue]\n")
        
        if not self.is_lmod_installed():
            self.console.print("[bold red]✗ Lmod is not installed[/bold red]")
            self.console.print("\nPlease install Lmod first:")
            self.console.print("  Ubuntu/Debian: [cyan]sudo apt-get install lmod[/cyan]")
            self.console.print("  CentOS/RHEL: [cyan]sudo yum install Lmod[/cyan]")
            self.console.print("  From source: [cyan]https://lmod.readthedocs.io[/cyan]")
            return False
        
        self.console.print("[green]✓ Lmod is installed[/green]")
        
        if not self.modulefiles_root.exists():
            self.console.print(f"[bold red]✗ Modulefiles directory not found: {self.modulefiles_root}[/bold red]")
            return False
        
        self.console.print(f"[green]✓ Modulefiles directory exists: {self.modulefiles_root}[/green]")
        
        # Check if modulefiles_root is already configured in bashrc
        bashrc_path = Path.home() / ".bashrc"
        modulepath_export = f"export MODULEPATH={self.modulefiles_root}"
        already_in_bashrc = False
        
        if bashrc_path.exists():
            with open(bashrc_path, 'r') as f:
                bashrc_content = f.read()
                if modulepath_export in bashrc_content or str(self.modulefiles_root) in bashrc_content:
                    already_in_bashrc = True
        
        if already_in_bashrc:
            self.console.print("[green]✓ COMPASS modulefiles already configured in ~/.bashrc[/green]")
            
            # Check if also in current MODULEPATH
            modulepath = os.environ.get("MODULEPATH", "")
            if str(self.modulefiles_root) in modulepath:
                self.console.print("[green]✓ COMPASS modulefiles in current MODULEPATH[/green]")
                self.console.print("\n[bold green]✓ Setup complete![/bold green]")
                self.console.print("\nYou can now load the module:")
                self.console.print("  [cyan]module load compass/local[/cyan]")
            else:
                self.console.print("[yellow]⚠ Not active in current session[/yellow]")
                self.console.print("\n[bold green]✓ Setup complete![/bold green]")
                self.console.print("\n[yellow]Next steps:[/yellow]")
                self.console.print("  1. [cyan]source ~/.bashrc[/cyan]  (or restart your shell)")
                self.console.print("  2. [cyan]module load compass/local[/cyan]")
            return True
        
        # Add to bashrc
        if self.add_to_bashrc():
            self.console.print("[green]✓ Added COMPASS modulefiles to ~/.bashrc[/green]")
            self.console.print("\n[bold green]✓ Setup complete![/bold green]")
            self.console.print("\n[yellow]Next steps:[/yellow]")
            self.console.print("  1. [cyan]source ~/.bashrc[/cyan]  (or restart your shell)")
            self.console.print("  2. [cyan]module load compass/local[/cyan]")
            return True
        
        return False
    
    def show_status(self):
        """Show Lmod status and available modules."""
        self.console.print("\n[bold cyan]Lmod Status[/bold cyan]\n")
        
        # Check if Lmod is installed
        if self.is_lmod_installed():
            self.console.print("[green]✓ Lmod is installed[/green]")
            
            # Show LMOD_CMD
            lmod_cmd = os.environ.get("LMOD_CMD")
            if lmod_cmd:
                self.console.print(f"  LMOD_CMD: {lmod_cmd}")
        else:
            self.console.print("[red]✗ Lmod is not installed[/red]")
        
        # Check modulefiles
        self.console.print()
        if self.modulefiles_root.exists():
            self.console.print(f"[green]✓ Modulefiles directory: {self.modulefiles_root}[/green]")
            
            modules = self.get_available_modules()
            if modules:
                self.console.print(f"\n[bold]Available modules ({len(modules)}):[/bold]")
                for module in modules:
                    self.console.print(f"  • {module}")
            else:
                self.console.print("[yellow]⚠ No module files found[/yellow]")
        else:
            self.console.print(f"[red]✗ Modulefiles directory not found: {self.modulefiles_root}[/red]")
        
        # Check MODULEPATH
        self.console.print()
        modulepath = os.environ.get("MODULEPATH", "")
        if str(self.modulefiles_root) in modulepath:
            self.console.print("[green]✓ COMPASS modulefiles in MODULEPATH[/green]")
        else:
            self.console.print("[yellow]⚠ COMPASS modulefiles NOT in MODULEPATH[/yellow]")
            self.console.print(f"  Run: [cyan]compass lmod setup[/cyan]")


def main():
    """Command-line entry point for Lmod management."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Manage COMPASS Lmod configuration")
    parser.add_argument("action", choices=["check", "setup", "status"], help="Action to perform")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    manager = LmodManager(verbose=args.verbose)
    
    if args.action == "check":
        success = manager.check_and_setup()
    elif args.action == "setup":
        success = manager.setup()
    elif args.action == "status":
        manager.show_status()
        success = True
    
    import sys
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
