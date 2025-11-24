"""
Build and compilation management module for COMPASS.

Handles building libcarma, libsutra, and Python wrappers.
"""

import os
import sys
import logging
import subprocess
from pathlib import Path
from typing import List, Optional
from dataclasses import dataclass, field
from enum import Enum
from rich.console import Console
from rich.table import Table

from compass.core.environment import EnvironmentChecker
from compass.core.logger import LoggerMixin


class BuildType(Enum):
    """Build type enumeration."""
    RELEASE = "Release"
    DEBUG = "Debug"
    RELWITHDEBINFO = "RelWithDebInfo"


@dataclass
class Component:
    """Represents a buildable COMPASS component."""
    
    name: str
    root_path: Path
    build_dir: Optional[Path] = None
    install_path: Optional[Path] = None
    required: bool = True
    dependencies: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        if self.build_dir is None:
            self.build_dir = self.root_path / "build"


class Builder(LoggerMixin):
    """
    Handles building and compilation of COMPASS components.
    
    Components:
    - libcarma: CUDA-based AO Real-time Modules
    - libsutra: Simulation Utilities for Tomographic Reconstruction
    - python_module: Python wrappers (shesha)
    """
    
    def __init__(self, verbose: bool = False):
        """
        Initialize the Builder.
        
        Args:
            verbose: Enable verbose output
        """
        self.verbose = verbose
        self.console = Console()
        self.env_checker = EnvironmentChecker(verbose=verbose)
        
        # Get compass root
        self.compass_root = Path(os.getenv("COMPASS_ROOT", Path.home() / "compass"))
        self.build_dir = self.compass_root / "build"
        self.install_dir = Path(os.getenv(
            "COMPASS_INSTALL_ROOT",
            self.compass_root / "local"
        ))
        
        # Setup logger via mixin
        if verbose:
            self.logger.setLevel(logging.DEBUG)
    
    def get_components(self) -> List[Component]:
        """
        Get list of components to build in correct dependency order.
        
        Returns:
            List of Component objects
        """
        components = [
            Component(
                name="libcarma",
                root_path=self.compass_root / "libcarma",
                install_path=self.install_dir,
                required=True,
                dependencies=[],
            ),
            Component(
                name="libsutra",
                root_path=self.compass_root / "libsutra",
                install_path=self.install_dir,
                required=True,
                dependencies=["libcarma"],
            ),
            Component(
                name="python_module",
                root_path=self.compass_root / "python_module",
                install_path=self.install_dir,
                required=True,
                dependencies=["libcarma", "libsutra"],
            ),
        ]
        
        return components
    
    def clean_build(self) -> bool:
        """
        Clean main build directory.
        
        Returns:
            True if successful
        """
        if self.build_dir.exists():
            import shutil
            try:
                self.console.print(f"[yellow]Cleaning build directory: {self.build_dir}[/yellow]")
                shutil.rmtree(self.build_dir)
                self.logger.info("Cleaned build directory")
                return True
            except Exception as e:
                self.logger.error(f"Failed to clean build directory: {e}")
                return False
        return True
    
    def clean_install(self) -> bool:
        """
        Clean install directory.
        
        Returns:
            True if successful
        """
        if self.install_dir.exists():
            import shutil
            try:
                self.console.print(f"[yellow]Cleaning install directory: {self.install_dir}[/yellow]")
                shutil.rmtree(self.install_dir)
                self.logger.info("Cleaned install directory")
                return True
            except Exception as e:
                self.logger.error(f"Failed to clean install directory: {e}")
                return False
        return True
    
    def _run_cmake_configure(self) -> bool:
        """Run CMake configuration."""
        self.console.print("[cyan]Configuring with CMake...[/cyan]")
        
        # Ensure build directory exists
        self.build_dir.mkdir(parents=True, exist_ok=True)
        
        cmake_args = [
            "cmake",
            "-S", str(self.compass_root),
            "-B", str(self.build_dir),
            f"-DCMAKE_INSTALL_PREFIX={self.install_dir}",
            "-DCMAKE_BUILD_TYPE=Release",
        ]
        
        # Add CUDA support if available
        cuda_root = os.getenv("CUDA_ROOT")
        if cuda_root:
            cmake_args.append(f"-DCUDA_TOOLKIT_ROOT_DIR={cuda_root}")

        # Add pybind11 support
        try:
            pybind11_dir = subprocess.check_output(
                [sys.executable, "-m", "pybind11", "--cmakedir"],
                text=True
            ).strip()
            cmake_args.append(f"-Dpybind11_DIR={pybind11_dir}")
        except subprocess.CalledProcessError:
            self.logger.warning("Could not determine pybind11 CMake directory")
        
        return self._run_command(cmake_args, cwd=self.compass_root)
    
    def _run_cmake_build(self, targets: Optional[List[str]] = None) -> bool:
        """
        Run CMake build.
        
        Args:
            targets: List of specific targets to build (None = all)
        
        Returns:
            True if build succeeded
        """
        if targets:
            self.console.print(f"[cyan]Building targets: {', '.join(targets)}...[/cyan]")
        else:
            self.console.print("[cyan]Building all targets...[/cyan]")
        
        # Detect number of cores
        import multiprocessing
        n_cores = multiprocessing.cpu_count()
        
        cmake_args = [
            "cmake",
            "--build", str(self.build_dir),
            "--parallel", str(n_cores),
        ]
        
        # Add specific targets if requested
        if targets:
            for target in targets:
                cmake_args.extend(["--target", target])
        
        return self._run_command(cmake_args, cwd=self.build_dir)
    
    def _run_cmake_install(self) -> bool:
        """Run CMake install."""
        self.console.print("[cyan]Installing...[/cyan]")
        
        cmake_args = [
            "cmake",
            "--install", str(self.build_dir),
        ]
        
        return self._run_command(cmake_args, cwd=self.build_dir)
    
    def _run_command(self, cmd: List[str], cwd: Path) -> bool:
        """
        Run a shell command.
        
        Args:
            cmd: Command and arguments
            cwd: Working directory
        
        Returns:
            True if command succeeded
        """
        if self.verbose:
            self.console.print(f"[dim]Running: {' '.join(cmd)}[/dim]")
        
        try:
            if self.verbose:
                # Stream output in real-time when verbose
                process = subprocess.Popen(
                    cmd,
                    cwd=cwd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                
                for line in process.stdout:
                    print(line, end='')
                
                process.wait()
                
                if process.returncode != 0:
                    raise subprocess.CalledProcessError(process.returncode, cmd)
            else:
                # Capture output silently when not verbose
                subprocess.run(
                    cmd,
                    cwd=cwd,
                    capture_output=True,
                    text=True,
                    check=True,
                )
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.console.print(f"[red]✗ Command failed with code {e.returncode}[/red]")
            if hasattr(e, 'stderr') and e.stderr:
                self.console.print(f"[red]{e.stderr}[/red]")
            return False
        except Exception as e:
            self.console.print(f"[red]✗ Error running command: {e}[/red]")
            return False
    
    def build_all(
        self,
        clean: bool = False,
        components: Optional[List[str]] = None,
    ) -> bool:
        """
        Build all or specified COMPASS components using CMake.
        
        Args:
            clean: Clean build directories before building
            components: List of specific components to build (None = all)
                       Available: libcarma, libsutra, python_module
        
        Returns:
            True if build succeeded
        """
        if components:
            self.console.print(f"\n[bold cyan]🔨 Building COMPASS Components: {', '.join(components)}[/bold cyan]\n")
        else:
            self.console.print("\n[bold cyan]🔨 Building COMPASS Components[/bold cyan]\n")
        
        # Validate component names
        valid_components = {"libcarma", "libsutra", "python_module"}
        if components:
            invalid = set(components) - valid_components
            if invalid:
                self.console.print(f"[red]✗ Invalid components: {', '.join(invalid)}[/red]")
                self.console.print(f"[yellow]Valid components: {', '.join(sorted(valid_components))}[/yellow]")
                return False
        
        # Clean if requested     
        if clean:
            self.clean_build()
            self.clean_install()
        
        # Configure (always needed, even for partial builds)
        if not self._run_cmake_configure():
            return False
        
        # Build with specific targets or all
        if not self._run_cmake_build(targets=components):
            return False
        
        # Install
        if not self._run_cmake_install():
            return False
        
        self.console.print("\n[bold green]✓ Build completed successfully![/bold green]")
        self.console.print(f"[cyan]Installation: {self.install_dir}[/cyan]")
        return True
    
    def show_status(self):
        """Show build status information."""
        console = Console()
        
        console.print("\n[bold cyan]Build Configuration Status[/bold cyan]\n")
        
        table = Table(title="COMPASS Build Paths")
        table.add_column("Component", style="cyan")
        table.add_column("Path", style="green")
        table.add_column("Status", style="white")
        
        # Check paths
        paths_to_check = [
            ("COMPASS Root", self.compass_root),
            ("Build Directory", self.build_dir),
            ("Install Directory", self.install_dir),
        ]
        
        for name, path in paths_to_check:
            status = "✓" if path.exists() else "✗"
            table.add_row(name, str(path), status)
        
        console.print(table)
        
        # Check for built libraries
        console.print("\n[bold]Built Libraries:[/bold]")
        lib_dir = self.install_dir / "lib"
        if lib_dir.exists():
            libs = list(lib_dir.glob("*.so")) + list(lib_dir.glob("*.a"))
            if libs:
                for lib in sorted(libs)[:10]:  # Show first 10
                    console.print(f"  [green]✓ {lib.name}[/green]")
                if len(libs) > 10:
                    console.print(f"  [dim]... and {len(libs) - 10} more[/dim]")
            else:
                console.print("  [yellow]⚠ No libraries found[/yellow]")
        else:
            console.print("  [red]✗ Library directory not found[/red]")


def main():
    """Command-line entry point for building."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Build COMPASS components")
    parser.add_argument("--clean", action="store_true", help="Clean before building")
    parser.add_argument("-c", "--components", nargs="*", help="Specific components to build")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    builder = Builder(verbose=args.verbose)
    success = builder.build_all(
        clean=args.clean,
        components=args.components,
    )
    
    import sys
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
