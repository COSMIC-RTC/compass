"""
Simulation runner for COMPASS.

Provides interface to run COMPASS/Shesha simulations.
"""

import os
import subprocess
import json
from pathlib import Path
from typing import Optional, List
from rich.console import Console
from rich.table import Table

from compass.core.logger import LoggerMixin


class SimulationRunner(LoggerMixin):
    """
    Manages running COMPASS simulations.
    
    Wraps Shesha simulation scripts and provides enhanced monitoring.
    """
    
    def __init__(self, verbose: bool = True):
        """
        Initialize the simulation runner.
        
        Args:
            verbose: Enable verbose output
        """
        self.verbose = verbose
        self.console = Console()
        self.compass_root = Path(os.getenv("COMPASS_ROOT", Path.home() / "compass"))
        self.shesha_root = Path(os.getenv("SHESHA_ROOT", self.compass_root / "shesha"))
        self.config_file = Path.home() / ".compass" / "sim_config.json"
        self.default_script = "closed_loop.py"
        self.default_gui_script = "compass_gui.py"
    
    def get_script_path(self, script_name: Optional[str] = None, script_type: str = "default") -> Path:
        """
        Get the full path to a script.
        
        Args:
            script_name: Name or path of the script (optional, uses configured default if not provided)
            script_type: Type of script - "default" or "gui"
        
        Returns:
            Path to the script
        """
        if script_name is None:
            # Load from config
            config = self.load_config()
            if script_type == "gui":
                script_name = config.get("default_gui_script", self.default_gui_script)
            else:
                script_name = config.get("default_script", self.default_script)
        
        # If it's already a full path, use it
        script_path = Path(script_name)
        if script_path.is_absolute():
            return script_path
        
        # Otherwise, look in appropriate directory
        # Both GUI and regular scripts are in shesha/scripts
        script_path = self.shesha_root / "shesha" / "scripts" / script_name
        
        if not script_path.exists():
            # For backwards compatibility, also check widgets directory for GUI scripts
            if script_type == "gui":
                widget_path = self.shesha_root / "shesha" / "widgets" / script_name
                if widget_path.exists():
                    return widget_path
            
            # Try without .py extension
            if not script_name.endswith(".py"):
                alt_path = self.shesha_root / "shesha" / "scripts" / f"{script_name}.py"
                if alt_path.exists():
                    return alt_path
                
                # Also check widgets for GUI scripts (backwards compatibility)
                if script_type == "gui":
                    widget_alt_path = self.shesha_root / "shesha" / "widgets" / f"{script_name}.py"
                    if widget_alt_path.exists():
                        return widget_alt_path
        
        return script_path
    
    def load_config(self) -> dict:
        """Load simulation configuration from file."""
        if self.config_file.exists():
            try:
                with open(self.config_file, 'r') as f:
                    return json.load(f)
            except Exception:
                pass
        return {}
    
    def save_config(self, config: dict):
        """Save simulation configuration to file."""
        self.config_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self.config_file, 'w') as f:
            json.dump(config, f, indent=2)
    
    def set_default_script(self, script_name: str, script_type: str = "default") -> bool:
        """
        Set the default script for simulations.
        
        Args:
            script_name: Name or absolute path of the script
            script_type: Type of script - "default" or "gui"
        
        Returns:
            True if successful
        """
        script_path = self.get_script_path(script_name, script_type=script_type)
        
        if not script_path.exists():
            self.console.print(f"[red]✗ Script not found: {script_path}[/red]")
            return False
        
        # Store the absolute path in config
        config = self.load_config()
        if script_type == "gui":
            config["default_gui_script"] = str(script_path.absolute())
            script_label = "GUI script"
        else:
            config["default_script"] = str(script_path.absolute())
            script_label = "script"
        self.save_config(config)
        
        self.console.print(f"[green]✓ Default {script_label} set to: {script_path.name}[/green]")
        self.console.print(f"  Full path: {script_path.absolute()}")
        
        return True
    
    def show_script_config(self):
        """Show current script configuration."""
        script_path = self.get_script_path()
        
        self.console.print("\n[bold cyan]Current Simulation Script Configuration[/bold cyan]\n")
        
        table = Table(title="Script Settings")
        table.add_column("Setting", style="cyan")
        table.add_column("Value", style="green")
        
        # Show just the name if it's in shesha/scripts, otherwise show full path
        if str(self.shesha_root / "shesha" / "scripts") in str(script_path):
            display_name = script_path.name
        else:
            display_name = str(script_path)
        
        table.add_row("Default Script", display_name)
        table.add_row("Full Path", str(script_path.absolute()))
        table.add_row("Exists", "✓" if script_path.exists() else "✗")
        
        self.console.print(table)
        
        # Show available scripts in shesha/scripts
        scripts_dir = self.shesha_root / "shesha" / "scripts"
        if scripts_dir.exists():
            scripts = sorted([s.name for s in scripts_dir.glob("*.py") if s.name != "__init__.py"])
            if scripts:
                self.console.print(f"\n[dim]Available scripts in {scripts_dir.relative_to(self.compass_root)}:[/dim]")
                for script in scripts:
                    full_script_path = scripts_dir / script
                    marker = " [default]" if str(full_script_path.absolute()) == str(script_path.absolute()) else ""
                    self.console.print(f"  • {script}{marker}")
    
    def list_available_scripts(self):
        """List all available scripts in shesha/scripts."""
        scripts_dir = self.shesha_root / "shesha" / "scripts"
        
        if not scripts_dir.exists():
            self.console.print(f"[red]✗ Scripts directory not found: {scripts_dir}[/red]")
            return
        
        scripts = sorted([s for s in scripts_dir.glob("*.py") if s.name != "__init__.py"])
        
        if not scripts:
            self.console.print(f"[yellow]No scripts found in {scripts_dir}[/yellow]")
            return
        
        default_script_path = str(self.get_script_path().absolute())
        
        self.console.print("\n[bold cyan]Available Simulation Scripts[/bold cyan]\n")
        
        table = Table(title=f"Scripts in {scripts_dir.relative_to(self.compass_root)}")
        table.add_column("Script", style="cyan")
        table.add_column("Status", style="green")
        
        for script in scripts:
            status = "✓ [default]" if str(script.absolute()) == default_script_path else ""
            table.add_row(script.name, status)
        
        self.console.print(table)
        self.console.print("\n[dim]Set default with: compass sim script set <script_name_or_path>[/dim]")
        self.console.print("[dim]You can also use absolute paths outside this directory[/dim]")

    
    def _run_simulation(
        self,
        param_file: Optional[str],
        script_type: str = "default",
        iterations: Optional[int] = None,
        devices: Optional[str] = None,
        script_args: Optional[List[str]] = None,
    ) -> bool:
        """
        Run a COMPASS simulation with configured script.
        
        Args:
            param_file: Path to parameter file (optional for GUI mode)
            script_type: Type of script ("default" for CLI or "gui" for GUI)
            iterations: Number of iterations (optional)
            devices: GPU devices to use (comma-separated, optional)
            script_args: Additional arguments to pass to the script (optional)
        
        Returns:
            True if simulation completed successfully
        """
        # Handle parameter file path (optional for GUI mode)
        param_path = None
        if param_file:
            param_path = Path(param_file)
            if not param_path.is_absolute():
                # Try relative to current directory first
                if not param_path.exists():
                    # Try relative to COMPASS root
                    param_path = self.compass_root / param_file
            
            if not param_path.exists():
                self.console.print(f"[red]✗ Parameter file not found: {param_file}[/red]")
                return False
        
        # Get the configured script based on type
        script_path = self.get_script_path(script_type=script_type)
        
        if not script_path.exists():
            script_cmd = "script" if script_type == "default" else "gui-script"
            self.console.print(f"[red]✗ Script not found: {script_path}[/red]")
            self.console.print(f"[yellow]Run 'compass sim {script_cmd} set <script_name_or_path>' to set a valid script[/yellow]")
            return False
        
        # Configure display based on script type
        is_gui = script_type == "gui"
        title = "🖥️  Running COMPASS GUI Simulation" if is_gui else "🔭 Running COMPASS Simulation (Interactive)"
        success_msg = "GUI session ended" if is_gui else "Simulation session ended"
        
        self.console.print(f"\n[bold cyan]{title}[/bold cyan]")
        if param_path:
            self.console.print(f"[cyan]Parameter file: {param_path.absolute()}[/cyan]")
        
        # Show script name if in shesha directory, otherwise show full path
        shesha_scripts_dir = str(self.shesha_root / "shesha" / "scripts")
        shesha_widgets_dir = str(self.shesha_root / "shesha" / "widgets")
        
        if shesha_scripts_dir in str(script_path) or shesha_widgets_dir in str(script_path):
            script_label = "GUI script" if is_gui else "script"
            self.console.print(f"[cyan]Using {script_label}: {script_path.name}[/cyan]\n")
        else:
            self.console.print(f"[cyan]Using script: {script_path.absolute()}[/cyan]\n")
        
        # Build command based on script type
        if is_gui:
            # GUI: run the script directly with python
            # compass_gui.py handles its own GUI setup and event loop
            cmd = ["python", str(script_path.absolute())]
            if param_path:
                cmd.append(str(param_path.absolute()))
        else:
            # CLI: use ipython for interactive session
            if not param_path:
                self.console.print("[red]✗ Parameter file required for CLI mode[/red]")
                return False
            cmd = ["ipython", "-i", str(script_path.absolute()), "--", str(param_path.absolute())]
        
        # Add any additional script arguments
        if script_args:
            cmd.extend(script_args)
            self.console.print(f"[cyan]Script args: {' '.join(script_args)}[/cyan]")
        
        # Set environment
        env = os.environ.copy()
        
        if devices:
            env["CUDA_VISIBLE_DEVICES"] = devices
            self.console.print(f"[cyan]GPU devices: {devices}[/cyan]")
        
        if iterations:
            env["SHESHA_ITERATIONS"] = str(iterations)
            self.console.print(f"[cyan]Iterations: {iterations}[/cyan]")
        
        self.console.print()
        
        # Run simulation
        try:
            result = subprocess.run(
                cmd,
                env=env,
                cwd=self.compass_root,
            )
            
            if result.returncode == 0:
                self.console.print(f"\n[bold green]✓ {success_msg}[/bold green]")
                return True
            else:
                # Non-zero exit is common for interactive sessions (Ctrl+D)
                if not is_gui:
                    return True
                else:
                    self.console.print(f"\n[red]✗ GUI failed with code {result.returncode}[/red]")
                    return False
                
        except KeyboardInterrupt:
            self.console.print("\n[yellow]⚠ Simulation interrupted by user[/yellow]")
            return False
        except Exception as e:
            self.console.print(f"\n[red]✗ Error running simulation: {e}[/red]")
            return False
    
    def run(
        self,
        param_file: str,
        iterations: Optional[int] = None,
        devices: Optional[str] = None,
        script_args: Optional[List[str]] = None,
    ) -> bool:
        """
        Run a COMPASS simulation with configured script in interactive IPython session.
        
        Args:
            param_file: Path to parameter file (relative or absolute)
            iterations: Number of iterations (optional)
            devices: GPU devices to use (comma-separated, optional)
            script_args: Additional arguments to pass to the script (optional)
        
        Returns:
            True if simulation completed successfully
        """
        return self._run_simulation(param_file, "default", iterations, devices, script_args)
    
    def gui(
        self,
        param_file: Optional[str] = None,
        param_file2: Optional[str] = None,
        frequency_ratio: int = 1,
        iterations: Optional[int] = None,
        devices: Optional[str] = None,
        script_args: Optional[List[str]] = None,
    ) -> bool:
        """
        Run a COMPASS simulation with GUI using configured GUI script.
        
        Args:
            param_file: Path to parameter file (optional)
            param_file2: Optional second parameter file for two-stages mode
            frequency_ratio: Frequency ratio for two-stages mode (default: 1)
            iterations: Number of iterations (optional)
            devices: GPU devices to use (comma-separated, optional)
            script_args: Additional arguments to pass to the script (optional)
        
        Returns:
            True if simulation completed successfully
        """
        # Build script args to include two-stages parameters if provided
        combined_args = script_args or []
        if param_file2:
            # For two-stages mode, pass both files and frequency ratio
            combined_args = [param_file2, str(frequency_ratio)] + combined_args
        
        return self._run_simulation(param_file, "gui", iterations, devices, combined_args)
    
    def gui_server(
        self,
        param_file: str,
        param_file2: Optional[str] = None,
        frequency_ratio: int = 1,
        command_port: int = 5555,
        telemetry_port: int = 5556,
        bind_address: str = "*",
        iterations: Optional[int] = None,
        devices: Optional[str] = None,
    ) -> bool:
        """
        Run a COMPASS simulation with GUI in server mode for remote control.
        
        Args:
            param_file: Path to parameter file (relative or absolute)
            param_file2: Optional second parameter file for two-stages mode
            frequency_ratio: Frequency ratio for two-stages mode (default: 1)
            command_port: Port for command socket (default: 5555)
            telemetry_port: Port for telemetry socket (default: 5556)
            bind_address: Address to bind to (default: "*" for all interfaces)
            iterations: Number of iterations (optional)
            devices: GPU devices to use (comma-separated, optional)
        
        Returns:
            True if server started successfully
        """
        # Handle both relative and absolute paths for parameter file
        param_path = Path(param_file)
        if not param_path.is_absolute():
            # Try relative to current directory first
            if not param_path.exists():
                # Try relative to COMPASS root
                param_path = self.compass_root / param_file
        
        if not param_path.exists():
            self.console.print(f"[red]✗ Parameter file not found: {param_file}[/red]")
            return False
        
        # Handle second parameter file for two-stages mode
        param_path2 = None
        is_two_stages = param_file2 is not None
        
        if is_two_stages:
            param_path2 = Path(param_file2)
            if not param_path2.is_absolute():
                if not param_path2.exists():
                    param_path2 = self.compass_root / param_file2
            
            if not param_path2.exists():
                self.console.print(f"[red]✗ Second parameter file not found: {param_file2}[/red]")
                return False
        
        # Get the remote server script
        server_script = self.shesha_root / "shesha" / "scripts" / "run_remote_server.py"
        
        if not server_script.exists():
            self.console.print(f"[red]✗ Remote server script not found: {server_script}[/red]")
            return False
        
        if is_two_stages:
            self.console.print("\n[bold cyan]🖥️  Starting COMPASS Remote GUI Server (Two-Stages Mode)[/bold cyan]")
            self.console.print(f"[cyan]First stage: {param_path.absolute()}[/cyan]")
            self.console.print(f"[cyan]Second stage: {param_path2.absolute()}[/cyan]")
            self.console.print(f"[cyan]Frequency ratio: {frequency_ratio}[/cyan]")
        else:
            self.console.print("\n[bold cyan]🖥️  Starting COMPASS Remote GUI Server[/bold cyan]")
            self.console.print(f"[cyan]Parameter file: {param_path.absolute()}[/cyan]")
        
        self.console.print(f"[cyan]Command port: {command_port}[/cyan]")
        self.console.print(f"[cyan]Telemetry port: {telemetry_port}[/cyan]")
        self.console.print(f"[cyan]Bind address: {bind_address}[/cyan]\n")
        
        # Build command
        cmd = [
            "python",
            str(server_script.absolute()),
            str(param_path.absolute()),
        ]
        
        # Add second param file and frequency ratio if two-stages
        if is_two_stages:
            cmd.append(str(param_path2.absolute()))
            cmd.append(str(frequency_ratio))
        
        cmd.extend([
            "--host", bind_address,
            "--cmd-port", str(command_port),
            "--tel-port", str(telemetry_port),
        ])
        
        # Set environment
        env = os.environ.copy()
        
        if devices:
            env["CUDA_VISIBLE_DEVICES"] = devices
            self.console.print(f"[cyan]GPU devices: {devices}[/cyan]")
        
        if iterations:
            env["SHESHA_ITERATIONS"] = str(iterations)
            self.console.print(f"[cyan]Iterations: {iterations}[/cyan]")
        
        self.console.print()
        
        # Run server
        try:
            result = subprocess.run(
                cmd,
                env=env,
                cwd=self.compass_root,
            )
            
            if result.returncode == 0:
                self.console.print("\n[bold green]✓ Server shut down successfully[/bold green]")
                return True
            else:
                self.console.print(f"\n[red]✗ Server failed with code {result.returncode}[/red]")
                return False
                
        except KeyboardInterrupt:
            self.console.print("\n[yellow]⚠ Server interrupted by user[/yellow]")
            return True
        except Exception as e:
            self.console.print(f"\n[red]✗ Error running server: {e}[/red]")
            return False
    
    def gui_remote(
        self,
        host: str,
        command_port: int = 5555,
        telemetry_port: int = 5556,
    ) -> bool:
        """
        Launch a remote GUI client that connects to a COMPASS server.
        
        Args:
            host: Hostname or IP address of the remote server
            command_port: Command port of the remote server (default: 5555)
            telemetry_port: Telemetry port of the remote server (default: 5556)
        
        Returns:
            True if GUI client started successfully
        """
        # Get the remote client script
        client_script = self.shesha_root / "shesha" / "scripts" / "run_remote_client.py"
        
        if not client_script.exists():
            self.console.print(f"[red]✗ Remote client script not found: {client_script}[/red]")
            return False
        
        self.console.print("\n[bold cyan]🖥️  Starting COMPASS Remote GUI Client[/bold cyan]")
        self.console.print(f"[cyan]Server: {host}:{command_port}[/cyan]")
        self.console.print(f"[cyan]Telemetry: {host}:{telemetry_port}[/cyan]\n")
        
        # Build command
        cmd = [
            "python",
            str(client_script.absolute()),
            "--server", host,
            "--cmd-port", str(command_port),
            "--tel-port", str(telemetry_port),
            "--auto-connect",
        ]
        
        # Run GUI client
        try:
            result = subprocess.run(
                cmd,
                cwd=self.compass_root,
            )
            
            if result.returncode == 0:
                self.console.print("\n[bold green]✓ GUI client ended[/bold green]")
                return True
            else:
                # Non-zero exit is common for GUI apps
                return True
                
        except KeyboardInterrupt:
            self.console.print("\n[yellow]⚠ GUI interrupted by user[/yellow]")
            return False
        except Exception as e:
            self.console.print(f"\n[red]✗ Error running GUI client: {e}[/red]")
            return False
    
    def list_examples(self, directory: Optional[str] = None):
        """
        List available example simulations.
        
        Args:
            directory: Specific directory to list parfiles from (optional)
        """
        par_dir = self.shesha_root / "data" / "par"
        
        if not par_dir.exists():
            self.console.print(f"[red]✗ Parameter directory not found: {par_dir}[/red]")
            self.console.print("[yellow]Make sure SHESHA_ROOT is set correctly[/yellow]")
            return
        
        # If no directory specified, list available directories
        if directory is None:
            self.console.print("[bold cyan]Available Parameter Directories:[/bold cyan]\n")
            
            subdirs = [d for d in par_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
            
            if subdirs:
                table = Table(title="Shesha Parameter Directories")
                table.add_column("Directory", style="cyan")
                table.add_column("Path", style="green")
                table.add_column("Files", style="yellow")
                
                for subdir in sorted(subdirs):
                    # Count .py files in directory
                    py_files = list(subdir.glob("*.py"))
                    table.add_row(
                        subdir.name,
                        str(subdir.relative_to(self.shesha_root)),
                        str(len(py_files))
                    )
                
                self.console.print(table)
                self.console.print("\n[dim]Use 'compass sim list <DIRECTORY>' to see parfiles in a specific directory[/dim]")
                self.console.print("[dim]Example: compass sim list MICADO[/dim]")
            else:
                self.console.print(f"[yellow]No subdirectories found in {par_dir}[/yellow]")
        
        # If directory specified, list parfiles in that directory
        else:
            target_dir = par_dir / directory
            
            if not target_dir.exists():
                self.console.print(f"[red]✗ Directory not found: {directory}[/red]")
                self.console.print("[yellow]Use 'compass sim list' to see available directories[/yellow]")
                return
            
            self.console.print(f"[bold cyan]Parameter Files in {directory}:[/bold cyan]\n")
            
            parfiles = sorted(target_dir.glob("*.py"))
            
            if parfiles:
                table = Table(title=f"{directory} Parameter Files")
                table.add_column("File", style="cyan")
                table.add_column("Path", style="green")
                
                for parfile in parfiles:
                    table.add_row(
                        parfile.name,
                        str(parfile.relative_to(self.shesha_root))
                    )
                
                self.console.print(table)
                self.console.print(f"\n[dim]Run with: compass sim run {parfiles[0].relative_to(self.compass_root)}[/dim]")
            else:
                self.console.print(f"[yellow]No .py files found in {target_dir}[/yellow]")
