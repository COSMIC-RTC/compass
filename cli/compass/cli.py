"""
Unified CLI interface for COMPASS.

Provides a single entry point for all COMPASS operations.
"""

import sys
import os
import argparse
from pathlib import Path
from rich.console import Console

from compass import __version__
from compass.core.environment import EnvironmentChecker
from compass.core.config import Config, Environment
from compass.deployment.build import Builder
from compass.simulation.runner import SimulationRunner
from compass.system.lmod import LmodManager


console = Console()


# Command handlers

def cmd_build(args):
    """Build COMPASS components."""
    verbose = not args.silent
    
    builder = Builder(verbose=verbose)
    success = builder.build_all(
        clean=args.clean,
        components=list(args.components) if args.components else None,
    )
    
    sys.exit(0 if success else 1)


def cmd_config_export(args):
    """Export environment configuration to shell script."""
    verbose = not args.silent
    
    env = Environment()
    output_path = Path(args.output)
    env.export_to_shell_script(output_path)
    
    if not args.silent:
        console.print(f"[green]✓ Environment exported to {output_path}[/green]")
        console.print(f"[cyan]Run: source {output_path}[/cyan]")
    
    sys.exit(0)


def cmd_config_show(args):
    """Show current COMPASS configuration."""
    verbose = not args.silent
    env_checker = EnvironmentChecker(verbose=verbose)
    env_checker.show_status()
    sys.exit(0)


def cmd_config_init(args):
    """Create a default configuration file."""
    verbose = not args.silent
    
    output_path = Path(args.output)
    
    # Create default configuration
    default_config = {
        'paths': {
            'compass_root': str(Path.home() / 'compass'),
            'shesha_root': str(Path.home() / 'compass' / 'shesha'),
            'compass_install_root': str(Path.home() / 'compass' / 'local'),
            'conda_root': str(Path.home() / 'miniconda3'),
        },
        'cuda': {
            'root': '/usr/local/cuda',
        },
        'build': {
            'type': 'Release',
            'parallel_jobs': 'auto',
        },
        'simulation': {
            'default_devices': '0',
            'default_iterations': 1000,
        }
    }
    
    import yaml
    with open(output_path, 'w') as f:
        f.write("# COMPASS Configuration File\n")
        f.write("# Edit this file and run: compass config export\n\n")
        yaml.dump(default_config, f, default_flow_style=False, sort_keys=False)
    
    if not args.silent:
        console.print(f"[green]✓ Configuration file created: {output_path}[/green]")
        console.print(f"[cyan]Edit and then run: compass config export[/cyan]")
    
    sys.exit(0)


def cmd_sim_run(args):
    """Run COMPASS simulation."""
    verbose = not args.silent
    
    runner = SimulationRunner(verbose=verbose)
    
    # Collect any extra arguments passed after --
    script_args = getattr(args, 'script_args', None)
    
    success = runner.run(
        args.param_file,
        iterations=args.iterations,
        devices=args.devices,
        script_args=script_args,
    )
    
    sys.exit(0 if success else 1)


def cmd_sim_list(args):
    """List available example simulations."""
    verbose = not args.silent
    
    runner = SimulationRunner(verbose=verbose)
    runner.list_examples(directory=args.directory)
    
    sys.exit(0)


def cmd_sim_script(args):
    """Manage default simulation script."""
    verbose = not args.silent
    
    runner = SimulationRunner(verbose=verbose)
    
    if args.script_action == "set":
        if not args.script_name:
            console.print("[red]✗ Script name required for 'set' action[/red]")
            console.print("[yellow]Usage: compass sim script set <script_name>[/yellow]")
            sys.exit(1)
        success = runner.set_default_script(args.script_name, script_type="default")
        sys.exit(0 if success else 1)
    elif args.script_action == "show":
        runner.show_script_config()
        sys.exit(0)
    elif args.script_action == "list":
        runner.list_available_scripts()
        sys.exit(0)
    else:
        console.print("[red]✗ Unknown action[/red]")
        sys.exit(1)


def cmd_sim_gui(args):
    """Run COMPASS GUI simulation."""
    verbose = not args.silent
    
    runner = SimulationRunner(verbose=verbose)
    
    # Collect any extra arguments passed after --
    script_args = getattr(args, 'script_args', None)
    
    success = runner.gui(
        args.param_file,
        iterations=args.iterations,
        devices=args.devices,
        script_args=script_args,
    )
    
    sys.exit(0 if success else 1)


def cmd_sim_gui_script(args):
    """Manage default GUI script."""
    verbose = not args.silent
    
    runner = SimulationRunner(verbose=verbose)
    
    if args.script_action == "set":
        if not args.script_name:
            console.print("[red]✗ Script name required for 'set' action[/red]")
            console.print("[yellow]Usage: compass sim gui-script set <script_name>[/yellow]")
            sys.exit(1)
        success = runner.set_default_script(args.script_name, script_type="gui")
        sys.exit(0 if success else 1)
    elif args.script_action == "show":
        # Show GUI script configuration
        config = runner.load_config()
        default_gui_script = config.get("default_gui_script", runner.default_gui_script)
        script_path = runner.get_script_path(script_type="gui")
        
        console.print("\n[bold cyan]Current GUI Script Configuration[/bold cyan]\n")
        from rich.table import Table
        table = Table(title="GUI Script Settings")
        table.add_column("Setting", style="cyan")
        table.add_column("Value", style="green")
        
        if str(runner.shesha_root / "shesha" / "widgets") in str(script_path):
            display_name = script_path.name
        else:
            display_name = str(script_path)
        
        table.add_row("Default GUI Script", display_name)
        table.add_row("Full Path", str(script_path.absolute()))
        table.add_row("Exists", "✓" if script_path.exists() else "✗")
        
        console.print(table)
        sys.exit(0)
    elif args.script_action == "list":
        # List available GUI scripts
        widgets_dir = runner.shesha_root / "shesha" / "widgets"
        
        if not widgets_dir.exists():
            console.print(f"[red]✗ Widgets directory not found: {widgets_dir}[/red]")
            sys.exit(1)
        
        scripts = sorted([s for s in widgets_dir.glob("*.py") if s.name != "__init__.py"])
        
        if not scripts:
            console.print(f"[yellow]No GUI scripts found in {widgets_dir}[/yellow]")
            sys.exit(0)
        
        config = runner.load_config()
        default_gui_path = str(runner.get_script_path(script_type="gui").absolute())
        
        console.print(f"\n[bold cyan]Available GUI Scripts[/bold cyan]\n")
        
        from rich.table import Table
        table = Table(title=f"Scripts in {widgets_dir.relative_to(runner.compass_root)}")
        table.add_column("Script", style="cyan")
        table.add_column("Status", style="green")
        
        for script in scripts:
            status = "✓ [default]" if str(script.absolute()) == default_gui_path else ""
            table.add_row(script.name, status)
        
        console.print(table)
        console.print(f"\n[dim]Set default with: compass sim gui-script set <script_name_or_path>[/dim]")
        sys.exit(0)
    else:
        console.print("[red]✗ Unknown action[/red]")
        sys.exit(1)



def cmd_check(args):
    """Run comprehensive system check."""
    verbose = not args.silent
    env_checker = EnvironmentChecker(verbose=verbose)
    
    # If specific check requested
    if args.config:
        env_checker.show_status()
        return
    
    if args.builds:
        builder = Builder(verbose=verbose)
        builder.show_status()
        return
    
    # Otherwise, run full check
    console.print("[bold blue]🔍 COMPASS Comprehensive System Check[/bold blue]\n")
    
    # Configuration Check
    env_checker.show_status()
    console.print()
    
    # Check if COMPASS module is loaded
    lmod_manager = LmodManager(verbose=verbose)
    if lmod_manager.is_compass_module_loaded():
        console.print("[green]✓ COMPASS module is loaded[/green]")
    else:
        console.print("[yellow]⚠ COMPASS module is not loaded[/yellow]")
        # Check if modulefiles are available
        modulepath = os.environ.get("MODULEPATH", "")
        if lmod_manager.is_lmod_installed() and str(lmod_manager.modulefiles_root) in modulepath:
            console.print("  [cyan]Run: module load compass/local[/cyan]")
        else:
            console.print("  [cyan]Run: compass init[/cyan] to setup modulefiles first")
    console.print()
    
    # Build Status
    builder = Builder(verbose=verbose)
    builder.show_status()
    
    console.print("\n[bold green]✓ System check complete![/bold green]")


def cmd_init(args):
    """Initialize COMPASS modulefiles setup."""
    verbose = not args.silent
    
    lmod_manager = LmodManager(verbose=verbose)
    success = lmod_manager.setup()
    
    sys.exit(0 if success else 1)


def create_parser():
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        prog="compass",
        description="🔭 COMPASS - COMputing Platform for Adaptive optics SystemS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  compass --version                     Show version
  compass init                          Setup COMPASS modulefiles
  compass check                         Check system configuration
  compass build --clean                 Build with clean
  compass config export                 Export environment variables
  compass sim script show               Show current default script
  compass sim script set closed_loop.py Set default simulation script
  compass sim run parfile.py            Run simulation (interactive IPython)
  compass sim gui parfile.py            Run GUI simulation
  compass sim gui-script show           Show current GUI script
  compass sim list                      List parameter directories
  compass sim list MICADO               List MICADO parameter files
        """
    )
    
    # Global options
    parser.add_argument(
        "--version",
        action="version",
        version=f"COMPASS {__version__}"
    )
    parser.add_argument(
        "-s", "--silent",
        action="store_true",
        help="Suppress output (silent mode)"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Init command
    init_parser = subparsers.add_parser("init", help="Initialize COMPASS modulefiles")
    init_parser.set_defaults(func=cmd_init)
    
    # Build command
    build_parser = subparsers.add_parser("build", help="Build COMPASS components")
    build_parser.add_argument("--clean", action="store_true", help="Clean build directories before building")
    build_parser.add_argument("-c", "--components", nargs="*", default=[], 
                             help="Specific components to build (libcarma, libsutra, python_module)")
    build_parser.set_defaults(func=cmd_build)
    
    # Config subcommand
    config_parser = subparsers.add_parser("config", help="Configuration management")
    config_subparsers = config_parser.add_subparsers(dest="config_command")
    
    # Config export
    export_parser = config_subparsers.add_parser("export", help="Export environment to shell script")
    export_parser.add_argument("--output", default="compass_env.sh", help="Output file path")
    export_parser.set_defaults(func=cmd_config_export)
    
    # Config show
    show_parser = config_subparsers.add_parser("show", help="Show current configuration")
    show_parser.set_defaults(func=cmd_config_show)
    
    # Config init
    init_parser = config_subparsers.add_parser("init", help="Create default configuration file")
    init_parser.add_argument("--output", default="compass-config.yaml", help="Output configuration file")
    init_parser.set_defaults(func=cmd_config_init)
    
    # Sim subcommand
    sim_parser = subparsers.add_parser("sim", help="Simulation control")
    sim_subparsers = sim_parser.add_subparsers(dest="sim_command")
    
    # Sim run
    run_parser = sim_subparsers.add_parser("run", help="Run COMPASS simulation (interactive)")
    run_parser.add_argument("param_file", help="Path to parameter file")
    run_parser.add_argument("--iterations", type=int, help="Number of iterations")
    run_parser.add_argument("--devices", help="GPU devices (comma-separated)")
    run_parser.add_argument("script_args", nargs="*", help="Additional arguments to pass to the script")
    run_parser.set_defaults(func=cmd_sim_run)
    
    # Sim list
    list_parser = sim_subparsers.add_parser("list", help="List example simulations")
    list_parser.add_argument("directory", nargs="?", help="Specific directory to list (optional)")
    list_parser.set_defaults(func=cmd_sim_list)
    
    # Sim gui
    gui_parser = sim_subparsers.add_parser("gui", help="Run COMPASS GUI simulation")
    gui_parser.add_argument("param_file", help="Path to parameter file")
    gui_parser.add_argument("--iterations", type=int, help="Number of iterations")
    gui_parser.add_argument("--devices", help="GPU devices (comma-separated)")
    gui_parser.add_argument("script_args", nargs="*", help="Additional arguments to pass to the GUI script")
    gui_parser.set_defaults(func=cmd_sim_gui)
    
    # Sim script
    script_parser = sim_subparsers.add_parser("script", help="Manage default simulation script")
    script_parser.add_argument("script_action", choices=["set", "show", "list"], help="Action: set, show, or list")
    script_parser.add_argument("script_name", nargs="?", help="Script name (for 'set' action)")
    script_parser.set_defaults(func=cmd_sim_script)
    
    # Sim gui-script
    gui_script_parser = sim_subparsers.add_parser("gui-script", help="Manage default GUI script")
    gui_script_parser.add_argument("script_action", choices=["set", "show", "list"], help="Action: set, show, or list")
    gui_script_parser.add_argument("script_name", nargs="?", help="Script name (for 'set' action)")
    gui_script_parser.set_defaults(func=cmd_sim_gui_script)


    
    # Check command
    check_parser = subparsers.add_parser("check", help="Run system check")
    check_parser.add_argument("--config", action="store_true", help="Show only configuration")
    check_parser.add_argument("--builds", action="store_true", help="Show only build status")
    check_parser.set_defaults(func=cmd_check)
    
    return parser


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    try:
        if hasattr(args, "func"):
            args.func(args)
        else:
            parser.print_help()
            sys.exit(0)
    except KeyboardInterrupt:
        console.print("\n[yellow]Interrupted by user[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"[bold red]Error: {e}[/bold red]")
        if "--verbose" in sys.argv or "-v" in sys.argv:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
