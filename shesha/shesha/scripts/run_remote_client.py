#!/usr/bin/env python
"""
COMPASS GUI Remote Client

Lightweight client for connecting to a remote COMPASS supervisor server.
This script does not require full COMPASS installation (no CUDA, no sutra/carma).

Requirements:
    - PyQt6
    - pyqtgraph
    - pyzmq
    - numpy
    - matplotlib (for colormaps)
    - qtconsole (optional, for IPython console)

Usage:
    python run_remote_client.py [options]
    
Options:
    --server HOST         Server address (default: localhost)
    --cmd-port PORT       Command port (default: 5555)
    --tel-port PORT       Telemetry port (default: 5556)
    --auto-connect        Automatically connect on startup
    --help                Show this help message
"""

import sys
import os
import argparse

# Check for required GUI dependencies
try:
    from PyQt6.QtWidgets import QApplication, QMessageBox
except ImportError:
    print("ERROR: PyQt6 is required but not installed.")
    print("Install with: pip install PyQt6")
    sys.exit(1)

try:
    import pyqtgraph as pg
    _ = pg
except ImportError:
    print("ERROR: pyqtgraph is required but not installed.")
    print("Install with: pip install pyqtgraph")
    sys.exit(1)

try:
    import zmq
    _ = zmq
except ImportError:
    print("ERROR: pyzmq is required but not installed.")
    print("Install with: pip install pyzmq")
    sys.exit(1)

try:
    import numpy as np
    _ = np
except ImportError:
    print("ERROR: numpy is required but not installed.")
    print("Install with: pip install numpy")
    sys.exit(1)

# Add shesha directory to path for GUI imports
# This allows running from any directory
script_dir = os.path.dirname(os.path.abspath(__file__))
shesha_root = os.path.dirname(script_dir)
if shesha_root not in sys.path:
    sys.path.insert(0, shesha_root)

# Now we can import the GUI modules (these don't require COMPASS core)
try:
    from shesha.gui.main_window import CompassMainWindow
except ImportError as e:
    print(f"ERROR: Failed to import GUI modules: {e}")
    print("Make sure you're running from the correct directory.")
    sys.exit(1)


def check_optional_dependencies():
    """Check for optional dependencies and warn if missing"""
    warnings = []
    
    try:
        import matplotlib
        _ = matplotlib
    except ImportError:
        warnings.append("matplotlib (colormaps will be limited)")
    
    try:
        import qtconsole
        _ = qtconsole
    except ImportError:
        warnings.append("qtconsole (IPython console will be disabled)")
    
    if warnings:
        print("\nOptional dependencies missing:")
        for warning in warnings:
            print(f"  - {warning}")
        print("\nInstall all with: pip install matplotlib qtconsole ipykernel")
        print()


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="COMPASS GUI Remote Client - Connect to remote supervisor server",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Launch GUI (will prompt for connection)
    python run_remote_client.py
    
    # Connect to specific server
    python run_remote_client.py --server 192.168.1.100
    
    # Custom ports with auto-connect
    python run_remote_client.py --server gpu-server --cmd-port 6000 --tel-port 6001 --auto-connect
    
    # Connect to localhost
    python run_remote_client.py --server localhost --auto-connect

Notes:
    - This client does NOT require full COMPASS installation (no CUDA, no sutra/carma)
    - Only GUI dependencies are needed: PyQt6, pyqtgraph, pyzmq, numpy
    - Optional: matplotlib (colormaps), qtconsole (IPython console)
    - The server must be running before connecting
        """
    )
    
    parser.add_argument('--server', default='localhost', 
                       help='Server address (default: localhost)')
    parser.add_argument('--cmd-port', type=int, default=5555,
                       help='Command port (default: 5555)')
    parser.add_argument('--tel-port', type=int, default=5556,
                       help='Telemetry port (default: 5556)')
    parser.add_argument('--auto-connect', action='store_true',
                       help='Automatically connect on startup')
    
    args = parser.parse_args()
    
    # Check optional dependencies
    check_optional_dependencies()
    
    # Create Qt application
    app = QApplication(sys.argv)
    app.setApplicationName("COMPASS Remote Client")
    
    # Create main window
    window = CompassMainWindow()
    
    # Enable remote mode automatically
    window.remote_mode_checkbox.setChecked(True)
    
    # If auto-connect requested, attempt connection
    if args.auto_connect:
        print(f"\nAuto-connecting to {args.server}:{args.cmd_port}...")
        
        # Store connection info
        window.connection_info = {
            'server_address': args.server,
            'command_port': args.cmd_port,
            'telemetry_port': args.tel_port
        }
        
        # Import here to avoid requiring it at module level
        from shesha.gui.remote_supervisor_client import RemoteSupervisorClient
        
        # Create remote client
        window.remote_client = RemoteSupervisorClient(
            server_address=args.server,
            command_port=args.cmd_port,
            telemetry_port=args.tel_port
        )
        
        # Try to connect
        window.status_label.setText("Connecting...")
        app.processEvents()
        
        if window.remote_client.connect():
            print(f"✓ Successfully connected to {args.server}")
            
            # Connection successful
            window.param_label.setText(f"Connected to {args.server}")
            window.status_label.setText("Connected")
            window.connect_remote_btn.setText("Disconnect")
            
            # Use remote client as supervisor
            window.supervisor = window.remote_client
            window.config = window.remote_client.config
            
            # Enable controls
            window._enable_simulation_controls()
            
            # Setup telemetry receiver
            window._setup_remote_telemetry()
            
            # Update IPython console if available
            try:
                from qtconsole.rich_jupyter_widget import RichJupyterWidget
                _ = RichJupyterWidget
                if window.kernel_manager:
                    window._update_console_namespace()
            except ImportError:
                pass
            
            QMessageBox.information(
                window,
                "Connected",
                f"Successfully connected to remote supervisor at {args.server}"
            )
        else:
            print(f"✗ Failed to connect to {args.server}:{args.cmd_port}")
            QMessageBox.critical(
                window,
                "Connection Failed",
                f"Failed to connect to {args.server}:{args.cmd_port}\n\n"
                "Make sure the remote server is running:\n"
                f"  python shesha/run_remote_server.py params.py\n\n"
                "Check:\n"
                "  - Server is running\n"
                "  - IP address is correct\n"
                "  - Ports are not blocked by firewall\n"
                f"  - Try: telnet {args.server} {args.cmd_port}"
            )
            window.status_label.setText("Connection failed")
            window.remote_client = None
    else:
        print("\nGUI launched in remote mode.")
        print("Click 'Connect...' to connect to a server.")
    
    # Show window
    window.show()
    
    print("\n" + "="*60)
    print("COMPASS Remote Client Ready")
    print("="*60)
    if not args.auto_connect:
        print("\nTo connect:")
        print("  1. Click 'Connect...' button")
        print("  2. Enter server address and ports")
        print("  3. Click 'Connect'")
    print("\nTo exit: Close window or press Ctrl+C")
    print("="*60 + "\n")
    
    # Run application
    return app.exec()


if __name__ == '__main__':
    sys.exit(main())
