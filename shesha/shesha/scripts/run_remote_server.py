#!/usr/bin/env python
"""
Run a COMPASS supervisor as a remote server

This script initializes a COMPASS supervisor from a parameter file
and runs it as a remote server that can be controlled via ZeroMQ.

Usage:
    # Single stage
    python run_remote_server.py <parameter_file> [options]
    
    # Two-stages mode
    python run_remote_server.py <param_file1> <param_file2> <frequency_ratio> [options]
    
Options:
    --host HOST           Host to bind to (default: "*" for all interfaces)
    --cmd-port PORT       Command port (default: 5555)
    --tel-port PORT       Telemetry port (default: 5556)
    --auto-start          Automatically start the loop on initialization
    --help                Show this help message
"""

import sys
import os
import argparse
import signal
import logging
import time

# Add parent directory to path to import shesha modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from shesha.config import ParamConfig
from shesha.supervisor.compassSupervisor import CompassSupervisor
from shesha.gui.remote_supervisor_server import RemoteSupervisorServer
from shesha.gui.supervisor_thread import SupervisorThread

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class RemoteServerRunner:
    """Manages the remote supervisor server"""
    
    def __init__(self, param_file, param_file2=None, frequency_ratio=1, host="*", cmd_port=5555, tel_port=5556, auto_start=False):
        self.param_file = param_file
        self.param_file2 = param_file2
        self.frequency_ratio = frequency_ratio
        self.host = host
        self.cmd_port = cmd_port
        self.tel_port = tel_port
        self.auto_start = auto_start
        
        self.config = None
        self.supervisor = None
        self.server = None
        self.thread = None
        self.running = False
        self.is_two_stages = param_file2 is not None
        self.first_stage = None
        self.second_stage = None
        
    def initialize(self):
        """Initialize supervisor and server"""
        if self.is_two_stages:
            logger.info(f"Loading two-stages parameters from: {self.param_file} and {self.param_file2}")
            
            # Load both configurations
            try:
                from shesha.supervisor.twoStagesManager import TwoStagesManager
                from shesha.supervisor.stageSupervisor import StageSupervisor
                
                config1 = ParamConfig(self.param_file)
                config2 = ParamConfig(self.param_file2)
            except Exception as e:
                logger.error(f"Failed to load parameter files: {e}")
                return False
            
            logger.info("Initializing two-stages supervisor...")
            
            # Create both supervisors
            try:
                self.first_stage = StageSupervisor(config1)
                self.second_stage = StageSupervisor(config2)
                self.supervisor = TwoStagesManager(self.first_stage, self.second_stage, self.frequency_ratio)
                self.config = config1  # Use first stage config for telemetry setup
            except Exception as e:
                logger.error(f"Failed to initialize two-stages supervisor: {e}")
                return False
            
            logger.info(f"Two-stages supervisor initialized (frequency ratio: {self.frequency_ratio})")
        else:
            logger.info(f"Loading parameters from: {self.param_file}")
            
            # Load configuration
            try:
                self.config = ParamConfig(self.param_file)
            except Exception as e:
                logger.error(f"Failed to load parameter file: {e}")
                return False
            
            logger.info("Initializing COMPASS supervisor...")
            
            # Create supervisor
            try:
                from shesha.supervisor.compassSupervisor import CompassSupervisor
                self.supervisor = CompassSupervisor(self.config)
            except Exception as e:
                logger.error(f"Failed to initialize supervisor: {e}")
                return False
            
            logger.info("Supervisor initialized successfully")
        
        # Create supervisor thread first
        self.thread = SupervisorThread(self.supervisor)
        
        # Create remote server with thread reference
        self.server = RemoteSupervisorServer(
            self.supervisor,
            command_port=self.cmd_port,
            telemetry_port=self.tel_port,
            bind_address=self.host,
            supervisor_thread=self.thread
        )
        
        # Set remote server on thread for telemetry publishing
        self.thread.set_remote_server(self.server)
        
        # Connect signals for logging
        self.thread.error_occurred.connect(self._on_error)
        self.thread.status_changed.connect(self._on_status_changed)
        self.thread.iteration_done.connect(self._on_iteration)
        
        # Request all telemetry for full data streaming
        self._request_all_telemetry()
        
        logger.info("Remote server initialized")
        return True
    
    def _request_all_telemetry(self):
        """Request all available telemetry data"""
        # In two-stages mode, request telemetry from second stage (by default)
        active_supervisor = self.second_stage if self.is_two_stages else self.supervisor
        
        # Atmosphere
        if active_supervisor.atmos is not None:
            self.thread.request_telemetry('atmos_phase', enable=True)
        
        # Targets
        if active_supervisor.target is not None and self.config.p_targets:
            for tar_idx in range(len(self.config.p_targets)):
                self.thread.request_telemetry('target_psf_se', tar_idx, True)
                self.thread.request_telemetry('target_psf_le', tar_idx, True)
                self.thread.request_telemetry('target_phase', tar_idx, True)
        
        # WFS
        if active_supervisor.wfs is not None and self.config.p_wfss:
            for wfs_idx in range(len(self.config.p_wfss)):
                self.thread.request_telemetry('wfs_image', wfs_idx, True)
                self.thread.request_telemetry('wfs_phase', wfs_idx, True)
        
        # DM
        if active_supervisor.dms is not None and self.config.p_dms:
            for dm_idx in range(len(self.config.p_dms)):
                self.thread.request_telemetry('dm_shape', dm_idx, True)
        
        # Coronagraph
        if active_supervisor.corono is not None and self.config.p_coronos:
            for coro_idx in range(len(self.config.p_coronos)):
                self.thread.request_telemetry('corono_image', coro_idx, True)
    
    def start(self):
        """Start the server"""
        if not self.supervisor or not self.server:
            logger.error("Server not initialized. Call initialize() first.")
            return False
        
        logger.info("Starting remote server...")
        
        # Start ZeroMQ server
        self.server.start()
        
        logger.info(f"Server listening on:")
        logger.info(f"  Command port:   {self.cmd_port}")
        logger.info(f"  Telemetry port: {self.tel_port}")
        logger.info("Server is ready for connections")
        
        # Start supervisor thread if auto_start enabled
        if self.auto_start:
            logger.info("Auto-starting supervisor loop...")
            self.thread.start()
        
        self.running = True
        return True
    
    def run(self):
        """Main server loop"""
        logger.info("Server running. Press Ctrl+C to stop.")
        
        # Get QCoreApplication if available
        try:
            from PyQt6.QtCore import QCoreApplication
            qt_app = QCoreApplication.instance()
        except:
            qt_app = None
        
        try:
            while self.running:
                # Handle commands (with timeout to allow checking running flag)
                # This MUST continue even when supervisor thread is running
                self.server.handle_command(timeout_ms=10)
                
                # Process Qt events if available
                # This allows the supervisor thread signals to be processed
                if qt_app:
                    qt_app.processEvents()
                
                # Small sleep to avoid busy waiting
                time.sleep(0.001)
                
        except KeyboardInterrupt:
            logger.info("Received interrupt signal")
        finally:
            self.stop()
    
    def stop(self):
        """Stop the server"""
        logger.info("Stopping remote server...")
        
        self.running = False
        
        # Stop supervisor thread if running
        if self.thread and self.thread.isRunning():
            logger.info("Stopping supervisor thread...")
            self.thread.stop_loop()
            self.thread.wait()
        
        # Stop server
        if self.server:
            self.server.stop()
        
        logger.info("Server stopped")
    
    def _on_error(self, error_msg):
        """Handle error from supervisor thread"""
        logger.error(f"Supervisor error: {error_msg}")
    
    def _on_status_changed(self, status):
        """Handle status change"""
        logger.info(f"Status: {status}")
    
    def _on_iteration(self, iter_num):
        """Handle iteration completion"""
        if iter_num % 100 == 0:  # Log every 100 iterations
            logger.info(f"Iteration: {iter_num}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Run COMPASS supervisor as a remote server",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run server on default ports
    python run_remote_server.py params.py
    
    # Run server on custom ports
    python run_remote_server.py params.py --cmd-port 6000 --tel-port 6001
    
    # Run server and auto-start the loop
    python run_remote_server.py params.py --auto-start
    
    # Bind to specific interface
    python run_remote_server.py params.py --host 192.168.1.100
        """
    )
    
    parser.add_argument('param_file', help='COMPASS parameter file (first stage for two-stages)')
    parser.add_argument('param_file2', nargs='?', default=None, help='Second parameter file (for two-stages mode)')
    parser.add_argument('frequency_ratio', nargs='?', type=int, default=1, help='Frequency ratio for two-stages mode (default: 1)')
    parser.add_argument('--host', default='*', help='Host to bind to (default: "*" for all)')
    parser.add_argument('--cmd-port', type=int, default=5555, help='Command port (default: 5555)')
    parser.add_argument('--tel-port', type=int, default=5556, help='Telemetry port (default: 5556)')
    parser.add_argument('--auto-start', action='store_true', help='Auto-start the loop')
    
    args = parser.parse_args()
    
    # Check parameter file(s) exist
    if not os.path.exists(args.param_file):
        logger.error(f"Parameter file not found: {args.param_file}")
        return 1
    
    if args.param_file2 and not os.path.exists(args.param_file2):
        logger.error(f"Second parameter file not found: {args.param_file2}")
        return 1
    
    # Create QCoreApplication for Qt event loop (required for SupervisorThread signals)
    try:
        from PyQt6.QtCore import QCoreApplication
        app = QCoreApplication(sys.argv)
        logger.info("Qt event loop initialized")
    except ImportError:
        logger.warning("PyQt6 not available - running without Qt event loop")
        app = None
    
    # Create and initialize server
    runner = RemoteServerRunner(
        args.param_file,
        param_file2=args.param_file2,
        frequency_ratio=args.frequency_ratio,
        host=args.host,
        cmd_port=args.cmd_port,
        tel_port=args.tel_port,
        auto_start=args.auto_start
    )
    
    if not runner.initialize():
        logger.error("Failed to initialize server")
        return 1
    
    if not runner.start():
        logger.error("Failed to start server")
        return 1
    
    # Setup signal handler for clean shutdown
    def signal_handler(sig, frame):
        logger.info("Received signal, stopping...")
        runner.stop()
        sys.exit(0)
    
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    # Run main loop
    runner.run()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
