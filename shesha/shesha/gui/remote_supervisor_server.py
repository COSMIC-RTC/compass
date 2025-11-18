"""
Remote Supervisor Server for COMPASS GUI

This module provides a ZeroMQ-based server that wraps a COMPASS supervisor
and exposes it for remote control. It handles command execution and streams
telemetry data to connected clients.

Architecture:
- Command Socket (REP): Receives commands from clients and sends responses
- Telemetry Socket (PUB): Publishes telemetry data to subscribed clients
"""

import zmq
import pickle
import logging
import traceback
from typing import Any, Dict, Optional
import numpy as np

logger = logging.getLogger(__name__)


class RemoteSupervisorServer:
    """
    ZeroMQ server that wraps a COMPASS supervisor for remote access.
    
    Provides two communication channels:
    1. Command channel (REP socket): Synchronous request-reply for commands
    2. Telemetry channel (PUB socket): Asynchronous publish-subscribe for data streaming
    """
    
    def __init__(self, supervisor, command_port: int = 5555, telemetry_port: int = 5556, 
                 bind_address: str = "*", supervisor_thread=None):
        """
        Initialize the remote supervisor server.
        
        Args:
            supervisor: COMPASS supervisor instance to wrap
            command_port: Port for command REP socket (default: 5555)
            telemetry_port: Port for telemetry PUB socket (default: 5556)
            bind_address: Address to bind to (default: "*" for all interfaces)
            supervisor_thread: Optional SupervisorThread for non-blocking loop execution
        """
        self.supervisor = supervisor
        self.supervisor_thread = supervisor_thread
        self.command_port = command_port
        self.telemetry_port = telemetry_port
        self.bind_address = bind_address
        
        # ZeroMQ context and sockets
        self.context = zmq.Context()
        self.command_socket: Optional[zmq.Socket] = None
        self.telemetry_socket: Optional[zmq.Socket] = None
        
        # Server state
        self.running = False
        
        # Detect two-stages mode
        self.is_two_stages = hasattr(supervisor, 'first_stage') and hasattr(supervisor, 'second_stage')
        
        logger.info(f"RemoteSupervisorServer initialized (cmd:{command_port}, tel:{telemetry_port}, two_stages:{self.is_two_stages})")
    
    def _get_active_supervisor(self):
        """Get the active supervisor (second stage for two-stages, or main supervisor)"""
        if self.is_two_stages:
            return self.supervisor.second_stage
        return self.supervisor
    
    def _get_config(self):
        """Get the config (from active supervisor or first stage for two-stages)"""
        if self.is_two_stages:
            return self.supervisor.first_stage.config
        return self.supervisor.config
    
    def start(self):
        """Start the server and bind sockets."""
        if self.running:
            logger.warning("Server already running")
            return
        
        # Create and bind command socket (REP pattern)
        self.command_socket = self.context.socket(zmq.REP)
        command_address = f"tcp://{self.bind_address}:{self.command_port}"
        self.command_socket.bind(command_address)
        logger.info(f"Command socket bound to {command_address}")
        
        # Create and bind telemetry socket (PUB pattern)
        self.telemetry_socket = self.context.socket(zmq.PUB)
        telemetry_address = f"tcp://{self.bind_address}:{self.telemetry_port}"
        self.telemetry_socket.bind(telemetry_address)
        logger.info(f"Telemetry socket bound to {telemetry_address}")
        
        self.running = True
        logger.info("RemoteSupervisorServer started")
    
    def stop(self):
        """Stop the server and close sockets."""
        if not self.running:
            return
        
        self.running = False
        
        if self.command_socket:
            self.command_socket.close()
            self.command_socket = None
        
        if self.telemetry_socket:
            self.telemetry_socket.close()
            self.telemetry_socket = None
        
        logger.info("RemoteSupervisorServer stopped")
    
    def __del__(self):
        """Cleanup on destruction."""
        self.stop()
        if hasattr(self, 'context'):
            self.context.term()
    
    def handle_command(self, timeout_ms: int = 100) -> bool:
        """
        Handle a single command request (non-blocking with timeout).
        
        Args:
            timeout_ms: Timeout in milliseconds for waiting for commands
            
        Returns:
            True if a command was processed, False if timeout
        """
        if not self.running or not self.command_socket:
            return False
        
        try:
            # Check for incoming commands (non-blocking)
            if self.command_socket.poll(timeout_ms, zmq.POLLIN):
                # Receive command
                message = self.command_socket.recv()
                request = pickle.loads(message)
                
                # Execute command
                response = self._execute_command(request)
                
                # Send response
                self.command_socket.send(pickle.dumps(response))
                return True
            
            return False
            
        except Exception as e:
            logger.error(f"Error handling command: {e}")
            logger.error(traceback.format_exc())
            # Send error response
            try:
                error_response = {
                    "status": "error",
                    "error": str(e),
                    "traceback": traceback.format_exc()
                }
                self.command_socket.send(pickle.dumps(error_response))
            except:
                pass
            return False
    
    def _execute_command(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute a command from a client request.
        
        Args:
            request: Dictionary with 'cmd' and optional 'args' keys
            
        Returns:
            Response dictionary with 'status' and optional 'result' or 'error'
        """
        cmd = request.get("cmd")
        args = request.get("args", {})
        
        logger.debug(f"Executing command: {cmd} with args: {args}")
        
        try:
            result = None
            
            # Thread-based loop control commands (non-blocking)
            if cmd == "start_thread":
                if self.supervisor_thread:
                    if not self.supervisor_thread.isRunning():
                        self.supervisor_thread.start()
                        result = "Thread started"
                    else:
                        result = "Thread already running"
                else:
                    return {"status": "error", "error": "No supervisor thread available"}
            
            elif cmd == "stop_thread":
                if self.supervisor_thread:
                    self.supervisor_thread.stop_loop()
                    result = "Stop requested"
                else:
                    return {"status": "error", "error": "No supervisor thread available"}
            
            elif cmd == "pause_thread":
                if self.supervisor_thread:
                    self.supervisor_thread.pause_loop()
                    result = "Thread paused"
                else:
                    return {"status": "error", "error": "No supervisor thread available"}
            
            elif cmd == "resume_thread":
                if self.supervisor_thread:
                    self.supervisor_thread.resume_loop()
                    result = "Thread resumed"
                else:
                    return {"status": "error", "error": "No supervisor thread available"}
            
            elif cmd == "step_thread":
                if self.supervisor_thread:
                    self.supervisor_thread.step_loop()
                    result = "Step requested"
                else:
                    return {"status": "error", "error": "No supervisor thread available"}
            
            # Direct supervisor commands (blocking - use with caution)
            elif cmd == "next":
                self.supervisor.next()
                result = "OK"
            
            elif cmd == "loop":
                iterations = args.get("iterations", 1)
                self.supervisor.loop(iterations)
                result = f"Completed {iterations} iterations"
            
            elif cmd == "reset":
                self.supervisor.reset()
                result = "Supervisor reset"
            
            # RTC commands
            elif cmd == "close_loop":
                if self.is_two_stages:
                    # Close loop on both stages
                    self.supervisor.first_stage.rtc.close_loop()
                    self.supervisor.second_stage.rtc.close_loop()
                else:
                    self.supervisor.rtc.close_loop()
                result = "Loop closed"
            
            elif cmd == "open_loop":
                if self.is_two_stages:
                    # Open loop on both stages
                    self.supervisor.first_stage.rtc.open_loop()
                    self.supervisor.second_stage.rtc.open_loop()
                else:
                    self.supervisor.rtc.open_loop()
                result = "Loop opened"
            
            # Target commands
            elif cmd == "reset_strehl":
                if self.is_two_stages:
                    # Reset both stages using TwoStagesManager method
                    self.supervisor.reset_exposure()
                else:
                    for target in self.supervisor.target:
                        target.reset_strehl()
                result = "Strehl reset on all targets"
            
            # Atmosphere commands
            elif cmd == "enable_atmos":
                enabled = args.get("enabled", True)
                if self.is_two_stages:
                    # Only enable/disable first stage atmos
                    # Second stage atmos is always disabled by TwoStagesManager
                    self.supervisor.first_stage.atmos.enable_atmos(enabled)
                else:
                    self.supervisor.atmos.enable_atmos(enabled)
                result = f"Atmosphere {'enabled' if enabled else 'disabled'}"
            
            # Two-stages stage selection
            elif cmd == "set_displayed_stage":
                stage_index = args.get("stage_index", 0)
                if self.supervisor_thread:
                    self.supervisor_thread.set_displayed_stage(stage_index)
                    stage_name = "Second Stage" if stage_index == 0 else "First Stage"
                    result = f"Now displaying {stage_name}"
                else:
                    result = "Stage selection only available with supervisor thread"
            
            elif cmd == "get_stage_config":
                if not self.is_two_stages:
                    return {
                        "status": "error",
                        "error": "Not in two-stages mode"
                    }
                
                stage_index = args.get("stage_index", 0)
                if stage_index == 0:
                    stage_supervisor = self.supervisor.second_stage
                else:
                    stage_supervisor = self.supervisor.first_stage
                
                stage_config = stage_supervisor.config
                config_info = {
                    "n_targets": len(stage_config.p_targets) if stage_config.p_targets else 0,
                    "n_wfs": len(stage_config.p_wfss) if stage_config.p_wfss else 0,
                    "n_dms": len(stage_config.p_dms) if stage_config.p_dms else 0,
                    "n_coronos": len(stage_config.p_coronos) if stage_config.p_coronos else 0,
                }
                result = config_info
            
            # Generic attribute access
            elif cmd == "get_attr":
                obj_path = args.get("obj_path")  # e.g., "rtc.delay"
                result = self._get_nested_attr(self.supervisor, obj_path)
            
            elif cmd == "set_attr":
                obj_path = args.get("obj_path")
                value = args.get("value")
                self._set_nested_attr(self.supervisor, obj_path, value)
                result = f"Set {obj_path} = {value}"
            
            elif cmd == "call_method":
                obj_path = args.get("obj_path")  # e.g., "rtc.close_loop"
                method_args = args.get("method_args", [])
                method_kwargs = args.get("method_kwargs", {})
                obj, method_name = self._get_nested_obj_and_method(self.supervisor, obj_path)
                method = getattr(obj, method_name)
                result = method(*method_args, **method_kwargs)
            
            # Configuration query
            elif cmd == "get_config":
                # Don't send the full config - it contains unpicklable objects
                # Send serializable info that clients need for UI setup
                active_supervisor = self._get_active_supervisor()
                config = self._get_config()
                
                config_info = {
                    "iter": self.supervisor.iter if hasattr(self.supervisor, 'iter') else 0,
                    "has_atmos": active_supervisor.atmos is not None,
                    "n_targets": len(config.p_targets) if config.p_targets else 0,
                    "n_wfs": len(config.p_wfss) if config.p_wfss else 0,
                    "n_dms": len(config.p_dms) if config.p_dms else 0,
                    "n_coronos": len(config.p_coronos) if config.p_coronos else 0,
                    "is_two_stages": self.is_two_stages,
                }
                result = config_info
            
            else:
                return {
                    "status": "error",
                    "error": f"Unknown command: {cmd}"
                }
            
            return {
                "status": "ok",
                "result": result
            }
            
        except Exception as e:
            logger.error(f"Error executing command {cmd}: {e}")
            return {
                "status": "error",
                "error": str(e),
                "traceback": traceback.format_exc()
            }
    
    def _get_nested_attr(self, obj: Any, path: str) -> Any:
        """Get nested attribute from object using dot notation."""
        parts = path.split('.')
        result = obj
        for part in parts:
            result = getattr(result, part)
        return result
    
    def _set_nested_attr(self, obj: Any, path: str, value: Any):
        """Set nested attribute on object using dot notation."""
        parts = path.split('.')
        target = obj
        for part in parts[:-1]:
            target = getattr(target, part)
        setattr(target, parts[-1], value)
    
    def _get_nested_obj_and_method(self, obj: Any, path: str):
        """Get object and method name from path like 'rtc.close_loop'."""
        parts = path.split('.')
        target = obj
        for part in parts[:-1]:
            target = getattr(target, part)
        return target, parts[-1]
    
    def publish_telemetry(self, topic: str, data: Any):
        """
        Publish telemetry data to all subscribed clients.
        
        Args:
            topic: Topic string for filtering (e.g., "image", "plot", "status")
            data: Data to publish (will be pickled)
        """
        if not self.running or not self.telemetry_socket:
            return
        
        try:
            # ZMQ PUB/SUB requires topic as bytes prefix
            message = pickle.dumps({
                "topic": topic,
                "data": data
            })
            self.telemetry_socket.send_multipart([topic.encode('utf-8'), message])
            
        except Exception as e:
            logger.error(f"Error publishing telemetry on topic '{topic}': {e}")
    
    def publish_image(self, category: str, title: str, image: np.ndarray, 
                     metadata: Optional[Dict] = None):
        """
        Publish an image for display.
        
        Args:
            category: Image category (e.g., "wfs", "dm", "target")
            title: Image title/identifier
            image: NumPy array containing image data
            metadata: Optional metadata dictionary
        """
        self.publish_telemetry("image", {
            "category": category,
            "title": title,
            "image": image,
            "metadata": metadata or {}
        })
    
    def publish_plot(self, category: str, title: str, x_data: np.ndarray, 
                    y_data: np.ndarray, metadata: Optional[Dict] = None):
        """
        Publish plot data.
        
        Args:
            category: Plot category
            title: Plot title/identifier
            x_data: X-axis data
            y_data: Y-axis data
            metadata: Optional metadata dictionary
        """
        self.publish_telemetry("plot", {
            "category": category,
            "title": title,
            "x_data": x_data,
            "y_data": y_data,
            "metadata": metadata or {}
        })
    
    def publish_status(self, status_dict: Dict[str, Any]):
        """
        Publish status information.
        
        Args:
            status_dict: Dictionary containing status information
        """
        self.publish_telemetry("status", status_dict)
