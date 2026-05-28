"""
Remote Supervisor Client for COMPASS GUI

This module provides a ZeroMQ-based client that connects to a RemoteSupervisorServer
and mimics the supervisor interface for transparent remote control.

The client can be used as a drop-in replacement for a local supervisor instance.
"""

import zmq
import pickle
import logging
import threading
from typing import Any, Dict, Optional, Callable
import queue

logger = logging.getLogger(__name__)


class RemoteSupervisorClient:
    """
    ZeroMQ client that connects to a remote supervisor server.
    
    Provides two communication channels:
    1. Command channel (REQ socket): Synchronous request-reply for commands
    2. Telemetry channel (SUB socket): Asynchronous receive for data streaming
    
    This class mimics the supervisor interface to allow transparent remote operation.
    """
    
    def __init__(self, server_address: str = "localhost", command_port: int = 5555, 
                 telemetry_port: int = 5556, timeout_ms: int = 5000):
        """
        Initialize the remote supervisor client.
        
        Args:
            server_address: Address of the remote server
            command_port: Port for command REQ socket (default: 5555)
            telemetry_port: Port for telemetry SUB socket (default: 5556)
            timeout_ms: Timeout for command requests in milliseconds (default: 5000)
        """
        self.server_address = server_address
        self.command_port = command_port
        self.telemetry_port = telemetry_port
        self.timeout_ms = timeout_ms
        
        # ZeroMQ context and sockets
        self.context = zmq.Context()
        self.command_socket: Optional[zmq.Socket] = None
        self.telemetry_socket: Optional[zmq.Socket] = None
        
        # Connection state
        self.connected = False
        
        # Lock to serialise command socket usage across threads
        self._cmd_lock = threading.Lock()
        
        # Telemetry handling
        self._telemetry_thread: Optional[threading.Thread] = None
        self._telemetry_running = False
        self._telemetry_callbacks: Dict[str, list] = {}
        self._telemetry_queue = queue.Queue()
        
        # Mock attributes to mimic supervisor interface
        self.rtc = RemoteRTCProxy(self)
        self.target = RemoteTargetListProxy(self)
        self.atmos = RemoteAtmosProxy(self)
        
        # Cached configuration
        self.iter = 0
        self.config = None
        self.config_info = {}  # Serializable config info from server
        
        logger.info(f"RemoteSupervisorClient initialized (server:{server_address}:{command_port})")
    
    def connect(self) -> bool:
        """
        Connect to the remote server.
        
        Returns:
            True if connection successful, False otherwise
        """
        if self.connected:
            logger.warning("Already connected")
            return True
        
        try:
            # Create and connect command socket (REQ pattern)
            self.command_socket = self.context.socket(zmq.REQ)
            self.command_socket.setsockopt(zmq.RCVTIMEO, self.timeout_ms)
            self.command_socket.setsockopt(zmq.SNDTIMEO, self.timeout_ms)
            command_address = f"tcp://{self.server_address}:{self.command_port}"
            self.command_socket.connect(command_address)
            logger.info(f"Command socket connected to {command_address}")
            
            # Create and connect telemetry socket (SUB pattern)
            self.telemetry_socket = self.context.socket(zmq.SUB)
            telemetry_address = f"tcp://{self.server_address}:{self.telemetry_port}"
            self.telemetry_socket.connect(telemetry_address)
            # Subscribe to all topics
            self.telemetry_socket.setsockopt(zmq.SUBSCRIBE, b"")
            logger.info(f"Telemetry socket connected to {telemetry_address}")
            
            # Set connected flag BEFORE testing so _send_command works
            self.connected = True
            
            # Test connection with a config query
            response = self._send_command("get_config", {})
            if response["status"] == "ok":
                result = response["result"]
                self.iter = result.get("iter", 0)
                
                # Store config info for UI setup
                self.config_info = result
                
                # Create a minimal mock config object for compatibility
                class MockConfig:
                    def __init__(self, config_info):
                        self.p_atmos = "mock" if config_info.get("has_atmos") else None
                        
                        # Create mock target list
                        n_targets = config_info.get("n_targets", 0)
                        self.p_targets = [f"target_{i}" for i in range(n_targets)] if n_targets > 0 else None
                        
                        # Create mock wfs list
                        n_wfs = config_info.get("n_wfs", 0)
                        self.p_wfss = [f"wfs_{i}" for i in range(n_wfs)] if n_wfs > 0 else None
                        
                        # Create mock dm list
                        n_dms = config_info.get("n_dms", 0)
                        self.p_dms = [f"dm_{i}" for i in range(n_dms)] if n_dms > 0 else None
                        
                        # Create mock corono list
                        n_coronos = config_info.get("n_coronos", 0)
                        self.p_coronos = [f"corono_{i}" for i in range(n_coronos)] if n_coronos > 0 else None
                
                self.config = MockConfig(self.config_info)
                
                # Start telemetry receiver thread
                self._start_telemetry_thread()
                
                logger.info("Successfully connected to remote supervisor")
                return True
            else:
                logger.error(f"Connection test failed: {response.get('error')}")
                self.connected = False
                self.disconnect()
                return False
                
        except Exception as e:
            logger.error(f"Failed to connect: {e}")
            self.disconnect()
            return False
    
    def disconnect(self):
        """Disconnect from the remote server."""
        if not self.connected:
            return
        
        # Stop telemetry thread
        self._stop_telemetry_thread()
        
        # Close sockets
        with self._cmd_lock:
            if self.command_socket:
                self.command_socket.close()
                self.command_socket = None
        
        if self.telemetry_socket:
            self.telemetry_socket.close()
            self.telemetry_socket = None
        
        self.connected = False
        logger.info("Disconnected from remote supervisor")
    
    def __del__(self):
        """Cleanup on destruction."""
        self.disconnect()
        if hasattr(self, 'context'):
            self.context.term()
    
    def _send_command(self, cmd: str, args: Dict[str, Any]) -> Dict[str, Any]:
        """
        Send a command to the server and wait for response.

        .. warning::
            Data is serialised with :mod:`pickle`.  Only connect to trusted
            servers on secured networks — a malicious server can execute
            arbitrary code on the client via crafted pickle payloads.

        Args:
            cmd: Command name
            args: Command arguments dictionary

        Returns:
            Response dictionary with 'status' and optional 'result' or 'error'
        """
        if not self.connected or not self.command_socket:
            return {
                "status": "error",
                "error": "Not connected to server"
            }

        with self._cmd_lock:
            try:
                request = {"cmd": cmd, "args": args}
                self.command_socket.send(pickle.dumps(request))

                # Wait for response
                message = self.command_socket.recv()
                response = pickle.loads(message)  # nosec B301 – trusted internal channel

                return response

            except zmq.Again:
                # ZMQ REQ socket is left in a broken state after a timeout.
                # The socket MUST be closed and recreated before the next send().
                logger.error(f"Command '{cmd}' timeout after {self.timeout_ms}ms — reconnecting socket")
                self._reconnect_command_socket()
                return {
                    "status": "error",
                    "error": f"Timeout after {self.timeout_ms}ms"
                }
            except zmq.ZMQError as e:
                logger.error(f"ZMQ error sending command '{cmd}': {e} — reconnecting socket")
                self._reconnect_command_socket()
                return {
                    "status": "error",
                    "error": str(e)
                }
            except Exception as e:
                logger.error(f"Error sending command '{cmd}': {e}")
                return {
                    "status": "error",
                    "error": str(e)
                }

    def _reconnect_command_socket(self) -> None:
        """Close and recreate the command socket to recover from a broken REQ state."""
        if self.command_socket:
            try:
                self.command_socket.close(linger=0)
            except Exception:
                pass
            self.command_socket = None

        try:
            self.command_socket = self.context.socket(zmq.REQ)
            self.command_socket.setsockopt(zmq.RCVTIMEO, self.timeout_ms)
            self.command_socket.setsockopt(zmq.SNDTIMEO, self.timeout_ms)
            command_address = f"tcp://{self.server_address}:{self.command_port}"
            self.command_socket.connect(command_address)
            logger.info(f"Command socket reconnected to {command_address}")
        except Exception as e:
            logger.error(f"Failed to reconnect command socket: {e}")
            self.connected = False
    
    # Supervisor interface methods
    
    def next(self):
        """Execute one iteration (mimics supervisor.next())."""
        response = self._send_command("next", {})
        if response["status"] != "ok":
            raise RuntimeError(f"next() failed: {response.get('error')}")
        self.iter += 1
    
    def loop(self, iterations: int = 1):
        """Execute multiple iterations (mimics supervisor.loop())."""
        response = self._send_command("loop", {"iterations": iterations})
        if response["status"] != "ok":
            raise RuntimeError(f"loop() failed: {response.get('error')}")
        self.iter += iterations
    
    def get_attr(self, obj_path: str) -> Any:
        """Get attribute from remote supervisor using dot notation."""
        response = self._send_command("get_attr", {"obj_path": obj_path})
        if response["status"] != "ok":
            raise AttributeError(f"Failed to get {obj_path}: {response.get('error')}")
        return response["result"]
    
    def set_attr(self, obj_path: str, value: Any):
        """Set attribute on remote supervisor using dot notation."""
        response = self._send_command("set_attr", {"obj_path": obj_path, "value": value})
        if response["status"] != "ok":
            raise AttributeError(f"Failed to set {obj_path}: {response.get('error')}")
    
    def call_method(self, obj_path: str, *args, **kwargs) -> Any:
        """Call method on remote supervisor using dot notation."""
        response = self._send_command("call_method", {
            "obj_path": obj_path,
            "method_args": args,
            "method_kwargs": kwargs
        })
        if response["status"] != "ok":
            raise RuntimeError(f"Failed to call {obj_path}: {response.get('error')}")
        return response["result"]
    
    # Telemetry handling
    
    def _start_telemetry_thread(self):
        """Start background thread for receiving telemetry."""
        if self._telemetry_running:
            return
        
        self._telemetry_running = True
        self._telemetry_thread = threading.Thread(target=self._telemetry_receiver, daemon=True)
        self._telemetry_thread.start()
        logger.info("Telemetry receiver thread started")
    
    def _stop_telemetry_thread(self):
        """Stop telemetry receiver thread."""
        if not self._telemetry_running:
            return
        
        self._telemetry_running = False
        if self._telemetry_thread:
            self._telemetry_thread.join(timeout=1.0)
            self._telemetry_thread = None
        logger.info("Telemetry receiver thread stopped")
    
    def _telemetry_receiver(self):
        """Background thread that receives telemetry from server."""
        while self._telemetry_running:
            try:
                if self.telemetry_socket and self.telemetry_socket.poll(100, zmq.POLLIN):
                    # Receive multipart message: [topic, data]
                    parts = self.telemetry_socket.recv_multipart()
                    if len(parts) >= 2:
                        topic = parts[0].decode('utf-8')
                        message = pickle.loads(parts[1])
                        
                        # Put in queue for processing in main thread
                        self._telemetry_queue.put((topic, message))
                        
                        # Also call registered callbacks
                        if topic in self._telemetry_callbacks:
                            for callback in self._telemetry_callbacks[topic]:
                                try:
                                    callback(message)
                                except Exception as e:
                                    logger.error(f"Error in telemetry callback: {e}")
                
            except Exception as e:
                if self._telemetry_running:
                    logger.error(f"Error receiving telemetry: {e}")
    
    def register_telemetry_callback(self, topic: str, callback: Callable):
        """
        Register a callback for telemetry on a specific topic.
        
        Args:
            topic: Topic to subscribe to (e.g., "image", "plot", "status")
            callback: Callable that takes the telemetry data as argument
        """
        if topic not in self._telemetry_callbacks:
            self._telemetry_callbacks[topic] = []
        self._telemetry_callbacks[topic].append(callback)
    
    def unregister_telemetry_callback(self, topic: str, callback: Callable):
        """Unregister a telemetry callback."""
        if topic in self._telemetry_callbacks:
            try:
                self._telemetry_callbacks[topic].remove(callback)
            except ValueError:
                pass
    
    def get_telemetry(self, block: bool = False, timeout: Optional[float] = None):
        """
        Get telemetry from queue.
        
        Args:
            block: Whether to block waiting for telemetry
            timeout: Timeout for blocking (None = infinite)
            
        Returns:
            Tuple of (topic, message) or None if no telemetry available
        """
        try:
            return self._telemetry_queue.get(block=block, timeout=timeout)
        except queue.Empty:
            return None


class RemoteRTCProxy:
    """Proxy for remote RTC object."""
    
    def __init__(self, client: RemoteSupervisorClient):
        self._client = client
    
    def close_loop(self):
        """Close the AO loop."""
        response = self._client._send_command("close_loop", {})
        if response["status"] != "ok":
            raise RuntimeError(f"close_loop failed: {response.get('error')}")
    
    def open_loop(self):
        """Open the AO loop."""
        response = self._client._send_command("open_loop", {})
        if response["status"] != "ok":
            raise RuntimeError(f"open_loop failed: {response.get('error')}")
    
    def __getattr__(self, name: str):
        """Forward attribute access to remote supervisor."""
        return self._client.get_attr(f"rtc.{name}")
    
    def __setattr__(self, name: str, value: Any):
        """Forward attribute setting to remote supervisor."""
        if name == "_client":
            object.__setattr__(self, name, value)
        else:
            self._client.set_attr(f"rtc.{name}", value)


class RemoteTargetListProxy:
    """Proxy for remote target list."""
    
    def __init__(self, client: RemoteSupervisorClient):
        self._client = client
    
    def reset_strehl(self):
        """Reset Strehl on all targets."""
        response = self._client._send_command("reset_strehl", {})
        if response["status"] != "ok":
            raise RuntimeError(f"reset_strehl failed: {response.get('error')}")
    
    def __iter__(self):
        """Iterate over targets - not fully supported in remote mode."""
        logger.warning("Target iteration not fully supported in remote mode")
        return iter([])
    
    def __getitem__(self, index: int):
        """Access target by index - limited support in remote mode."""
        logger.warning("Direct target access limited in remote mode")
        return None


class RemoteAtmosProxy:
    """Proxy for remote atmosphere object."""
    
    def __init__(self, client: RemoteSupervisorClient):
        self._client = client
    
    def enable_atmos(self, enabled: bool = True):
        """Enable or disable atmosphere."""
        response = self._client._send_command("enable_atmos", {"enabled": enabled})
        if response["status"] != "ok":
            raise RuntimeError(f"enable_atmos failed: {response.get('error')}")
    
    def __getattr__(self, name: str):
        """Forward attribute access to remote supervisor."""
        return self._client.get_attr(f"atmos.{name}")
    
    def __setattr__(self, name: str, value: Any):
        """Forward attribute setting to remote supervisor."""
        if name == "_client":
            object.__setattr__(self, name, value)
        else:
            self._client.set_attr(f"atmos.{name}", value)
