#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
#
# COMPASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COMPASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COMPASS. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2024 COSMIC Team

"""
Supervisor worker thread for running the AO loop independently from the GUI
"""

from PyQt6.QtCore import QThread, pyqtSignal
import numpy as np
import time
from typing import Optional, TYPE_CHECKING
from queue import Queue, Empty

if TYPE_CHECKING:
    from .remote_supervisor_server import RemoteSupervisorServer


class SupervisorThread(QThread):
    """
    Worker thread that runs the COMPASS supervisor loop independently from the GUI.
    Communicates via signals for thread-safe updates and accepts commands via a queue.
    """
    
    # Signals for thread-safe communication with GUI
    iteration_done = pyqtSignal(int)  # Emitted after each iteration with iteration number
    strehl_updated = pyqtSignal(float, float)  # SE Strehl, LE Strehl
    framerate_updated = pyqtSignal(float)  # Current framerate
    
    # Display data signals - only emitted when requested
    atmos_phase_updated = pyqtSignal(np.ndarray)  # Atmospheric phase screen
    target_psf_se_updated = pyqtSignal(np.ndarray, int)  # SE PSF and target index
    target_psf_le_updated = pyqtSignal(np.ndarray, int)  # LE PSF and target index
    target_phase_updated = pyqtSignal(np.ndarray, int)  # Target phase and target index
    wfs_image_updated = pyqtSignal(np.ndarray, int)  # WFS image and WFS index
    wfs_phase_updated = pyqtSignal(np.ndarray, int)  # WFS phase and WFS index
    dm_shape_updated = pyqtSignal(np.ndarray, int)  # DM shape and DM index
    corono_image_updated = pyqtSignal(np.ndarray, int)  # Coronagraph image and index
    
    error_occurred = pyqtSignal(str)  # Error message
    loop_finished = pyqtSignal()  # Emitted when loop completes
    status_changed = pyqtSignal(str)  # Status message
    
    def __init__(self, supervisor=None, remote_server=None):
        """
        Initialize the supervisor thread.
        
        Args:
            supervisor: COMPASS supervisor instance (can be set later)
            remote_server: Optional RemoteSupervisorServer for publishing telemetry
        """
        super().__init__()
        self.supervisor = supervisor
        self.remote_server = remote_server
        self.command_queue = Queue()
        self._running = False
        self._paused = False
        self._stop_requested = False
        self._single_step = False
        
        # Performance tracking
        self._last_time = 0
        self._frame_count = 0
        self._monitoring_freq = 10  # Update GUI every N frames
        
        # Data emission control - only emit what's being displayed
        self.telemetry_requests = {
            'atmos_phase': False,
            'target_psf_se': set(),  # Set of target indices
            'target_psf_le': set(),  # Set of target indices
            'target_phase': set(),  # Set of target indices
            'wfs_image': set(),  # Set of WFS indices
            'wfs_phase': set(),  # Set of WFS indices
            'dm_shape': set(),  # Set of DM indices
            'corono_image': set(),  # Set of coronagraph indices
        }
        
        # Two-stages mode: which stage to display (0=second stage, 1=first stage)
        self.displayed_stage_index = 0
        
    def set_supervisor(self, supervisor):
        """Set or update the supervisor instance"""
        self.supervisor = supervisor
    
    def set_remote_server(self, remote_server):
        """Set or update the remote server for telemetry publishing"""
        self.remote_server = remote_server
    
    def set_displayed_stage(self, stage_index):
        """Set which stage to display in two-stages mode (0=second, 1=first)"""
        self.displayed_stage_index = stage_index
        
    def run(self):
        """Main loop running in the worker thread"""
        if self.supervisor is None:
            self.error_occurred.emit("No supervisor initialized")
            return
            
        self._running = True
        self._stop_requested = False
        self.status_changed.emit("Running")
        
        self._last_time = time.time()
        self._frame_count = 0
        
        try:
            while self._running and not self._stop_requested:
                # Process commands from the queue
                self._process_commands()
                
                # Wait if paused
                if self._paused and not self._single_step:
                    time.sleep(0.01)  # Small sleep to avoid busy waiting
                    continue
                    
                # Execute one iteration
                try:
                    self._execute_iteration()
                except Exception as e:
                    self.error_occurred.emit(f"Error in iteration: {str(e)}")
                    self._running = False
                    break
                    
                # Reset single step flag
                if self._single_step:
                    self._single_step = False
                    self._paused = True
                    self.status_changed.emit("Paused")
                
                # CRITICAL: Allow Qt event loop to process events
                # This prevents GUI from becoming unresponsive
                self.msleep(1)  # Sleep for 1ms to yield to GUI thread
                    
        except Exception as e:
            self.error_occurred.emit(f"Fatal error in loop: {str(e)}")
        finally:
            self._running = False
            self.loop_finished.emit()
            self.status_changed.emit("Stopped")
            
    def _execute_iteration(self):
        """Execute a single supervisor iteration and emit updates."""
        # Check if this is a TwoStagesManager
        is_two_stages = hasattr(self.supervisor, 'first_stage') and hasattr(self.supervisor, 'second_stage')
        
        if is_two_stages:
            # Two-stages mode: call manager's next method
            self.supervisor.next(do_control=True)
            # Get iteration from first stage
            iter_num = self.supervisor.first_stage.get_frame_counter()
        else:
            # Standard mode: run one iteration of the supervisor
            self.supervisor.next(compute_tar_psf=True)
            iter_num = self.supervisor.get_frame_counter()
        self.iteration_done.emit(iter_num)
        
        # Update GUI periodically to avoid overwhelming it
        # OR if this is a single step (force immediate update)
        self._frame_count += 1
        if self._frame_count >= self._monitoring_freq or self._single_step:
            self._emit_telemetry()
            self._frame_count = 0
            
    def _emit_telemetry(self):
        """Emit telemetry data to update GUI displays - only requested data"""
        try:
            # Check if this is a TwoStagesManager
            is_two_stages = hasattr(self.supervisor, 'first_stage') and hasattr(self.supervisor, 'second_stage')
            
            # Select which supervisor to use for telemetry
            if is_two_stages:
                # Use displayed_stage_index to determine which stage to show
                if self.displayed_stage_index == 0:
                    active_supervisor = self.supervisor.second_stage
                else:
                    active_supervisor = self.supervisor.first_stage
                iter_num = self.supervisor.first_stage.get_frame_counter()
            else:
                active_supervisor = self.supervisor
                iter_num = self.supervisor.get_frame_counter()
            
            # Calculate framerate
            current_time = time.time()
            elapsed = current_time - self._last_time
            if elapsed > 0:
                framerate = self._frame_count / elapsed
                self.framerate_updated.emit(framerate)
            self._last_time = current_time
            
            # Always emit Strehl ratio if target exists
            if active_supervisor.target is not None:
                strehl = active_supervisor.target.get_strehl(0)
                self.strehl_updated.emit(float(strehl[0]), float(strehl[1]))
                
                # Publish to remote server if available
                if self.remote_server:
                    self.remote_server.publish_status({
                        'strehl_se': float(strehl[0]),
                        'strehl_le': float(strehl[1]),
                        'iter': iter_num,
                        'framerate': framerate if 'framerate' in locals() else 0,
                        'two_stages': is_two_stages
                    })
            
            # Emit atmospheric phase if requested
            if self.telemetry_requests['atmos_phase'] and active_supervisor.atmos is not None:
                atmos_phase = active_supervisor.atmos.get_atmos_layer(0)
                if atmos_phase is not None:
                    self.atmos_phase_updated.emit(atmos_phase.copy())
                    if self.remote_server:
                        self.remote_server.publish_image('atmos', 'Phase Screen', atmos_phase)
            
            # Emit target data if requested
            if active_supervisor.target is not None:
                for tar_idx in self.telemetry_requests['target_psf_se']:
                    psf_se = active_supervisor.target.get_tar_image(tar_idx, expo_type='se')
                    if psf_se is not None:
                        self.target_psf_se_updated.emit(psf_se.copy(), tar_idx)
                        if self.remote_server:
                            self.remote_server.publish_image('target', f'PSF SE {tar_idx}', psf_se)
                
                for tar_idx in self.telemetry_requests['target_psf_le']:
                    psf_le = active_supervisor.target.get_tar_image(tar_idx, expo_type='le')
                    if psf_le is not None:
                        self.target_psf_le_updated.emit(psf_le.copy(), tar_idx)
                        if self.remote_server:
                            self.remote_server.publish_image('target', f'PSF LE {tar_idx}', psf_le)
                
                for tar_idx in self.telemetry_requests['target_phase']:
                    phase = active_supervisor.target.get_tar_phase(tar_idx)
                    if phase is not None:
                        self.target_phase_updated.emit(phase.copy(), tar_idx)
                        if self.remote_server:
                            self.remote_server.publish_image('target', f'Phase {tar_idx}', phase)
            
            # Emit WFS data if requested
            if active_supervisor.wfs is not None:
                for wfs_idx in self.telemetry_requests['wfs_image']:
                    wfs_image = active_supervisor.wfs.get_wfs_image(wfs_idx)
                    if wfs_image is not None:
                        self.wfs_image_updated.emit(wfs_image.copy(), wfs_idx)
                        if self.remote_server:
                            self.remote_server.publish_image('wfs', f'Image {wfs_idx}', wfs_image)
                
                for wfs_idx in self.telemetry_requests['wfs_phase']:
                    wfs_phase = active_supervisor.wfs.get_wfs_phase(wfs_idx)
                    if wfs_phase is not None:
                        self.wfs_phase_updated.emit(wfs_phase.copy(), wfs_idx)
                        if self.remote_server:
                            self.remote_server.publish_image('wfs', f'Phase {wfs_idx}', wfs_phase)
            
            # Emit DM data if requested
            if active_supervisor.dms is not None:
                for dm_idx in self.telemetry_requests['dm_shape']:
                    dm_shape = active_supervisor.dms.get_dm_shape(dm_idx)
                    if dm_shape is not None:
                        self.dm_shape_updated.emit(dm_shape.copy(), dm_idx)
                        if self.remote_server:
                            self.remote_server.publish_image('dm', f'Shape {dm_idx}', dm_shape)
            
            # Emit coronagraph data if requested
            if active_supervisor.corono is not None:
                for coro_idx in self.telemetry_requests['corono_image']:
                    coro_image = active_supervisor.corono.get_image(coro_idx)
                    if coro_image is not None:
                        self.corono_image_updated.emit(coro_image.copy(), coro_idx)
                        if self.remote_server:
                            self.remote_server.publish_image('corono', f'Image {coro_idx}', coro_image)
                    
        except Exception as e:
            # Don't stop the loop for telemetry errors
            self.error_occurred.emit(f"Telemetry error: {str(e)}")
            
    def _process_commands(self):
        """Process commands from the command queue"""
        try:
            # Process only a few commands per iteration to avoid blocking
            for _ in range(10):  # Limit to 10 commands per iteration
                if self.command_queue.empty():
                    break
                command = self.command_queue.get_nowait()
                self._execute_command(command)
        except Empty:
            pass
        except Exception as e:
            self.error_occurred.emit(f"Command processing error: {str(e)}")
            
    def _execute_command(self, command):
        """Execute a command on the supervisor"""
        cmd_type = command.get('type')
        
        if cmd_type == 'pause':
            self._paused = True
            self.status_changed.emit("Paused")
            
        elif cmd_type == 'resume':
            self._paused = False
            self._last_time = time.time()
            self._frame_count = 0
            self.status_changed.emit("Running")
            
        elif cmd_type == 'stop':
            self._stop_requested = True
            self._running = False
            
        elif cmd_type == 'step':
            self._single_step = True
            
        elif cmd_type == 'reset':
            self.supervisor.reset()
            self.status_changed.emit("Reset complete")
            
        elif cmd_type == 'set_monitoring_freq':
            self._monitoring_freq = command.get('value', 10)
            
        elif cmd_type == 'custom':
            # Execute custom commands on supervisor
            func = command.get('function')
            args = command.get('args', ())
            kwargs = command.get('kwargs', {})
            if func and hasattr(self.supervisor, func):
                getattr(self.supervisor, func)(*args, **kwargs)
                
    # Public methods to send commands (thread-safe)
    def pause_loop(self):
        """Pause the loop"""
        self.command_queue.put({'type': 'pause'})
        
    def resume_loop(self):
        """Resume the loop"""
        self.command_queue.put({'type': 'resume'})
        
    def stop_loop(self):
        """Stop the loop"""
        self.command_queue.put({'type': 'stop'})
        
    def step_loop(self):
        """Execute one step"""
        self.command_queue.put({'type': 'step'})
        
    def reset_simulation(self):
        """Reset the simulation"""
        self.command_queue.put({'type': 'reset'})
        
    def set_monitoring_frequency(self, freq: int):
        """Set how often GUI updates are sent (every N frames)"""
        self.command_queue.put({'type': 'set_monitoring_freq', 'value': freq})
        
    def request_telemetry(self, data_type: str, index: int = None, enable: bool = True):
        """
        Request specific telemetry data to be emitted.
        Only requested data will be emitted, improving performance.
        
        Args:
            data_type: Type of data ('atmos_phase', 'target_psf_se', 'target_psf_le', 
                       'target_phase', 'wfs_image', 'wfs_phase', 'dm_shape', 'corono_image')
            index: Index for indexed data (target, WFS, DM, corono), None for atmos_phase
            enable: True to enable emission, False to disable
        """
        if data_type == 'atmos_phase':
            self.telemetry_requests[data_type] = enable
        elif data_type in self.telemetry_requests:
            if enable:
                if index is not None:
                    self.telemetry_requests[data_type].add(index)
            else:
                if index is not None:
                    self.telemetry_requests[data_type].discard(index)
        
    def send_custom_command(self, function_name: str, *args, **kwargs):
        """Send a custom command to execute on the supervisor"""
        self.command_queue.put({
            'type': 'custom',
            'function': function_name,
            'args': args,
            'kwargs': kwargs
        })
