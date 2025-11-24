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
# Copyright (C) 2011-2025 COSMIC Team

"""
Main window for the COMPASS GUI application
"""

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QLabel, QFileDialog, QMessageBox, QGroupBox, QGridLayout,
    QSpinBox, QStatusBar, QTabWidget, QSplitter,
    QApplication, QCheckBox, QDialog, QComboBox
)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QAction

try:
    from qtconsole.rich_jupyter_widget import RichJupyterWidget
    from qtconsole.inprocess import QtInProcessKernelManager
    IPYTHON_AVAILABLE = True
except ImportError:
    IPYTHON_AVAILABLE = False
    print("Warning: qtconsole not available. IPython console will be disabled.")
    print("Install with: pip install qtconsole")

from .supervisor_thread import SupervisorThread
from .display_widgets import PlotWidget, DualPlotWidget, SelectableImageDisplayWidget
from .layout_manager import LayoutManager
from .connection_dialog import ConnectionDialog
from .remote_supervisor_client import RemoteSupervisorClient

class CompassMainWindow(QMainWindow):
    """Main window for COMPASS GUI"""
    
    def __init__(self):
        super().__init__()
        self.supervisor = None
        self.supervisor_thread = None
        self.config = None
        self.param_file = None
        
        # Two-stages mode support
        self.is_two_stages = False
        self.two_stages_manager = None
        self.first_stage_supervisor = None
        self.second_stage_supervisor = None
        
        # Remote mode support
        self.remote_mode = False
        self.remote_client = None
        self.connection_info = None
        
        # IPython console
        self.ipython_widget = None
        self.kernel_manager = None
        self.kernel_client = None
        
        # Timer to ensure GUI responsiveness
        self.gui_update_timer = QTimer()
        self.gui_update_timer.timeout.connect(self._process_gui_events)
        self.gui_update_timer.setInterval(10)  # Process events every 10ms
        
        self._init_ui()
        self._create_menu_bar()
        
    def _init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle("COMPASS - AO Simulation GUI")
        self.setGeometry(100, 100, 1600, 900)
        
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout(central_widget)
        
        # Top section: Control panel
        control_group = self._create_control_panel()
        main_layout.addWidget(control_group)
        
        # Middle section: Split between displays and plots
        splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # Left: Image displays (tabs)
        display_tabs = self._create_display_tabs()
        splitter.addWidget(display_tabs)
        
        # Right: Performance plots
        plot_widget = self._create_plot_panel()
        splitter.addWidget(plot_widget)
        
        splitter.setSizes([800, 800])
        main_layout.addWidget(splitter)
        
        # Bottom section: IPython console
        if IPYTHON_AVAILABLE:
            console_widget = self._create_ipython_console()
            main_layout.addWidget(console_widget)
        
        # Bottom: Status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready - Load a parameter file to begin")
        
    def _create_menu_bar(self):
        """Create the menu bar"""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu('&File')
        
        load_action = QAction('&Load Parameters...', self)
        load_action.setShortcut('Ctrl+O')
        load_action.triggered.connect(self._load_parameters)
        file_menu.addAction(load_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction('E&xit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Simulation menu
        sim_menu = menubar.addMenu('&Simulation')
        
        reset_action = QAction('&Reset', self)
        reset_action.setShortcut('Ctrl+R')
        reset_action.triggered.connect(self._reset_simulation)
        sim_menu.addAction(reset_action)
        
        # View menu
        view_menu = menubar.addMenu('&View')
        
        toggle_layout_action = QAction('&Custom Layout', self)
        toggle_layout_action.setShortcut('Ctrl+L')
        toggle_layout_action.setCheckable(True)
        toggle_layout_action.triggered.connect(self._toggle_layout_mode)
        view_menu.addAction(toggle_layout_action)
        self.toggle_layout_action = toggle_layout_action
        
        # Help menu
        help_menu = menubar.addMenu('&Help')
        
        about_action = QAction('&About', self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)
        
    def _create_control_panel(self):
        """Create the control panel with buttons and parameters"""
        group = QGroupBox("Control Panel")
        layout = QGridLayout()
        
        # Parameter file section
        self.param_label = QLabel("No parameter file loaded")
        layout.addWidget(QLabel("Parameter File:"), 0, 0)
        layout.addWidget(self.param_label, 0, 1, 1, 3)
        
        self.load_btn = QPushButton("Load Parameters")
        self.load_btn.clicked.connect(self._load_parameters)
        layout.addWidget(self.load_btn, 0, 4)
        
        self.init_btn = QPushButton("Initialize Supervisor")
        self.init_btn.clicked.connect(self._initialize_supervisor)
        self.init_btn.setEnabled(False)
        layout.addWidget(self.init_btn, 0, 5)
        
        # Remote mode controls
        self.remote_mode_checkbox = QCheckBox("Remote Mode")
        self.remote_mode_checkbox.stateChanged.connect(self._toggle_remote_mode)
        layout.addWidget(self.remote_mode_checkbox, 0, 6)
        
        self.connect_remote_btn = QPushButton("Connect...")
        self.connect_remote_btn.clicked.connect(self._connect_remote)
        self.connect_remote_btn.setEnabled(False)
        self.connect_remote_btn.setVisible(False)
        layout.addWidget(self.connect_remote_btn, 0, 7)
        
        # Loop control buttons
        self.start_btn = QPushButton("▶ Start Loop")
        self.start_btn.clicked.connect(self._start_loop)
        self.start_btn.setEnabled(False)
        self.start_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold;")
        layout.addWidget(self.start_btn, 1, 0)
        
        self.pause_btn = QPushButton("⏸ Pause")
        self.pause_btn.clicked.connect(self._pause_loop)
        self.pause_btn.setEnabled(False)
        layout.addWidget(self.pause_btn, 1, 1)
        
        self.resume_btn = QPushButton("▶ Resume")
        self.resume_btn.clicked.connect(self._resume_loop)
        self.resume_btn.setEnabled(False)
        layout.addWidget(self.resume_btn, 1, 2)
        
        self.step_btn = QPushButton("⏭ Step")
        self.step_btn.clicked.connect(self._step_loop)
        self.step_btn.setEnabled(False)
        layout.addWidget(self.step_btn, 1, 3)
        
        self.stop_btn = QPushButton("⏹ Stop")
        self.stop_btn.clicked.connect(self._stop_loop)
        self.stop_btn.setEnabled(False)
        self.stop_btn.setStyleSheet("background-color: #f44336; color: white; font-weight: bold;")
        layout.addWidget(self.stop_btn, 1, 4)
        
        self.reset_btn = QPushButton("↻ Reset")
        self.reset_btn.clicked.connect(self._reset_simulation)
        self.reset_btn.setEnabled(False)
        layout.addWidget(self.reset_btn, 1, 5)
        
        # Monitoring frequency control
        layout.addWidget(QLabel("GUI Update Freq (frames):"), 2, 0)
        self.monitoring_freq_spin = QSpinBox()
        self.monitoring_freq_spin.setMinimum(1)
        self.monitoring_freq_spin.setMaximum(1000)
        self.monitoring_freq_spin.setValue(50)
        self.monitoring_freq_spin.valueChanged.connect(self._update_monitoring_freq)
        layout.addWidget(self.monitoring_freq_spin, 2, 1)
        
        # Additional control buttons
        self.close_loop_btn = QPushButton("Loop Closed")
        self.close_loop_btn.setCheckable(True)
        self.close_loop_btn.clicked[bool].connect(self._toggle_loop_state)
        self.close_loop_btn.setEnabled(False)
        layout.addWidget(self.close_loop_btn, 2, 2)
        
        # Two-stages mode selector
        layout.addWidget(QLabel("Display Stage:"), 2, 3)
        self.stage_selector = QComboBox()
        self.stage_selector.addItem("Second Stage (Default)")
        self.stage_selector.addItem("First Stage")
        self.stage_selector.currentIndexChanged.connect(self._change_displayed_stage)
        self.stage_selector.setEnabled(False)
        self.stage_selector.setVisible(False)
        layout.addWidget(self.stage_selector, 2, 5)
        
        self.reset_strehl_btn = QPushButton("Reset Strehl")
        self.reset_strehl_btn.clicked.connect(self._reset_strehl)
        self.reset_strehl_btn.setEnabled(False)
        layout.addWidget(self.reset_strehl_btn, 2, 3)
        
        self.enable_atmos_btn = QPushButton("Atmos Enabled")
        self.enable_atmos_btn.setCheckable(True)
        self.enable_atmos_btn.setChecked(True)  # Atmosphere enabled by default
        self.enable_atmos_btn.clicked[bool].connect(self._toggle_atmosphere)
        self.enable_atmos_btn.setEnabled(False)
        layout.addWidget(self.enable_atmos_btn, 2, 4)
        
        # Status indicators
        layout.addWidget(QLabel("Status:"), 3, 0)
        self.status_label = QLabel("Not initialized")
        self.status_label.setStyleSheet("font-weight: bold;")
        layout.addWidget(self.status_label, 3, 1, 1, 2)  # Span 2 columns
        
        layout.addWidget(QLabel("Iteration:"), 3, 3)
        self.iter_label = QLabel("0")
        self.iter_label.setStyleSheet("font-weight: bold;")
        layout.addWidget(self.iter_label, 3, 4)
        
        group.setLayout(layout)
        return group
        
    def _create_display_tabs(self):
        """Create tabbed display for different component types"""
        tabs = QTabWidget()
        
        # Store reference to tabs widget
        self.display_tabs = tabs
        
        # Atmosphere display
        self.atmos_display = SelectableImageDisplayWidget("Atmosphere")
        self.atmos_display.selection_changed.connect(self._on_display_selection_changed)
        tabs.addTab(self.atmos_display, "Atmosphere")
        
        # Target display
        self.target_display = SelectableImageDisplayWidget("Target")
        self.target_display.selection_changed.connect(self._on_display_selection_changed)
        tabs.addTab(self.target_display, "Target")
        
        # WFS display
        self.wfs_display = SelectableImageDisplayWidget("WFS")
        self.wfs_display.selection_changed.connect(self._on_display_selection_changed)
        tabs.addTab(self.wfs_display, "WFS")
        
        # DM display
        self.dm_display = SelectableImageDisplayWidget("DM")
        self.dm_display.selection_changed.connect(self._on_display_selection_changed)
        tabs.addTab(self.dm_display, "DM")
        
        # Coronagraph display
        self.corono_display = SelectableImageDisplayWidget("Coronagraph")
        self.corono_display.selection_changed.connect(self._on_display_selection_changed)
        tabs.addTab(self.corono_display, "Corona")
        
        # Create custom layout manager
        self.layout_manager = LayoutManager()
        self.layout_manager.display_selection_changed.connect(self._on_custom_display_selection_changed)
        tabs.addTab(self.layout_manager, "⚙ Custom")
        
        return tabs
        
    def _create_plot_panel(self):
        """Create the plotting panel"""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        
        # Strehl ratio plot (dual: SE and LE)
        self.strehl_plot = DualPlotWidget(
            title="Strehl Ratio",
            ylabel="Strehl Ratio",
            legend1="Short Exposure",
            legend2="Long Exposure",
            max_points=1000
        )
        layout.addWidget(self.strehl_plot)
        
        # Framerate plot
        self.framerate_plot = PlotWidget(
            title="Framerate",
            ylabel="FPS",
            max_points=500
        )
        layout.addWidget(self.framerate_plot)
        
        return widget
    
    def _create_ipython_console(self):
        """Create an embedded IPython console for interactive control"""
        # Create a group box for the console
        console_group = QGroupBox("IPython Console")
        console_layout = QVBoxLayout()
        
        # Add show/hide checkbox
        header_layout = QHBoxLayout()
        self.console_visible_checkbox = QCheckBox("Show Console")
        self.console_visible_checkbox.setChecked(False)
        self.console_visible_checkbox.stateChanged.connect(self._toggle_console_visibility)
        header_layout.addWidget(self.console_visible_checkbox)
        header_layout.addStretch()
        
        info_label = QLabel("Access supervisor, config, and GUI via IPython. Try: supervisor.get_frame_counter()")
        info_label.setStyleSheet("color: gray; font-style: italic;")
        header_layout.addWidget(info_label)
        
        console_layout.addLayout(header_layout)
        
        # Create the IPython widget
        self.kernel_manager = QtInProcessKernelManager()
        self.kernel_manager.start_kernel()
        self.kernel_client = self.kernel_manager.client()
        self.kernel_client.start_channels()
        
        self.ipython_widget = RichJupyterWidget()
        self.ipython_widget.kernel_manager = self.kernel_manager
        self.ipython_widget.kernel_client = self.kernel_client
        self.ipython_widget.setMinimumHeight(200)
        self.ipython_widget.setMaximumHeight(400)
        
        # Push initial variables to the console namespace
        self.kernel_manager.kernel.shell.push({
            'gui': self,
            'np': __import__('numpy'),
        })
        
        # Initially hide the console
        self.ipython_widget.setVisible(False)
        
        console_layout.addWidget(self.ipython_widget)
        console_group.setLayout(console_layout)
        
        return console_group
    
    def _toggle_console_visibility(self, state):
        """Toggle IPython console visibility"""
        self.ipython_widget.setVisible(state == Qt.CheckState.Checked.value)
        
    def _update_console_namespace(self):
        """Update the IPython console namespace with current supervisor and config"""
        if self.ipython_widget and self.kernel_manager:
            namespace = {
                'gui': self,
                'supervisor': self.supervisor,
                'config': self.config,
                'np': __import__('numpy'),
            }
            
            # Add two-stages specific objects if in two-stages mode
            if self.is_two_stages:
                namespace['manager'] = self.two_stages_manager
                namespace['first_stage'] = self.first_stage_supervisor
                namespace['second_stage'] = self.second_stage_supervisor
            
            # Add supervisor_thread if available
            if self.supervisor_thread:
                namespace['supervisor_thread'] = self.supervisor_thread
            
            self.kernel_manager.kernel.shell.push(namespace)
            
            # Print welcome message in console
            if self.supervisor:
                if self.is_two_stages:
                    self.kernel_manager.kernel.shell.run_cell(
                        "print('\\n=== COMPASS GUI Console (Two-Stages Mode) ===')\n"
                        "print('Available objects:')\n"
                        "print('  - manager: TwoStagesManager instance')\n"
                        "print('  - first_stage: First stage StageSupervisor')\n"
                        "print('  - second_stage: Second stage StageSupervisor')\n"
                        "print('  - config: Parameter configuration (first stage)')\n"
                        "print('  - gui: Main window instance')\n"
                        "print('  - np: NumPy module')"
                    )
                else:
                    self.kernel_manager.kernel.shell.run_cell(
                        "print('\\n=== COMPASS GUI Console ===')\n"
                        "print('Available objects:')\n"
                        "print('  - supervisor: CompassSupervisor instance')\n"
                    "print('  - config: Parameter configuration')\n"
                    "print('  - supervisor_thread: Worker thread')\n"
                    "print('  - gui: Main window instance')\n"
                    "print('  - np: numpy')\n"
                    "print('\\nExample commands:')\n"
                    "print('  supervisor.get_frame_counter()')\n"
                    "print('  supervisor.rtc.get_err()')\n"
                    "print('  config.p_wfss[0].nxsub')\n"
                    "print('========================\\n')"
                )
        
    def _load_parameters(self):
        """Load a parameter file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Parameter File",
            "",
            "Python Files (*.py);;All Files (*)"
        )
        
        if file_path:
            try:
                from shesha.config import ParamConfig
                self.param_file = file_path
                self.config = ParamConfig(file_path)
                self.param_label.setText(file_path.split('/')[-1])
                self.init_btn.setEnabled(True)
                self.status_bar.showMessage(f"Loaded parameters: {file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load parameters:\n{str(e)}")
                self.param_file = None
                self.config = None
                
    def _initialize_supervisor(self):
        """Initialize the COMPASS supervisor"""
        if self.config is None:
            QMessageBox.warning(self, "Warning", "No parameters loaded!")
            return
            
        try:
            self.status_bar.showMessage("Initializing supervisor...")
            self.status_label.setText("Initializing...")
            
            # Check if this is a two-stages configuration
            if hasattr(self, 'param_file2') and self.param_file2 is not None:
                # Two-stages mode: use TwoStagesManager with StageSupervisor
                from shesha.supervisor.twoStagesManager import TwoStagesManager
                from shesha.supervisor.stageSupervisor import StageSupervisor
                from shesha.config import ParamConfig
                
                self.status_bar.showMessage("Initializing two-stages supervisor...")
                
                # Load both configurations
                config1 = ParamConfig(self.param_file)
                config2 = ParamConfig(self.param_file2)
                
                # Create both StageSupervisors (NOT CompassSupervisor!)
                supervisor1 = StageSupervisor(config1)
                supervisor2 = StageSupervisor(config2)
                
                # Create two-stages manager
                freq_ratio = getattr(self, 'frequency_ratio', 1)
                self.two_stages_manager = TwoStagesManager(supervisor1, supervisor2, freq_ratio)
                
                # Store references
                self.first_stage_supervisor = supervisor1
                self.second_stage_supervisor = supervisor2
                self.supervisor = self.two_stages_manager  # Main reference
                self.config = config1  # Use first stage config for display options
                self.is_two_stages = True
                
                self.status_bar.showMessage("Two-stages supervisor initialized successfully")
            else:
                # Single stage mode: use standard CompassSupervisor
                from shesha.supervisor.compassSupervisor import CompassSupervisor
                self.supervisor = CompassSupervisor(self.config)
                self.is_two_stages = False
            
            # Create and setup supervisor thread
            self.supervisor_thread = SupervisorThread(self.supervisor)
            self._connect_signals()
            
            # Populate display options based on config
            self._populate_display_options()
            
            # Enable controls
            self.start_btn.setEnabled(True)
            self.step_btn.setEnabled(True)
            self.reset_btn.setEnabled(True)
            self.close_loop_btn.setEnabled(True)
            self.close_loop_btn.setChecked(True)  # Start in closed-loop
            self.reset_strehl_btn.setEnabled(True)
            
            # Enable stage selector if two-stages mode
            if self.is_two_stages:
                self.stage_selector.setEnabled(True)
                self.stage_selector.setVisible(True)
            self.enable_atmos_btn.setEnabled(True)
            self.enable_atmos_btn.setChecked(True)
            self.init_btn.setEnabled(False)
            
            self.status_label.setText("Initialized")
            self.status_bar.showMessage("Supervisor initialized successfully")
            
            # Update IPython console namespace
            if IPYTHON_AVAILABLE:
                self._update_console_namespace()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to initialize supervisor:\n{str(e)}")
            self.supervisor = None
            self.supervisor_thread = None
            
    def _connect_signals(self):
        """Connect signals from supervisor thread to GUI slots"""
        self.supervisor_thread.iteration_done.connect(self._on_iteration_done)
        self.supervisor_thread.strehl_updated.connect(self._on_strehl_updated)
        self.supervisor_thread.framerate_updated.connect(self._on_framerate_updated)
        
        # Connect new selective data signals
        self.supervisor_thread.atmos_phase_updated.connect(self._on_atmos_phase_updated)
        self.supervisor_thread.target_psf_se_updated.connect(self._on_target_psf_se_updated)
        self.supervisor_thread.target_psf_le_updated.connect(self._on_target_psf_le_updated)
        self.supervisor_thread.target_phase_updated.connect(self._on_target_phase_updated)
        self.supervisor_thread.wfs_image_updated.connect(self._on_wfs_image_updated)
        self.supervisor_thread.wfs_phase_updated.connect(self._on_wfs_phase_updated)
        self.supervisor_thread.dm_shape_updated.connect(self._on_dm_shape_updated)
        self.supervisor_thread.corono_image_updated.connect(self._on_corono_image_updated)
        
        self.supervisor_thread.error_occurred.connect(self._on_error)
        self.supervisor_thread.loop_finished.connect(self._on_loop_finished)
        self.supervisor_thread.status_changed.connect(self._on_status_changed)
    
    def _populate_display_options(self, stage_config=None):
        """Populate dropdown menus based on supervisor configuration
        
        Args:
            stage_config: Optional specific stage config (for two-stages mode)
        """
        if not self.supervisor:
            return
        
        # Determine which config to use
        if stage_config is not None:
            config = stage_config
        elif self.config is not None:
            config = self.config
        else:
            return
        
        # Build options dictionary for all categories
        options_dict = {}
        
        # Atmosphere options
        atmos_options = []
        if config.p_atmos is not None:
            atmos_options.append(("Atmospheric Phase", "atmos_phase", 0))
        self.atmos_display.set_options(atmos_options)
        if atmos_options:
            options_dict['Atmosphere'] = atmos_options
        
        # Target options
        target_options = []
        if config.p_targets is not None:
            for i in range(len(config.p_targets)):
                target_options.append((f"Target {i} - PSF SE", "target_psf_se", i))
                target_options.append((f"Target {i} - PSF LE", "target_psf_le", i))
                target_options.append((f"Target {i} - Phase", "target_phase", i))
        self.target_display.set_options(target_options)
        if target_options:
            options_dict['Target'] = target_options
        
        # WFS options
        wfs_options = []
        if config.p_wfss is not None:
            for i in range(len(config.p_wfss)):
                wfs_options.append((f"WFS {i} - Image", "wfs_image", i))
                wfs_options.append((f"WFS {i} - Phase", "wfs_phase", i))
        self.wfs_display.set_options(wfs_options)
        if wfs_options:
            options_dict['WFS'] = wfs_options
        
        # DM options
        dm_options = []
        if config.p_dms is not None:
            for i in range(len(config.p_dms)):
                dm_options.append((f"DM {i} - Shape", "dm_shape", i))
        self.dm_display.set_options(dm_options)
        if dm_options:
            options_dict['DM'] = dm_options
        
        # Coronagraph options
        corono_options = []
        if config.p_coronos is not None:
            for i in range(len(config.p_coronos)):
                corono_options.append((f"Coronagraph {i} - Image", "corono_image", i))
        self.corono_display.set_options(corono_options)
        if corono_options:
            options_dict['Coronagraph'] = corono_options
        
        # Set options for custom layout manager
        self.layout_manager.set_available_options(options_dict)
    
    def _on_display_selection_changed(self, data_type, index):
        """Handle display selection change - update telemetry requests"""
        if self.remote_mode:
            # In remote mode, all telemetry is streamed - just update display
            return
        
        if not self.supervisor_thread:
            return
        
        # Disable all telemetry first
        self.supervisor_thread.request_telemetry('atmos_phase', enable=False)
        
        # Clear all indexed telemetry
        for dt in ['target_psf_se', 'target_psf_le', 'target_phase', 
                   'wfs_image', 'wfs_phase', 'dm_shape', 'corono_image']:
            self.supervisor_thread.telemetry_requests[dt].clear()
        
        # Enable only the currently selected displays
        for display in [self.atmos_display, self.target_display, self.wfs_display,
                       self.dm_display, self.corono_display]:
            current_type, current_idx = display.get_current_selection()
            if current_type:
                self.supervisor_thread.request_telemetry(current_type, current_idx, enable=True)
        
        # Enable displays in custom layout
        for display in self.layout_manager.get_displays():
            current_type, current_idx = display.get_current_selection()
            if current_type:
                self.supervisor_thread.request_telemetry(current_type, current_idx, enable=True)
    
    def _on_custom_display_selection_changed(self, data_type, index, display_widget):
        """Handle display selection change in custom layout - update telemetry requests"""
        if self.remote_mode:
            # In remote mode, all telemetry is streamed - just update display
            return
        
        if not self.supervisor_thread:
            return
        
        # Disable all telemetry first
        self.supervisor_thread.request_telemetry('atmos_phase', enable=False)
        
        # Clear all indexed telemetry
        for dt in ['target_psf_se', 'target_psf_le', 'target_phase', 
                   'wfs_image', 'wfs_phase', 'dm_shape', 'corono_image']:
            self.supervisor_thread.telemetry_requests[dt].clear()
        
        # Enable only the currently selected displays (tabs)
        for display in [self.atmos_display, self.target_display, self.wfs_display,
                       self.dm_display, self.corono_display]:
            current_type, current_idx = display.get_current_selection()
            if current_type:
                self.supervisor_thread.request_telemetry(current_type, current_idx, enable=True)
        
        # Enable displays in custom layout
        for display in self.layout_manager.get_displays():
            current_type, current_idx = display.get_current_selection()
            if current_type:
                self.supervisor_thread.request_telemetry(current_type, current_idx, enable=True)
    
    def _toggle_layout_mode(self, checked):
        """Toggle between tab view and custom layout"""
        if checked:
            # Switch to custom layout tab
            self.display_tabs.setCurrentWidget(self.layout_manager)
        
    def _start_loop(self):
        """Start the simulation loop"""
        if not self.remote_mode and self.supervisor_thread is None:
            QMessageBox.warning(self, "Warning", "Supervisor not initialized!")
            return
        
        if self.remote_mode and not self.remote_client:
            QMessageBox.warning(self, "Warning", "Not connected to remote server!")
            return
            
        # Update button states
        self.start_btn.setEnabled(False)
        self.pause_btn.setEnabled(True)
        self.stop_btn.setEnabled(True)
        self.step_btn.setEnabled(False)
        self.init_btn.setEnabled(False)
        self.load_btn.setEnabled(False)
        
        # Clear plots
        self.strehl_plot.clear()
        self.framerate_plot.clear()
        
        # Start GUI responsiveness timer
        self.gui_update_timer.start()
        
        if self.remote_mode:
            # In remote mode, send start_thread command to server
            try:
                # Start supervisor thread on remote server
                response = self.remote_client._send_command('start_thread', {})
                if response.get('status') != 'ok':
                    raise Exception(response.get('error', 'Unknown error'))
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to start remote loop:\n{str(e)}")
                # Revert button states
                self.start_btn.setEnabled(True)
                self.pause_btn.setEnabled(False)
                self.stop_btn.setEnabled(False)
                return
        else:
            # Local mode - start thread
            self.supervisor_thread.start()
        
    def _pause_loop(self):
        """Pause the simulation loop"""
        if self.remote_mode:
            try:
                response = self.remote_client._send_command('pause_thread', {})
                if response.get('status') == 'ok':
                    self.pause_btn.setEnabled(False)
                    self.resume_btn.setEnabled(True)
                    self.step_btn.setEnabled(True)
                else:
                    raise Exception(response.get('error', 'Unknown error'))
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to pause:\n{str(e)}")
        elif self.supervisor_thread:
            self.supervisor_thread.pause_loop()
            self.pause_btn.setEnabled(False)
            self.resume_btn.setEnabled(True)
            self.step_btn.setEnabled(True)
        # Process events immediately to ensure button updates
        QApplication.processEvents()
            
    def _resume_loop(self):
        """Resume the simulation loop"""
        if self.remote_mode:
            try:
                response = self.remote_client._send_command('resume_thread', {})
                if response.get('status') == 'ok':
                    self.resume_btn.setEnabled(False)
                    self.pause_btn.setEnabled(True)
                    self.step_btn.setEnabled(False)
                else:
                    raise Exception(response.get('error', 'Unknown error'))
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to resume:\n{str(e)}")
        elif self.supervisor_thread:
            self.supervisor_thread.resume_loop()
            self.resume_btn.setEnabled(False)
            self.pause_btn.setEnabled(True)
            self.step_btn.setEnabled(False)
        # Process events immediately to ensure button updates
        QApplication.processEvents()
            
    def _step_loop(self):
        """Execute one iteration"""
        if self.remote_mode:
            try:
                response = self.remote_client._send_command('step_thread', {})
                if response.get('status') == 'ok':
                    self.status_bar.showMessage("Step requested")
                else:
                    raise Exception(response.get('error', 'Unknown error'))
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to step:\n{str(e)}")
        elif self.supervisor_thread:
            self.supervisor_thread.step_loop()
        # Process events immediately to ensure button updates
        QApplication.processEvents()
            
    def _stop_loop(self):
        """Stop the simulation loop"""
        if self.remote_mode:
            try:
                self.remote_client._send_command('stop_thread', {})
                # Don't wait for thread to finish on remote side
                # Just update UI immediately
                self.start_btn.setEnabled(True)
                self.pause_btn.setEnabled(False)
                self.resume_btn.setEnabled(False)
                self.stop_btn.setEnabled(False)
                self.step_btn.setEnabled(True)
                self.gui_update_timer.stop()
                self.status_bar.showMessage("Stop requested")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to stop:\n{str(e)}")
        elif self.supervisor_thread:
            self.supervisor_thread.stop_loop()
            self.supervisor_thread.wait()  # Wait for thread to finish
            # Stop GUI timer
            self.gui_update_timer.stop()
        # Process events immediately to ensure button updates
        QApplication.processEvents()
            
    def _reset_simulation(self):
        """Reset the simulation"""
        if not self.remote_mode and self.supervisor_thread and self.supervisor_thread.isRunning():
            QMessageBox.warning(self, "Warning", "Stop the loop before resetting!")
            return
            
        if self.supervisor:
            try:
                if self.remote_mode:
                    # Remote mode - call reset via client
                    response = self.remote_client._send_command('reset', {})
                    if response.get('status') != 'ok':
                        raise Exception(response.get('error', 'Unknown error'))
                else:
                    # Local mode
                    self.supervisor.reset()
                
                self.iter_label.setText("0")
                self.strehl_plot.clear()
                self.framerate_plot.clear()
                self.status_bar.showMessage("Simulation reset")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Reset failed:\n{str(e)}")
    
    def _toggle_loop_state(self, closed):
        """Toggle between open and closed loop."""
        if not self.supervisor and not self.remote_mode:
            return
        
        try:
            if self.remote_mode:
                cmd = 'close_loop' if closed else 'open_loop'
                response = self.remote_client._send_command(cmd, {})
                if response.get('status') != 'ok':
                    raise Exception(response.get('error', 'Unknown error'))
            else:
                # Handle two-stages mode
                if self.is_two_stages:
                    # Apply to both stages
                    if closed:
                        self.first_stage_supervisor.rtc.close_loop()
                        self.second_stage_supervisor.rtc.close_loop()
                    else:
                        self.first_stage_supervisor.rtc.open_loop()
                        self.second_stage_supervisor.rtc.open_loop()
                else:
                    if closed:
                        self.supervisor.rtc.close_loop()
                    else:
                        self.supervisor.rtc.open_loop()
            
            if closed:
                self.close_loop_btn.setText("Loop Closed")
                self.status_bar.showMessage("Loop closed")
            else:
                self.close_loop_btn.setText("Loop Opened")
                self.status_bar.showMessage("Loop opened")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to toggle loop:\n{str(e)}")
            # Revert button state
            self.close_loop_btn.setChecked(not closed)
    
    def _change_displayed_stage(self, index):
        """Change which stage is displayed in two-stages mode."""
        if not self.is_two_stages:
            return
        
        if self.remote_mode:
            # Remote mode: send command to server and update local display options
            try:
                # Send command to remote server to change displayed stage
                response = self.remote_client._send_command('set_displayed_stage', {'stage_index': index})
                if response.get('status') != 'ok':
                    raise Exception(response.get('error', 'Unknown error'))
                
                # Get stage config from server
                config_response = self.remote_client._send_command('get_stage_config', {'stage_index': index})
                if config_response.get('status') == 'ok':
                    stage_config_info = config_response.get('result')
                    
                    # Create mock config for display options
                    class MockStageConfig:
                        def __init__(self, config_info):
                            n_targets = config_info.get('n_targets', 0)
                            self.p_targets = [f"target_{i}" for i in range(n_targets)] if n_targets > 0 else None
                            
                            n_wfs = config_info.get('n_wfs', 0)
                            self.p_wfss = [f"wfs_{i}" for i in range(n_wfs)] if n_wfs > 0 else None
                            
                            n_dms = config_info.get('n_dms', 0)
                            self.p_dms = [f"dm_{i}" for i in range(n_dms)] if n_dms > 0 else None
                            
                            n_coronos = config_info.get('n_coronos', 0)
                            self.p_coronos = [f"corono_{i}" for i in range(n_coronos)] if n_coronos > 0 else None
                            
                            self.p_atmos = "mock"  # Always present for remote
                    
                    stage_config = MockStageConfig(stage_config_info)
                    self._populate_display_options(stage_config)
                
                stage_name = "Second Stage" if index == 0 else "First Stage"
                self.status_bar.showMessage(f"Now displaying {stage_name}")
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to change stage:\n{str(e)}")
                return
        else:
            # Local mode: update supervisor thread and display options
            if not self.supervisor_thread:
                return
            
            # Update which supervisor is used for telemetry
            self.supervisor_thread.set_displayed_stage(index)
            
            # Update display options to match the selected stage
            if index == 0:
                # Second stage
                stage_config = self.second_stage_supervisor.config
                stage_name = "Second Stage"
            else:
                # First stage
                stage_config = self.first_stage_supervisor.config
                stage_name = "First Stage"
            
            # Repopulate display options with the selected stage's configuration
            self._populate_display_options(stage_config)
            
            # Show status message
            self.status_bar.showMessage(f"Now displaying {stage_name}")
            
            # Force immediate telemetry update to show new stage
            if hasattr(self.supervisor_thread, '_emit_telemetry'):
                self.supervisor_thread._emit_telemetry()
    
    def _reset_strehl(self):
        """Reset Strehl ratio on all targets."""
        if not self.supervisor and not self.remote_mode:
            return
        
        try:
            if self.remote_mode:
                response = self.remote_client._send_command('reset_strehl', {})
                if response.get('status') != 'ok':
                    raise Exception(response.get('error', 'Unknown error'))
            else:
                if self.is_two_stages:
                    # Two-stages mode: reset both stages
                    self.two_stages_manager.reset_exposure()
                else:
                    # Reset Strehl for all targets
                    for tar_idx in range(len(self.config.p_targets)):
                        self.supervisor.target.reset_strehl(tar_idx)
            
            # Clear the Strehl plot
            self.strehl_plot.clear()
            
            self.status_bar.showMessage("Strehl ratio reset on all targets")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to reset Strehl:\n{str(e)}")
    
    def _toggle_atmosphere(self, enabled):
        """Enable or disable atmosphere."""
        if not self.supervisor and not self.remote_mode:
            return
        
        # Check if atmos exists (handle two-stages mode)
        has_atmos = False
        if self.is_two_stages:
            has_atmos = self.first_stage_supervisor.atmos is not None
        elif self.supervisor:
            has_atmos = self.supervisor.atmos is not None
        
        if not has_atmos and not self.remote_mode:
            return
        
        try:
            if self.remote_mode:
                response = self.remote_client._send_command('enable_atmos', {'enabled': enabled})
                if response.get('status') != 'ok':
                    raise Exception(response.get('error', 'Unknown error'))
            else:
                if self.is_two_stages:
                    # Two-stages: only enable/disable first stage atmos
                    # Second stage atmos is always disabled by TwoStagesManager
                    self.first_stage_supervisor.atmos.enable_atmos(enabled)
                else:
                    self.supervisor.atmos.enable_atmos(enabled)
            
            if enabled:
                self.enable_atmos_btn.setText("Atmos Enabled")
                self.status_bar.showMessage("Atmosphere enabled")
            else:
                self.enable_atmos_btn.setText("Atmos Disabled")
                self.status_bar.showMessage("Atmosphere disabled")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to toggle atmosphere:\n{str(e)}")
            # Revert button state
            self.enable_atmos_btn.setChecked(not enabled)
                
    def _update_monitoring_freq(self, value):
        """Update the monitoring frequency"""
        if self.remote_mode:
            # In remote mode, monitoring frequency is controlled by server
            self.status_bar.showMessage("Note: Monitoring frequency controlled by server in remote mode")
        elif self.supervisor_thread:
            self.supervisor_thread.set_monitoring_frequency(value)
            
    # Slots for supervisor thread signals
    def _on_iteration_done(self, iter_num):
        """Update iteration counter"""
        self.iter_label.setText(str(iter_num))
        
    def _on_strehl_updated(self, se_strehl, le_strehl):
        """Update Strehl ratio plot"""
        iter_num = int(self.iter_label.text())
        self.strehl_plot.add_point(iter_num, se_strehl, le_strehl)
        self.status_bar.showMessage(f"SE Strehl: {se_strehl:.4f} | LE Strehl: {le_strehl:.4f}")
        
    def _on_framerate_updated(self, framerate):
        """Update framerate plot"""
        iter_num = int(self.iter_label.text())
        self.framerate_plot.add_point(iter_num, framerate)
    
    # New signal handlers for selective data display
    def _update_display_if_matching(self, data_type, data, index=None):
        """Helper method to update displays if they match the data type and index."""
        # Update main tabs
        for display in [self.atmos_display, self.target_display, self.wfs_display,
                       self.dm_display, self.corono_display]:
            if display.current_data_type == data_type:
                if index is None or display.current_index == index:
                    display.update_image(data)
        
        # Update custom layout displays
        for display in self.layout_manager.get_displays():
            if display.current_data_type == data_type:
                if index is None or display.current_index == index:
                    self.layout_manager.update_display(display, data)
    
    def _on_atmos_phase_updated(self, phase):
        """Update atmosphere display"""
        self._update_display_if_matching('atmos_phase', phase)
    
    def _on_target_psf_se_updated(self, psf_image, tar_index):
        """Update target PSF SE display"""
        self._update_display_if_matching('target_psf_se', psf_image, tar_index)
    
    def _on_target_psf_le_updated(self, psf_image, tar_index):
        """Update target PSF LE display"""
        self._update_display_if_matching('target_psf_le', psf_image, tar_index)
    
    def _on_target_phase_updated(self, phase, tar_index):
        """Update target phase display"""
        self._update_display_if_matching('target_phase', phase, tar_index)
        
    def _on_wfs_image_updated(self, wfs_image, wfs_index):
        """Update WFS image display"""
        self._update_display_if_matching('wfs_image', wfs_image, wfs_index)
    
    def _on_wfs_phase_updated(self, wfs_phase, wfs_index):
        """Update WFS phase display"""
        self._update_display_if_matching('wfs_phase', wfs_phase, wfs_index)
    
    def _on_dm_shape_updated(self, dm_shape, dm_index):
        """Update DM shape display"""
        self._update_display_if_matching('dm_shape', dm_shape, dm_index)
    
    def _on_corono_image_updated(self, corono_image, coro_index):
        """Update coronagraph image display"""
        self._update_display_if_matching('corono_image', corono_image, coro_index)
        
    def _on_error(self, error_msg):
        """Handle error from supervisor thread"""
        QMessageBox.critical(self, "Simulation Error", error_msg)
        self.status_bar.showMessage(f"Error: {error_msg}")
        
    def _on_loop_finished(self):
        """Handle loop finished"""
        self.start_btn.setEnabled(True)
        self.pause_btn.setEnabled(False)
        self.resume_btn.setEnabled(False)
        self.stop_btn.setEnabled(False)
        self.step_btn.setEnabled(True)
        self.load_btn.setEnabled(True)
        # Stop GUI timer
        self.gui_update_timer.stop()
        
    def _on_status_changed(self, status):
        """Update status label"""
        self.status_label.setText(status)
        
    def _process_gui_events(self):
        """Periodically process GUI events to ensure responsiveness during simulation."""
        QApplication.processEvents()
        
    def _show_about(self):
        """Show about dialog"""
        QMessageBox.about(
            self,
            "About COMPASS GUI",
            "<h3>COMPASS GUI</h3>"
            "<p>Graphical User Interface for COMPASS AO Simulation</p>"
            "<p>Copyright (C) 2011-2025 COSMIC Team</p>"
            "<p><a href='https://github.com/COSMIC-RTC/compass'>https://github.com/COSMIC-RTC/compass</a></p>"
        )
    
    def _toggle_remote_mode(self, state):
        """Toggle between local and remote supervisor mode"""
        self.remote_mode = (state == Qt.CheckState.Checked.value)
        
        if self.remote_mode:
            # Entering remote mode
            self.connect_remote_btn.setVisible(True)
            self.connect_remote_btn.setEnabled(True)
            self.load_btn.setEnabled(False)  # Can't load local params in remote mode
            self.init_btn.setEnabled(False)
            self.param_label.setText("Remote mode - connect to server")
            self.status_label.setText("Remote mode (not connected)")
        else:
            # Entering local mode
            self.connect_remote_btn.setVisible(False)
            self.load_btn.setEnabled(True)
            if self.remote_client:
                self.remote_client.disconnect()
                self.remote_client = None
            self.param_label.setText("No parameter file loaded" if not self.param_file else self.param_file)
            self.status_label.setText("Local mode")
            self.init_btn.setEnabled(self.param_file is not None)
    
    def _connect_remote(self):
        """Open connection dialog and connect to remote supervisor server."""
        dialog = ConnectionDialog(self)
        
        if dialog.exec() == QDialog.DialogCode.Accepted:
            conn_info = dialog.get_connection_info()
            self.connection_info = conn_info
            
            # Create remote client
            self.remote_client = RemoteSupervisorClient(
                server_address=conn_info['server_address'],
                command_port=conn_info['command_port'],
                telemetry_port=conn_info['telemetry_port']
            )
            
            # Try to connect
            self.status_label.setText("Connecting...")
            QApplication.processEvents()
            
            if self.remote_client.connect():
                # Connection successful
                self.param_label.setText(f"Connected to {conn_info['server_address']}")
                self.status_label.setText("Connected")
                self.connect_remote_btn.setText("Disconnect")
                self.connect_remote_btn.clicked.disconnect()
                self.connect_remote_btn.clicked.connect(self._disconnect_remote)
                
                # Use remote client as supervisor
                self.supervisor = self.remote_client
                self.config = self.remote_client.config
                
                # Check if remote server is in two-stages mode
                config_info = self.remote_client.config_info
                self.is_two_stages = config_info.get('is_two_stages', False)
                
                # Enable stage selector if two-stages mode
                if self.is_two_stages:
                    self.stage_selector.setEnabled(True)
                    self.stage_selector.setVisible(True)
                    self.param_label.setText(f"Connected to {conn_info['server_address']} (Two-Stages)")
                
                # Populate display options based on remote config
                self._populate_display_options()
                
                # Enable controls
                self._enable_simulation_controls()
                
                # Setup telemetry receiver
                self._setup_remote_telemetry()
                
                # Update IPython console namespace if available
                if IPYTHON_AVAILABLE and self.kernel_manager:
                    self._update_console_namespace()
                
                QMessageBox.information(
                    self,
                    "Connected",
                    f"Successfully connected to remote supervisor at {conn_info['server_address']}"
                )
            else:
                # Connection failed
                QMessageBox.critical(
                    self,
                    "Connection Failed",
                    f"Failed to connect to {conn_info['server_address']}:{conn_info['command_port']}\n\n"
                    "Make sure the remote server is running."
                )
                self.status_label.setText("Connection failed")
                self.remote_client = None
    
    def _disconnect_remote(self):
        """Disconnect from remote server"""
        if self.remote_client:
            self.remote_client.disconnect()
            self.remote_client = None
            self.supervisor = None
            self.config = None
            
            self.param_label.setText("Remote mode - not connected")
            self.status_label.setText("Disconnected")
            self.connect_remote_btn.setText("Connect...")
            self.connect_remote_btn.clicked.disconnect()
            self.connect_remote_btn.clicked.connect(self._connect_remote)
            
            # Disable controls
            self._disable_simulation_controls()
            
            QMessageBox.information(self, "Disconnected", "Disconnected from remote supervisor")
    
    def _setup_remote_telemetry(self):
        """Setup telemetry callbacks for remote supervisor"""
        if not self.remote_client:
            return
        
        # Register callbacks for telemetry topics
        self.remote_client.register_telemetry_callback('image', self._on_remote_image)
        self.remote_client.register_telemetry_callback('status', self._on_remote_status)
        
        # Start a timer to poll telemetry queue
        self.remote_telemetry_timer = QTimer()
        self.remote_telemetry_timer.timeout.connect(self._poll_remote_telemetry)
        self.remote_telemetry_timer.setInterval(10)  # Poll every 10ms
        self.remote_telemetry_timer.start()
    
    def _poll_remote_telemetry(self):
        """Poll telemetry queue from remote client"""
        if not self.remote_client:
            if hasattr(self, 'remote_telemetry_timer'):
                self.remote_telemetry_timer.stop()
            return
        
        # Process up to 10 telemetry messages per poll to avoid blocking
        for _ in range(10):
            telemetry = self.remote_client.get_telemetry(block=False)
            if telemetry is None:
                break
            
            topic, message = telemetry
            data = message.get('data', {})
            
            # Handle different telemetry types
            if topic == 'image':
                self._handle_remote_image(data)
            elif topic == 'status':
                self._handle_remote_status(data)
    
    def _on_remote_image(self, message):
        """Callback for remote image telemetry (runs in telemetry thread)"""
        # Just pass - actual processing happens in poll
        pass
    
    def _on_remote_status(self, message):
        """Callback for remote status telemetry (runs in telemetry thread)"""
        # Just pass - actual processing happens in poll
        pass
    
    def _handle_remote_image(self, data):
        """Handle remote image data in GUI thread"""
        category = data.get('category')
        title = data.get('title')
        image = data.get('image')
        
        if image is None:
            return
        
        # Parse index from title if present (e.g., "Target 0 - PSF SE" -> index=0)
        import re
        index_match = re.search(r'\b(\d+)\b', title)
        data_index = int(index_match.group(1)) if index_match else 0
        
        # Update appropriate display based on category and title
        # Only update if this is the currently selected data type and index
        if category == 'atmos':
            current_type, _ = self.atmos_display.get_current_selection()
            if current_type == 'atmos_phase':
                self.atmos_display.update_image(image)
        
        elif category == 'target':
            current_type, current_idx = self.target_display.get_current_selection()
            # Check if this matches the currently selected target display
            if 'PSF SE' in title and current_type == 'target_psf_se' and current_idx == data_index:
                self.target_display.update_image(image)
            elif 'PSF LE' in title and current_type == 'target_psf_le' and current_idx == data_index:
                self.target_display.update_image(image)
            elif 'Phase' in title and current_type == 'target_phase' and current_idx == data_index:
                self.target_display.update_image(image)
        
        elif category == 'wfs':
            current_type, current_idx = self.wfs_display.get_current_selection()
            # Check if this matches the currently selected WFS display
            if 'Image' in title and current_type == 'wfs_image' and current_idx == data_index:
                self.wfs_display.update_image(image)
            elif 'Phase' in title and current_type == 'wfs_phase' and current_idx == data_index:
                self.wfs_display.update_image(image)
        
        elif category == 'dm':
            current_type, current_idx = self.dm_display.get_current_selection()
            # Check if this matches the currently selected DM
            if current_type == 'dm_shape' and current_idx == data_index:
                self.dm_display.update_image(image)
        
        elif category == 'corono':
            current_type, current_idx = self.corono_display.get_current_selection()
            # Check if this matches the currently selected coronagraph
            if current_type == 'corono_image' and current_idx == data_index:
                self.corono_display.update_image(image)
        
        # Update custom layout displays
        for display in self.layout_manager.get_displays():
            current_type, current_idx = display.get_current_selection()
            # Match based on category and data type
            should_update = False
            
            if category == 'atmos' and current_type == 'atmos_phase':
                should_update = True
            elif category == 'target':
                if ('PSF SE' in title and current_type == 'target_psf_se' and current_idx == data_index):
                    should_update = True
                elif ('PSF LE' in title and current_type == 'target_psf_le' and current_idx == data_index):
                    should_update = True
                elif ('Phase' in title and current_type == 'target_phase' and current_idx == data_index):
                    should_update = True
            elif category == 'wfs':
                if ('Image' in title and current_type == 'wfs_image' and current_idx == data_index):
                    should_update = True
                elif ('Phase' in title and current_type == 'wfs_phase' and current_idx == data_index):
                    should_update = True
            elif category == 'dm' and current_type == 'dm_shape' and current_idx == data_index:
                should_update = True
            elif category == 'corono' and current_type == 'corono_image' and current_idx == data_index:
                should_update = True
            
            if should_update:
                self.layout_manager.update_display(display, image)
    
    def _handle_remote_status(self, data):
        """Handle remote status data in GUI thread"""
        if 'iter' in data:
            self.iter_label.setText(str(data['iter']))
        if 'strehl_se' in data and 'strehl_le' in data:
            self._on_strehl_updated(data['strehl_se'], data['strehl_le'])
        if 'framerate' in data:
            self._on_framerate_updated(data['framerate'])
    
    def _enable_simulation_controls(self):
        """Enable simulation control buttons"""
        self.start_btn.setEnabled(True)
        self.step_btn.setEnabled(True)
        self.reset_btn.setEnabled(True)
        self.close_loop_btn.setEnabled(True)
        self.reset_strehl_btn.setEnabled(True)
        self.enable_atmos_btn.setEnabled(True)
    
    def _disable_simulation_controls(self):
        """Disable simulation control buttons"""
        self.start_btn.setEnabled(False)
        self.pause_btn.setEnabled(False)
        self.resume_btn.setEnabled(False)
        self.step_btn.setEnabled(False)
        self.stop_btn.setEnabled(False)
        self.reset_btn.setEnabled(False)
        self.close_loop_btn.setEnabled(False)
        self.reset_strehl_btn.setEnabled(False)
        self.enable_atmos_btn.setEnabled(False)
        
    def closeEvent(self, event):
        """Handle window close event"""
        # Stop GUI timer
        self.gui_update_timer.stop()
        
        # Stop remote telemetry timer if exists
        if hasattr(self, 'remote_telemetry_timer'):
            self.remote_telemetry_timer.stop()
        
        # Disconnect remote client if connected
        if self.remote_client:
            self.remote_client.disconnect()
        
        if self.supervisor_thread and self.supervisor_thread.isRunning():
            reply = QMessageBox.question(
                self,
                "Confirm Exit",
                "Simulation is still running. Do you want to stop it and exit?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
            )
            
            if reply == QMessageBox.StandardButton.Yes:
                self.supervisor_thread.stop_loop()
                self.supervisor_thread.wait()
                
                # Cleanup IPython console
                if IPYTHON_AVAILABLE and self.kernel_client:
                    self.kernel_client.stop_channels()
                if IPYTHON_AVAILABLE and self.kernel_manager:
                    self.kernel_manager.shutdown_kernel()
                
                event.accept()
            else:
                event.ignore()
        else:
            # Cleanup IPython console
            if IPYTHON_AVAILABLE and self.kernel_client:
                self.kernel_client.stop_channels()
            if IPYTHON_AVAILABLE and self.kernel_manager:
                self.kernel_manager.shutdown_kernel()
            
            event.accept()
