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
Display widgets for real-time visualization of simulation data
"""

from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, QCheckBox
from PyQt6.QtCore import Qt, pyqtSignal
import numpy as np
import pyqtgraph as pg


class SelectableImageDisplayWidget(QWidget):
    """
    Widget for displaying 2D images with a dropdown selector.
    Allows user to choose what data to display.
    """
    
    # Signal emitted when selection changes (data_type, index)
    selection_changed = pyqtSignal(str, int)
    
    def __init__(self, title="Display", parent=None):
        super().__init__(parent)
        self.title = title
        self.category = title  # Store category for layout saving
        self.current_data_type = None
        self.current_index = 0
        self.log_scale = False
        self.current_image_data = None  # Store current image for re-processing
        self.current_colormap = 'viridis'  # Default colormap
        self._init_ui()
        
    def _init_ui(self):
        """Initialize the UI"""
        layout = QVBoxLayout()
        
        # Top bar with title and dropdown
        top_layout = QHBoxLayout()
        
        self.title_label = QLabel(self.title)
        self.title_label.setStyleSheet("font-weight: bold; font-size: 12pt;")
        top_layout.addWidget(self.title_label)
        
        top_layout.addStretch()
        
        self.selector = QComboBox()
        self.selector.currentIndexChanged.connect(self._on_selection_changed)
        top_layout.addWidget(self.selector)
        
        # Log scale checkbox
        self.log_checkbox = QCheckBox("Log Scale")
        self.log_checkbox.setChecked(False)
        self.log_checkbox.stateChanged.connect(self._on_log_scale_changed)
        top_layout.addWidget(self.log_checkbox)
        
        # Colormap selector
        top_layout.addWidget(QLabel("Color:"))
        self.colormap_selector = QComboBox()
        self.colormap_selector.addItems([
            'viridis', 'plasma', 'inferno', 'magma', 'cividis',
            'gray', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter',
            'jet', 'rainbow', 'turbo', 'hsv'
        ])
        self.colormap_selector.setCurrentText('viridis')
        self.colormap_selector.currentTextChanged.connect(self._on_colormap_changed)
        top_layout.addWidget(self.colormap_selector)
        
        layout.addLayout(top_layout)
        
        # Create pyqtgraph ImageView for fast image display
        self.image_view = pg.ImageView()
        self.image_view.ui.roiBtn.hide()  # Hide ROI button
        self.image_view.ui.menuBtn.hide()  # Hide menu button
        layout.addWidget(self.image_view)
        self._update_colormap()

        self.setLayout(layout)
        
    def set_options(self, options):
        """
        Set available display options.
        
        Args:
            options: List of tuples (display_name, data_type, index)
                    e.g., [("PSF SE", "target_psf_se", 0), ("PSF LE", "target_psf_le", 0)]
        """
        self.selector.blockSignals(True)
        self.selector.clear()
        
        self.options = options
        for display_name, data_type, index in options:
            self.selector.addItem(display_name)
        
        self.selector.blockSignals(False)
        
        if options:
            self._on_selection_changed(0)
    
    def _on_selection_changed(self, index):
        """Handle selection change"""
        if index < 0 or index >= len(self.options):
            return
            
        display_name, data_type, data_index = self.options[index]
        self.current_data_type = data_type
        self.current_index = data_index
        
        # Emit signal to request this data
        self.selection_changed.emit(data_type, data_index)
    
    def _on_log_scale_changed(self, state):
        """Handle log scale checkbox change"""
        self.log_scale = (state == Qt.CheckState.Checked.value)
        # Re-display current image with new scale
        if self.current_image_data is not None:
            self._display_image(self.current_image_data)
    
    def _on_colormap_changed(self, colormap_name):
        """Handle colormap selection change"""
        self.current_colormap = colormap_name
        self._update_colormap()
        # Re-display current image with new colormap
        if self.current_image_data is not None:
            self._display_image(self.current_image_data)
    
    def _update_colormap(self):
        """Update the ImageView colormap"""
        # Get the colormap from pyqtgraph or matplotlib
        try:
            # Try to get matplotlib colormap
            import matplotlib.cm as cm
            
            cmap = cm.get_cmap(self.current_colormap)
            # Convert to pyqtgraph format (Nx3 array of RGB values)
            colors = cmap(np.linspace(0, 1, 256))[:, :3] * 255
            colormap = pg.ColorMap(pos=np.linspace(0.0, 1.0, 256), color=colors)
            self.image_view.setColorMap(colormap)
        except Exception:
            # Fallback to pyqtgraph built-in colormaps
            pass
        
    def update_image(self, image_data):
        """Update the displayed image"""
        if image_data is not None:
            self.current_image_data = image_data.copy()
            self._display_image(image_data)
    
    def _display_image(self, image_data):
        """Display image with current scale setting"""
        if image_data is None:
            return
        
        if self.log_scale:
            # Apply log scale with safety for zero/negative values
            # Add small epsilon to avoid log(0)
            epsilon = 1e-10
            display_data = np.log10(np.abs(image_data) + epsilon)
        else:
            display_data = image_data
            
        self.image_view.setImage(display_data, autoRange=False, autoLevels=True, autoHistogramRange=True)
    
    def get_current_selection(self):
        """Get current data type and index"""
        return self.current_data_type, self.current_index


class ImageDisplayWidget(QWidget):
    """Widget for displaying 2D images (PSF, WFS, phase screens) - Simple version"""
    
    def __init__(self, title="Image", parent=None):
        super().__init__(parent)
        self.title = title
        self._init_ui()
        
    def _init_ui(self):
        """Initialize the UI"""
        layout = QVBoxLayout()
        
        # Title label
        self.title_label = QLabel(self.title)
        self.title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.title_label)
        
        # Create pyqtgraph ImageView for fast image display
        self.image_view = pg.ImageView()
        self.image_view.ui.roiBtn.hide()  # Hide ROI button
        self.image_view.ui.menuBtn.hide()  # Hide menu button
        layout.addWidget(self.image_view)
        
        self.setLayout(layout)
        
    def update_image(self, image_data):
        """Update the displayed image"""
        if image_data is not None:
            self.image_view.setImage(image_data, autoRange=False, autoLevels=True, autoHistogramRange=True)
            
    def set_title(self, title):
        """Update the title"""
        self.title = title
        self.title_label.setText(title)


class PlotWidget(QWidget):
    """Widget for plotting time-series data"""
    
    def __init__(self, title="Plot", ylabel="Value", max_points=1000, parent=None):
        super().__init__(parent)
        self.title = title
        self.ylabel = ylabel
        self.max_points = max_points
        self.data_x = []
        self.data_y = []
        self._init_ui()
        
    def _init_ui(self):
        """Initialize the UI"""
        layout = QVBoxLayout()
        
        # Create pyqtgraph PlotWidget
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')
        self.plot_widget.setTitle(self.title)
        self.plot_widget.setLabel('left', self.ylabel)
        self.plot_widget.setLabel('bottom', 'Iteration')
        self.plot_widget.showGrid(x=True, y=True, alpha=0.3)
        
        # Create plot curve
        self.curve = self.plot_widget.plot(pen=pg.mkPen(color='b', width=2))
        
        layout.addWidget(self.plot_widget)
        self.setLayout(layout)
        
    def add_point(self, x, y):
        """Add a data point to the plot"""
        self.data_x.append(x)
        self.data_y.append(y)
        
        # Limit the number of points to avoid memory issues
        if len(self.data_x) > self.max_points:
            self.data_x = self.data_x[-self.max_points:]
            self.data_y = self.data_y[-self.max_points:]
            
        self.curve.setData(self.data_x, self.data_y)
        
    def clear(self):
        """Clear the plot"""
        self.data_x = []
        self.data_y = []
        self.curve.setData([], [])


class DualPlotWidget(QWidget):
    """Widget for plotting two time-series on the same plot"""
    
    def __init__(self, title="Plot", ylabel="Value", legend1="Series 1", 
                 legend2="Series 2", max_points=1000, parent=None):
        super().__init__(parent)
        self.title = title
        self.ylabel = ylabel
        self.max_points = max_points
        self.data_x = []
        self.data_y1 = []
        self.data_y2 = []
        self.legend1 = legend1
        self.legend2 = legend2
        self._init_ui()
        
    def _init_ui(self):
        """Initialize the UI"""
        layout = QVBoxLayout()
        
        # Create pyqtgraph PlotWidget
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')
        self.plot_widget.setTitle(self.title)
        self.plot_widget.setLabel('left', self.ylabel)
        self.plot_widget.setLabel('bottom', 'Iteration')
        self.plot_widget.showGrid(x=True, y=True, alpha=0.3)
        self.plot_widget.addLegend()
        
        # Create plot curves
        self.curve1 = self.plot_widget.plot(pen=pg.mkPen(color='b', width=2), name=self.legend1)
        self.curve2 = self.plot_widget.plot(pen=pg.mkPen(color='r', width=2), name=self.legend2)
        
        layout.addWidget(self.plot_widget)
        self.setLayout(layout)
        
    def add_point(self, x, y1, y2):
        """Add data points to both plots"""
        self.data_x.append(x)
        self.data_y1.append(y1)
        self.data_y2.append(y2)
        
        # Limit the number of points
        if len(self.data_x) > self.max_points:
            self.data_x = self.data_x[-self.max_points:]
            self.data_y1 = self.data_y1[-self.max_points:]
            self.data_y2 = self.data_y2[-self.max_points:]
            
        self.curve1.setData(self.data_x, self.data_y1)
        self.curve2.setData(self.data_x, self.data_y2)
        
    def clear(self):
        """Clear both plots"""
        self.data_x = []
        self.data_y1 = []
        self.data_y2 = []
        self.curve1.setData([], [])
        self.curve2.setData([], [])
