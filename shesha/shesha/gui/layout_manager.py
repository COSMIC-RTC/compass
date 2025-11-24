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
Custom layout manager using pyqtgraph DockArea for flexible display arrangement
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QDialog,
    QListWidget, QDialogButtonBox, QLabel,
    QToolBar, QComboBox, QFileDialog, QMessageBox
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QAction
from pyqtgraph.dockarea import DockArea, Dock
from .display_widgets import SelectableImageDisplayWidget
import json


class LayoutManager(QWidget):
    """
    Flexible layout manager using pyqtgraph DockArea.
    Provides drag-and-drop, resizable docks, and floating windows.
    """
    
    # Signal emitted when a display selection changes (data_type, index, display_widget)
    display_selection_changed = pyqtSignal(str, int, object)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.display_options = {}  # Available options per category
        self.docks = {}  # Dict mapping dock names to (dock, display_widget)
        self.dock_counter = 0
        self._init_ui()
        
    def _init_ui(self):
        """Initialize the UI"""
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Toolbar with controls
        toolbar = QToolBar()
        toolbar.setMovable(False)
        
        add_action = QAction("➕ Add Display", self)
        add_action.triggered.connect(self._show_add_display_dialog)
        toolbar.addAction(add_action)
        
        remove_action = QAction("➖ Remove Display", self)
        remove_action.triggered.connect(self._show_remove_display_dialog)
        toolbar.addAction(remove_action)
        
        toolbar.addSeparator()
        
        # Quick add preset layouts
        preset_action = QAction("📐 Add Preset", self)
        preset_action.triggered.connect(self._show_preset_dialog)
        toolbar.addAction(preset_action)
        
        clear_action = QAction("🗑 Clear All", self)
        clear_action.triggered.connect(self._clear_all_displays)
        toolbar.addAction(clear_action)
        
        toolbar.addSeparator()
        
        # Save/Load layout
        save_action = QAction("💾 Save Layout", self)
        save_action.triggered.connect(self._save_layout)
        toolbar.addAction(save_action)
        
        load_action = QAction("📁 Load Layout", self)
        load_action.triggered.connect(self._load_layout)
        toolbar.addAction(load_action)
        
        toolbar.addSeparator()
        
        # Help text
        help_label = QLabel("  💡 Drag dock titles to rearrange | Double-click title to float")
        help_label.setStyleSheet("color: gray;")
        toolbar.addWidget(help_label)
        
        layout.addWidget(toolbar)
        
        # DockArea for flexible layout
        self.dock_area = DockArea()
        layout.addWidget(self.dock_area)
        
        self.setLayout(layout)
        
    def set_available_options(self, options_dict):
        """
        Set available display options.
        
        Args:
            options_dict: Dict mapping category names to list of (display_name, data_type, index)
        """
        self.display_options = options_dict
        
    def _show_add_display_dialog(self):
        """Show dialog to add a new display"""
        dialog = AddDisplayDialog(self.display_options, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            category, option_idx, position = dialog.get_selection()
            if category and option_idx >= 0:
                self._add_display(category, option_idx, position)
                
    def _add_display(self, category, option_idx, position='bottom', relative_to=None):
        """
        Add a new display widget in a dock.
        
        Args:
            category: Category name (Atmosphere, Target, WFS, DM, Coronagraph)
            option_idx: Index in the category options list
            position: Where to place dock ('left', 'right', 'top', 'bottom', 'above', 'below')
            relative_to: Dock to position relative to (None for first dock)
        """
        if category not in self.display_options:
            return
            
        options = self.display_options[category]
        if option_idx >= len(options):
            return
            
        display_name, data_type, data_index = options[option_idx]
        
        # Create display widget
        display = SelectableImageDisplayWidget(category)
        display.set_options(options)
        display.selector.setCurrentIndex(option_idx)
        
        # Connect selection change signal
        display.selection_changed.connect(
            lambda dt, idx: self.display_selection_changed.emit(dt, idx, display)
        )
        
        # Create dock
        self.dock_counter += 1
        dock_name = f"{category}_{self.dock_counter}"
        dock = Dock(display_name, size=(400, 400), closable=True)
        dock.addWidget(display)
        
        # Connect dock closed signal to cleanup
        dock.sigClosed.connect(lambda: self._on_dock_closed(dock_name))
        
        # Store dock and display
        self.docks[dock_name] = (dock, display)
        
        # Add to dock area
        if not self.docks or relative_to is None:
            # First dock or no relative position
            self.dock_area.addDock(dock)
        else:
            # Position relative to another dock
            if relative_to in self.docks:
                ref_dock = self.docks[relative_to][0]
                self.dock_area.addDock(dock, position, ref_dock)
            else:
                # Fallback to last added dock
                last_dock = list(self.docks.values())[-2][0]  # -2 because we just added current
                self.dock_area.addDock(dock, position, last_dock)
        
        # Emit initial selection
        current_type, current_index = display.get_current_selection()
        if current_type:
            self.display_selection_changed.emit(current_type, current_index, display)
        
    def _on_dock_closed(self, dock_name):
        """Handle dock closure"""
        if dock_name in self.docks:
            dock, display = self.docks[dock_name]
            display.setParent(None)
            display.deleteLater()
            del self.docks[dock_name]
            
    def _show_remove_display_dialog(self):
        """Show dialog to remove a display"""
        if not self.docks:
            return
            
        dialog = RemoveDisplayDialog(self.docks, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            dock_names = dialog.get_selected_dock_names()
            for dock_name in dock_names:
                self._remove_display(dock_name)
                
    def _remove_display(self, dock_name):
        """Remove a display by dock name"""
        if dock_name in self.docks:
            dock, display = self.docks[dock_name]
            self.dock_area.removeDock(dock)
            display.setParent(None)
            display.deleteLater()
            del self.docks[dock_name]
    
    def _show_preset_dialog(self):
        """Show dialog to add preset layout"""
        dialog = PresetLayoutDialog(self.display_options, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            preset_name = dialog.get_selected_preset()
            if preset_name:
                self._add_preset_layout(preset_name)
    
    def _add_preset_layout(self, preset_name):
        """Add a preset layout"""
        if preset_name == "2x2_overview":
            self._add_2x2_overview()
        elif preset_name == "wfs_comparison":
            self._add_wfs_comparison()
        elif preset_name == "target_analysis":
            self._add_target_analysis()
        elif preset_name == "ao_pipeline":
            self._add_ao_pipeline()
            
    def _add_2x2_overview(self):
        """Add 2x2 overview layout"""
        # Clear existing
        self._clear_all_displays()
        
        # Add 4 displays if available
        positions = [
            ('Target', 0, None, None),  # First
            ('WFS', 0, 'right', list(self.docks.keys())[0] if self.docks else None),
            ('DM', 0, 'bottom', list(self.docks.keys())[0] if len(self.docks) > 0 else None),
            ('Atmosphere', 0, 'right', list(self.docks.keys())[2] if len(self.docks) > 2 else None),
        ]
        
        for category, idx, pos, rel in positions:
            if category in self.display_options and len(self.display_options[category]) > idx:
                if rel and rel < len(self.docks):
                    rel_name = list(self.docks.keys())[rel] if isinstance(rel, int) else rel
                    self._add_display(category, idx, pos if pos else 'bottom', rel_name)
                else:
                    self._add_display(category, idx, pos if pos else 'bottom', None)
    
    def _add_wfs_comparison(self):
        """Add WFS comparison layout"""
        self._clear_all_displays()
        
        # Add all available WFS
        if 'WFS' in self.display_options:
            wfs_opts = self.display_options['WFS']
            # Add images in top row
            for i in range(0, len(wfs_opts), 2):  # Every 2nd is image (0, 2, 4...)
                pos = 'right' if i > 0 else None
                rel = list(self.docks.keys())[-1] if self.docks else None
                self._add_display('WFS', i, pos, rel)
    
    def _add_target_analysis(self):
        """Add target analysis layout"""
        self._clear_all_displays()
        
        # Add PSF SE, PSF LE, Phase for first target
        if 'Target' in self.display_options:
            target_opts = self.display_options['Target']
            if len(target_opts) >= 3:
                self._add_display('Target', 0, None, None)  # PSF SE
                self._add_display('Target', 1, 'right', list(self.docks.keys())[0])  # PSF LE
                self._add_display('Target', 2, 'right', list(self.docks.keys())[1])  # Phase
    
    def _add_ao_pipeline(self):
        """Add AO pipeline layout"""
        self._clear_all_displays()
        
        # Atmos -> WFS -> DM -> Target
        order = [('Atmosphere', 0), ('WFS', 0), ('DM', 0), ('Target', 0)]
        
        for i, (category, idx) in enumerate(order):
            if category in self.display_options and len(self.display_options[category]) > idx:
                pos = 'right' if i > 0 else None
                rel = list(self.docks.keys())[-1] if self.docks else None
                self._add_display(category, idx, pos, rel)
            
    def _clear_all_displays(self):
        """Remove all displays"""
        # Make a copy of keys since we'll modify during iteration
        dock_names = list(self.docks.keys())
        for dock_name in dock_names:
            self._remove_display(dock_name)
            
    def get_displays(self):
        """Get list of all display widgets"""
        return [display for _, display in self.docks.values()]
    
    def update_display(self, display_widget, image_data):
        """Update a specific display widget with new image data"""
        # Check if display is in our docks
        for dock, display in self.docks.values():
            if display is display_widget:
                display.update_image(image_data)
                break
    
    def _save_layout(self):
        """Save current layout configuration to a JSON file"""
        if not self.docks:
            QMessageBox.information(self, "Save Layout", "No displays to save!")
            return
        
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save Layout Configuration",
            "",
            "JSON Files (*.json);;All Files (*)"
        )
        
        if not filename:
            return
        
        try:
            # Build layout configuration
            layout_config = {
                'version': '1.0',
                'docks': []
            }
            
            # Save dock configuration
            for dock_name, (dock, display) in self.docks.items():
                category = display.category
                current_idx = display.selector.currentIndex()
                
                if current_idx >= 0 and category in self.display_options:
                    options = self.display_options[category]
                    if current_idx < len(options):
                        display_name, data_type, data_index = options[current_idx]
                        
                        dock_config = {
                            'dock_name': dock_name,
                            'category': category,
                            'option_idx': current_idx,
                            'display_name': display_name,
                            'data_type': data_type,
                            'data_index': data_index,
                            'log_scale': display.log_checkbox.isChecked(),
                            'colormap': display.colormap_selector.currentText()
                        }
                        layout_config['docks'].append(dock_config)
            
            # Save DockArea state (positions, sizes, etc.)
            dock_state = self.dock_area.saveState()
            layout_config['dock_area_state'] = dock_state
            
            # Write to file
            with open(filename, 'w') as f:
                json.dump(layout_config, f, indent=2)
            
            QMessageBox.information(
                self,
                "Layout Saved",
                f"Layout configuration saved to:\n{filename}"
            )
            
        except Exception as e:
            QMessageBox.critical(
                self,
                "Save Error",
                f"Failed to save layout:\n{str(e)}"
            )
    
    def _load_layout(self):
        """Load layout configuration from a JSON file"""
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Load Layout Configuration",
            "",
            "JSON Files (*.json);;All Files (*)"
        )
        
        if not filename:
            return
        
        try:
            # Read configuration file
            with open(filename, 'r') as f:
                layout_config = json.load(f)
            
            # Validate version
            if layout_config.get('version') != '1.0':
                QMessageBox.warning(
                    self,
                    "Version Mismatch",
                    "Layout file version not supported. Attempting to load anyway..."
                )
            
            # Clear existing layout
            self._clear_all_displays()
            
            # Restore docks
            docks_config = layout_config.get('docks', [])
            
            for dock_config in docks_config:
                category = dock_config.get('category')
                option_idx = dock_config.get('option_idx')
                
                if category and option_idx is not None:
                    # Verify option still exists
                    if category in self.display_options:
                        options = self.display_options[category]
                        if option_idx < len(options):
                            # Add display
                            self._add_display(category, option_idx)
                            
                            # Restore display settings
                            if self.docks:
                                last_dock_name = list(self.docks.keys())[-1]
                                _, display = self.docks[last_dock_name]
                                
                                # Restore log scale
                                if 'log_scale' in dock_config:
                                    display.log_checkbox.setChecked(dock_config['log_scale'])
                                
                                # Restore colormap
                                if 'colormap' in dock_config:
                                    colormap_name = dock_config['colormap']
                                    idx = display.colormap_selector.findText(colormap_name)
                                    if idx >= 0:
                                        display.colormap_selector.setCurrentIndex(idx)
            
            # Restore DockArea state (positions, sizes)
            if 'dock_area_state' in layout_config:
                try:
                    self.dock_area.restoreState(layout_config['dock_area_state'])
                except Exception as e:
                    # State restoration can fail if displays don't match exactly
                    print(f"Warning: Could not fully restore dock positions: {e}")
            
            QMessageBox.information(
                self,
                "Layout Loaded",
                f"Layout configuration loaded from:\n{filename}\n\n"
                f"Restored {len(docks_config)} displays."
            )
            
        except FileNotFoundError:
            QMessageBox.critical(
                self,
                "Load Error",
                f"File not found:\n{filename}"
            )
        except json.JSONDecodeError as e:
            QMessageBox.critical(
                self,
                "Load Error",
                f"Invalid JSON file:\n{str(e)}"
            )
        except Exception as e:
            QMessageBox.critical(
                self,
                "Load Error",
                f"Failed to load layout:\n{str(e)}"
            )


class AddDisplayDialog(QDialog):
    """Dialog to select which display to add"""
    
    def __init__(self, display_options, parent=None):
        super().__init__(parent)
        self.display_options = display_options
        self.selected_category = None
        self.selected_option_idx = -1
        self.selected_position = 'bottom'
        self._init_ui()
        
    def _init_ui(self):
        self.setWindowTitle("Add Display")
        self.setMinimumWidth(400)
        
        layout = QVBoxLayout()
        
        # Category selection
        layout.addWidget(QLabel("Select Category:"))
        self.category_list = QListWidget()
        self.category_list.addItems(sorted(self.display_options.keys()))
        self.category_list.currentItemChanged.connect(self._on_category_changed)
        layout.addWidget(self.category_list)
        
        # Display type selection
        layout.addWidget(QLabel("Select Display:"))
        self.option_list = QListWidget()
        layout.addWidget(self.option_list)
        
        # Position selection
        layout.addWidget(QLabel("Position:"))
        self.position_combo = QComboBox()
        self.position_combo.addItems(['bottom', 'right', 'top', 'left', 'above', 'below'])
        layout.addWidget(self.position_combo)
        
        # Buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        
        self.setLayout(layout)
        
        # Select first category by default
        if self.category_list.count() > 0:
            self.category_list.setCurrentRow(0)
            
    def _on_category_changed(self, current, previous):
        """Update options when category changes"""
        self.option_list.clear()
        if current:
            category = current.text()
            if category in self.display_options:
                options = self.display_options[category]
                for display_name, _, _ in options:
                    self.option_list.addItem(display_name)
                    
    def get_selection(self):
        """Get the selected category, option index, and position"""
        category_item = self.category_list.currentItem()
        if category_item:
            self.selected_category = category_item.text()
            self.selected_option_idx = self.option_list.currentRow()
            self.selected_position = self.position_combo.currentText()
            
        return self.selected_category, self.selected_option_idx, self.selected_position


class RemoveDisplayDialog(QDialog):
    """Dialog to select which displays to remove"""
    
    def __init__(self, docks, parent=None):
        super().__init__(parent)
        self.docks = docks
        self._init_ui()
        
    def _init_ui(self):
        self.setWindowTitle("Remove Displays")
        self.setMinimumWidth(300)
        
        layout = QVBoxLayout()
        
        layout.addWidget(QLabel("Select displays to remove:"))
        
        self.display_list = QListWidget()
        self.display_list.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
        
        # Add dock names to list
        for dock_name, (dock, _) in self.docks.items():
            self.display_list.addItem(f"{dock.label()} ({dock_name})")
            
        layout.addWidget(self.display_list)
        
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        
        self.setLayout(layout)
        
    def get_selected_dock_names(self):
        """Get list of selected dock names"""
        selected_items = self.display_list.selectedItems()
        dock_names = []
        for item in selected_items:
            # Extract dock_name from "Label (dock_name)" format
            text = item.text()
            if '(' in text and ')' in text:
                dock_name = text.split('(')[1].split(')')[0]
                dock_names.append(dock_name)
        return dock_names


class PresetLayoutDialog(QDialog):
    """Dialog to select a preset layout"""
    
    def __init__(self, display_options, parent=None):
        super().__init__(parent)
        self.display_options = display_options
        self.presets = {
            "2x2_overview": "2×2 Overview (Target, WFS, DM, Atmosphere)",
            "wfs_comparison": "WFS Comparison (All WFS side by side)",
            "target_analysis": "Target Analysis (PSF SE, PSF LE, Phase)",
            "ao_pipeline": "AO Pipeline (Atmos → WFS → DM → Target)"
        }
        self._init_ui()
        
    def _init_ui(self):
        self.setWindowTitle("Add Preset Layout")
        self.setMinimumWidth(400)
        
        layout = QVBoxLayout()
        
        layout.addWidget(QLabel("Select Preset:"))
        
        self.preset_list = QListWidget()
        for key, description in self.presets.items():
            self.preset_list.addItem(description)
            self.preset_list.item(self.preset_list.count() - 1).setData(Qt.ItemDataRole.UserRole, key)
            
        layout.addWidget(self.preset_list)
        
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        
        self.setLayout(layout)
        
        # Select first preset by default
        if self.preset_list.count() > 0:
            self.preset_list.setCurrentRow(0)
            
    def get_selected_preset(self):
        """Get the selected preset key"""
        current = self.preset_list.currentItem()
        if current:
            return current.data(Qt.ItemDataRole.UserRole)
        return None
