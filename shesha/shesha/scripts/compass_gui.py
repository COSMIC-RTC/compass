#!/usr/bin/env python

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
COMPASS GUI launcher

Launch the COMPASS graphical user interface for adaptive optics simulation.

Usage:
  compass_gui.py [options] [<parameters_filename>] [<parameters_filename2>] [<frequency_ratio>]

with:
  - 'parameters_filename' the optional path to the parameters file to load on startup
  - 'parameters_filename2' the optional second parameters file for two-stages AO system
  - 'frequency_ratio' the frequency ratio between first and second stage (integer, default=1)

Options:
  -h --help          Show this help message and exit
"""

import sys
from PyQt6.QtWidgets import QApplication
from shesha.gui.main_window import CompassMainWindow
from shesha.config import ParamConfig
from docopt import docopt


def main():
    """Main entry point for COMPASS GUI"""
    arguments = docopt(__doc__)
    
    # Create Qt application
    app = QApplication(sys.argv)
    app.setApplicationName("COMPASS GUI")
    app.setOrganizationName("COSMIC Team")
    
    # Create main window
    window = CompassMainWindow()
    
    # Load parameter file(s) if provided
    param_file = arguments.get("<parameters_filename>")
    param_file2 = arguments.get("<parameters_filename2>")
    freq_ratio = arguments.get("<frequency_ratio>")
    
    # Check for two-stages mode
    is_two_stages = param_file2 is not None
    
    if param_file:
        try:
            window.config = ParamConfig(param_file)
            window.param_file = param_file
            
            if is_two_stages:
                # Two-stages mode
                window.param_file2 = param_file2
                window.frequency_ratio = int(freq_ratio) if freq_ratio else 1
                
                display_name = f"{param_file.split('/')[-1]} + {param_file2.split('/')[-1]} (2-stages)"
                window.param_label.setText(display_name)
                window.status_bar.showMessage(
                    f"Loaded two-stages parameters: {param_file} + {param_file2} (ratio={window.frequency_ratio})"
                )
            else:
                # Standard single-stage mode
                window.param_label.setText(param_file.split('/')[-1])
                window.status_bar.showMessage(f"Loaded parameters: {param_file}")
            
            window.init_btn.setEnabled(True)
        except Exception as e:
            print(f"Warning: Failed to load parameter file: {e}")
    
    # Show window and run application
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()