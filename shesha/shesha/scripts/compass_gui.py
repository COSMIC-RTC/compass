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
# Copyright (C) 2011-2024 COSMIC Team

"""
COMPASS GUI launcher

Launch the COMPASS graphical user interface for adaptive optics simulation.

Usage:
  compass_gui.py [options] [<parameters_filename>]

with 'parameters_filename' the optional path to the parameters file to load on startup

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
    
    # Load parameter file if provided
    param_file = arguments.get("<parameters_filename>")
    if param_file:
        try:
            window.config = ParamConfig(param_file)
            window.param_file = param_file
            window.param_label.setText(param_file.split('/')[-1])
            window.init_btn.setEnabled(True)
            window.status_bar.showMessage(f"Loaded parameters: {param_file}")
        except Exception as e:
            print(f"Warning: Failed to load parameter file: {e}")
    
    # Show window and run application
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()