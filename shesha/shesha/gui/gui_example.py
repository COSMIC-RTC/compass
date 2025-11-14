#!/usr/bin/env python

"""
Example script showing how to use COMPASS GUI programmatically
and how to extend it with custom functionality
"""

from PyQt6.QtWidgets import QApplication, QPushButton, QVBoxLayout, QWidget
from shesha.gui.main_window import CompassMainWindow
from shesha.config import ParamConfig
import sys


def create_custom_gui():
    """Example: Create a custom GUI with additional controls"""
    
    app = QApplication(sys.argv)
    
    # Create main window
    window = CompassMainWindow()
    
    # Example: Add custom button to control panel
    # This shows how you can extend the GUI
    def custom_action():
        """Custom action to perform on the supervisor"""
        if window.supervisor:
            print("Executing custom action on supervisor...")
            # Example: Send custom command to supervisor thread
            window.supervisor_thread.send_custom_command(
                'reset',  # Any supervisor method
            )
    
    # You can add custom widgets programmatically
    # custom_btn = QPushButton("Custom Action")
    # custom_btn.clicked.connect(custom_action)
    # window.centralWidget().layout().addWidget(custom_btn)
    
    window.show()
    return app.exec()


def load_and_initialize(param_file):
    """Example: Automatically load parameters and initialize"""
    
    app = QApplication(sys.argv)
    window = CompassMainWindow()
    
    try:
        # Load parameters
        window.config = ParamConfig(param_file)
        window.param_file = param_file
        window.param_label.setText(param_file.split('/')[-1])
        window.init_btn.setEnabled(True)
        
        # Auto-initialize (optional)
        # window._initialize_supervisor()
        
        print(f"Loaded: {param_file}")
        print("GUI ready. Click 'Initialize Supervisor' to begin.")
        
    except Exception as e:
        print(f"Error loading parameters: {e}")
    
    window.show()
    return app.exec()


def monitor_simulation(param_file):
    """
    Example: Set up GUI with custom monitoring
    This shows how to connect to supervisor signals
    """
    
    app = QApplication(sys.argv)
    window = CompassMainWindow()
    
    # Custom handler for Strehl updates
    def on_strehl_update(se, le):
        print(f"Custom monitoring - SE: {se:.4f}, LE: {le:.4f}")
        # You could log to file, send to database, etc.
    
    # Load parameters
    window.config = ParamConfig(param_file)
    window.param_file = param_file
    window.param_label.setText(param_file.split('/')[-1])
    window.init_btn.setEnabled(True)
    
    # After initialization, you can connect to signals
    # Note: This should be done after supervisor thread is created
    original_init = window._initialize_supervisor
    
    def custom_init():
        original_init()
        # Now connect our custom handler
        if window.supervisor_thread:
            window.supervisor_thread.strehl_updated.connect(on_strehl_update)
            print("Custom monitoring enabled")
    
    window._initialize_supervisor = custom_init
    
    window.show()
    return app.exec()


if __name__ == "__main__":
    # Simple usage examples:
    
    # 1. Basic GUI
    # create_custom_gui()
    
    # 2. Load specific parameter file
    # if len(sys.argv) > 1:
    #     load_and_initialize(sys.argv[1])
    # else:
    #     print("Usage: python gui_example.py <parameter_file>")
    
    # 3. With custom monitoring
    if len(sys.argv) > 1:
        monitor_simulation(sys.argv[1])
    else:
        print("Usage: python gui_example.py <parameter_file>")
        print("\nOr edit this file to try different examples.")
