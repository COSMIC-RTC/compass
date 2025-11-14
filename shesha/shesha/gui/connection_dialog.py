"""
Connection dialog for configuring remote supervisor connection

Allows users to input server address and ports for remote operation.
"""

from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QSpinBox, QGroupBox, QFormLayout, QMessageBox
)
from PyQt6.QtCore import Qt


class ConnectionDialog(QDialog):
    """Dialog for configuring remote supervisor connection"""
    
    def __init__(self, parent=None, default_address="localhost", 
                 default_cmd_port=5555, default_tel_port=5556):
        super().__init__(parent)
        
        self.setWindowTitle("Remote Supervisor Connection")
        self.setModal(True)
        self.setMinimumWidth(400)
        
        # Connection settings
        self.server_address = default_address
        self.command_port = default_cmd_port
        self.telemetry_port = default_tel_port
        
        self._init_ui()
        
    def _init_ui(self):
        """Initialize the user interface"""
        layout = QVBoxLayout()
        
        # Server settings group
        server_group = QGroupBox("Server Settings")
        server_layout = QFormLayout()
        
        # Server address
        self.address_edit = QLineEdit(self.server_address)
        self.address_edit.setPlaceholderText("e.g., localhost or 192.168.1.100")
        server_layout.addRow("Server Address:", self.address_edit)
        
        # Command port
        self.cmd_port_spin = QSpinBox()
        self.cmd_port_spin.setRange(1024, 65535)
        self.cmd_port_spin.setValue(self.command_port)
        server_layout.addRow("Command Port:", self.cmd_port_spin)
        
        # Telemetry port
        self.tel_port_spin = QSpinBox()
        self.tel_port_spin.setRange(1024, 65535)
        self.tel_port_spin.setValue(self.telemetry_port)
        server_layout.addRow("Telemetry Port:", self.tel_port_spin)
        
        server_group.setLayout(server_layout)
        layout.addWidget(server_group)
        
        # Info label
        info_label = QLabel(
            "Connect to a remote COMPASS supervisor server.\n"
            "Make sure the server is running before connecting."
        )
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(info_label)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.connect_btn = QPushButton("Connect")
        self.connect_btn.clicked.connect(self._on_connect)
        self.connect_btn.setDefault(True)
        button_layout.addWidget(self.connect_btn)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(self.cancel_btn)
        
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
    
    def _on_connect(self):
        """Validate and accept connection settings"""
        address = self.address_edit.text().strip()
        
        if not address:
            QMessageBox.warning(
                self,
                "Invalid Address",
                "Please enter a server address."
            )
            return
        
        cmd_port = self.cmd_port_spin.value()
        tel_port = self.tel_port_spin.value()
        
        if cmd_port == tel_port:
            QMessageBox.warning(
                self,
                "Port Conflict",
                "Command and telemetry ports must be different."
            )
            return
        
        # Store values
        self.server_address = address
        self.command_port = cmd_port
        self.telemetry_port = tel_port
        
        self.accept()
    
    def get_connection_info(self):
        """
        Get connection information as a dictionary.
        
        Returns:
            dict with keys: server_address, command_port, telemetry_port
        """
        return {
            'server_address': self.server_address,
            'command_port': self.command_port,
            'telemetry_port': self.telemetry_port
        }
