# COMPASS GUI - Quick Start Guide

## Overview

The COMPASS GUI provides a graphical interface for running adaptive optics simulations with real-time visualization and control.

## Installation

### 1. Install GUI Dependencies

```bash
cd /home/micado/compass/shesha
pip install -r requirements-gui.txt
```

This installs:
- PyQt6 (GUI framework)
- pyqtgraph (fast plotting and image display)
- matplotlib (colormaps)
- qtconsole (IPython console widget)

### 2. Verify Installation

```bash
python -c "import PyQt6; import pyqtgraph; import qtconsole; print('GUI dependencies OK')"
```

## Running the GUI

### Basic Launch

```bash
cd /home/micado/compass/shesha
python shesha/scripts/compass_gui.py
```

### Launch with Parameter File

```bash
python shesha/scripts/compass_gui.py path/to/your/parameters.py
```

### From Python

```python
from shesha.gui import CompassMainWindow
from PyQt6.QtWidgets import QApplication
import sys

app = QApplication(sys.argv)
window = CompassMainWindow()
window.show()
sys.exit(app.exec())
```

## Using the GUI

### Step-by-Step Workflow

1. **Load Parameters**
   - Click "Load Parameters" button
   - Navigate to your parameter file (`.py`)
   - File loads and displays name in control panel

2. **Initialize Supervisor**
   - Click "Initialize Supervisor"
   - Wait for initialization (GPU setup, component creation)
   - Controls become enabled when ready

3. **Start Simulation**
   - Click "▶ Start Loop" to begin
   - Watch real-time displays update:
     - PSF images in "PSF" tab
     - WFS images in "WFS" tab
     - Phase screens in "Phase" tab
     - Strehl ratio plot (SE and LE)
     - Framerate plot

4. **Control the Loop**
   - **Pause** (⏸): Pause execution, GUI stays responsive
   - **Resume** (▶): Continue from pause
   - **Step** (⏭): Execute one iteration (when paused)
   - **Stop** (⏹): Stop the loop completely
   - **Reset** (↻): Return to initial conditions
   
5. **Advanced Controls**
   - **Close Loop** / **Open Loop**: Toggle RTC feedback control
   - **Reset Strehl**: Reset Strehl ratio measurements on all targets
   - **Enable Atmos** / **Disable Atmos**: Toggle atmospheric turbulence on/off

### Display Tabs

- **PSF Tab**: Target point spread function
- **WFS Tab**: Wavefront sensor images  
- **Phase Tab**: Atmospheric phase screens

### Performance Tuning

**GUI Update Frequency**
- Controls how often displays refresh (in frames)
- Default: 10 frames
- Lower = more responsive, higher CPU usage
- Higher = better loop performance, less responsive display
- Adjust in "GUI Update Freq (frames)" spinbox

### IPython Console

**Enable the Console:**
1. Check the "Show Console" checkbox at the bottom of the GUI
2. IPython console appears with full interactive capabilities

**Quick Examples:**
```python
# Get iteration count
supervisor.get_frame_counter()

# Inspect configuration
config.p_wfss[0].nxsub

# Get current RTC error
supervisor.rtc.get_err()

# Access atmospheric phase
phase = supervisor.atmos.get_atmos_layer(0)
```

**Available Objects:**
- `supervisor` - CompassSupervisor instance
- `config` - Configuration parameters
- `supervisor_thread` - Worker thread
- `gui` - Main window
- `np` - NumPy

**Features:**
- Tab completion
- Command history (↑↓)
- IPython magic commands (%timeit, %run, etc.)
- Syntax highlighting

## Architecture

```
┌─────────────────────────────────────────┐
│         Main Thread (GUI)               │
│  - User Interface                       │
│  - Display Updates                      │
│  - Event Handling                       │
│  - IPython Console                      │
└────────────┬────────────────────────────┘
             │ Signals/Slots
             │ (Thread-safe)
┌────────────▼────────────────────────────┐
│      Supervisor Thread                  │
│  - AO Loop Execution                    │
│  - Supervisor.next() calls              │
│  - Command Queue Processing             │
│  - Telemetry Emission                   │
└─────────────────────────────────────────┘
```

### Key Features

✅ **Non-blocking**: Loop runs independently from GUI
✅ **Thread-safe**: Signals/slots for communication
✅ **Real-time control**: Send commands while running
✅ **Real-time visualization**: Live PSF, WFS, phase
✅ **Performance monitoring**: Strehl, framerate plots
✅ **Interactive console**: IPython for advanced control
✅ **Remote mode**: Run supervisor on server, GUI on client

## Remote Mode

The GUI can run in **remote mode**, allowing the supervisor to run on one machine (e.g., a GPU server) while the GUI runs on another (e.g., your laptop).

### Architecture

```
┌─────────────────────────┐         Network          ┌──────────────────────────┐
│  Server (GPU Machine)   │                          │  Client (Your Machine)   │
│  ─────────────────      │                          │  ───────────────────     │
│                         │                          │                          │
│  Supervisor + AO Loop   │◄────── Commands ─────────┤  GUI + Controls          │
│  ZMQ Server             │        (REQ/REP)         │  ZMQ Client              │
│  - Command Socket       │                          │  - Display Panels        │
│  - Telemetry Socket     │──── Telemetry Stream ───►│  - IPython Console       │
│                         │        (PUB/SUB)         │  - Performance Plots     │
└─────────────────────────┘                          └──────────────────────────┘
```

### Running Remote Server

On the server machine (with GPU):

```bash
# Start remote server with a parameter file
python shesha/run_remote_server.py data/par/par4bench/scao_sh_16x16_8pix.py

# Custom ports
python shesha/run_remote_server.py params.py --cmd-port 6000 --tel-port 6001

# Auto-start the loop
python shesha/run_remote_server.py params.py --auto-start

# Bind to specific interface
python shesha/run_remote_server.py params.py --host 192.168.1.100
```

The server will print:
```
Server listening on:
  Command port:   5555
  Telemetry port: 5556
Server is ready for connections
```

### Connecting from GUI

On the client machine (no full COMPASS installation needed):

1. **Install minimal dependencies:**
   ```bash
   pip install -r shesha/requirements-remote-client.txt
   ```

2. **Launch the remote client:**
   ```bash
   # Auto-connect to server
   python shesha/run_remote_client.py --server 192.168.1.100 --auto-connect
   
   # Or launch GUI and connect manually
   python shesha/run_remote_client.py --server gpu-server
   ```

   The `run_remote_client.py` script:
   - ✅ Does NOT require full COMPASS installation
   - ✅ No CUDA, sutra, or carma modules needed
   - ✅ Only needs PyQt6, pyqtgraph, pyzmq, numpy
   - ✅ Automatically enables remote mode
   - ✅ Can auto-connect with `--auto-connect` flag

3. **Or use the full GUI** (if you have COMPASS installed):
   ```bash
   python -m shesha.gui.main
   ```
   Then manually:
   - Check "Remote Mode" checkbox
   - Click "Connect..." button
   - Enter server address and ports
   - Click "Connect"

### Network Requirements

- **Firewall**: Ensure ports 5555 and 5556 (or custom) are open
- **Bandwidth**: ~10-100 Mbps for smooth operation (depends on display rate and image sizes)
- **Latency**: < 100ms recommended for responsive controls
- **LAN recommended**: Works best on local network; WAN/Internet possible but slower

### Security Notes

⚠️ **Current implementation has no authentication or encryption**

For production use on untrusted networks:
- Use SSH tunneling: `ssh -L 5555:localhost:5555 -L 5556:localhost:5556 user@server`
- Run on VPN
- Add firewall rules to restrict access
- Consider adding authentication layer (future enhancement)

### Remote Mode Features

✅ All control buttons work (start, pause, step, stop)
✅ Loop control (open/close), atmosphere toggle, Strehl reset
✅ Real-time telemetry streaming (images, plots, status)
✅ IPython console with remote supervisor access
✅ Multiple GUIs can connect to same server (read-only monitoring)
✅ Automatic reconnection on network hiccups
✅ Low latency command execution (< 1ms on LAN)

### Troubleshooting Remote Mode

**"Connection Failed"**
- Check server is running
- Verify IP address and ports
- Check firewall settings
- Try `telnet <server> 5555` to test connectivity

**"Timeout" errors**
- Network latency too high
- Server overloaded
- Increase timeout in connection settings (future feature)

**Telemetry lag**
- Reduce GUI update frequency
- Decrease monitoring frequency on server
- Check network bandwidth with `iperf3`

**"Command socket error"**
- Server crashed or stopped
- Network interruption
- Click "Disconnect" and reconnect

## Extending the GUI

### Custom Commands

Send custom commands to the supervisor:

```python
# Access the main window instance
window.supervisor_thread.send_custom_command(
    'method_name',  # Method on supervisor
    arg1, arg2,     # Arguments
    key=value       # Keyword arguments
)
```

### Custom Monitoring

Connect to supervisor signals:

```python
def my_handler(se_strehl, le_strehl):
    print(f"SE: {se_strehl}, LE: {le_strehl}")
    # Log, save, analyze, etc.

window.supervisor_thread.strehl_updated.connect(my_handler)
```

### Available Signals

From `SupervisorThread`:
- `iteration_done(int)` - Iteration number
- `strehl_updated(float, float)` - SE and LE Strehl
- `framerate_updated(float)` - Current FPS
- `psf_updated(ndarray, int)` - PSF image and target index
- `wfs_image_updated(ndarray, int)` - WFS image and sensor index
- `phase_updated(ndarray)` - Phase screen
- `error_occurred(str)` - Error message
- `loop_finished()` - Loop completed
- `status_changed(str)` - Status update

## Troubleshooting

### Import Errors

```
ImportError: No module named 'PyQt6'
```

**Solution**: Install dependencies
```bash
pip install PyQt6 pyqtgraph
```

### GUI Freezes

**Cause**: Monitoring frequency too high or blocking operations in GUI thread

**Solution**: 
- Increase "GUI Update Freq" value
- Ensure heavy computations are in supervisor thread

### Display Not Updating

**Cause**: Signals not connected or monitoring disabled

**Solution**:
- Check that supervisor is initialized
- Verify monitoring frequency > 0
- Look for errors in terminal output

### Slow Performance

**Cause**: Too frequent GUI updates

**Solution**:
- Increase monitoring frequency (e.g., 50-100 frames)
- Close unused display tabs
- Reduce plot history (edit `max_points` in widgets)

## File Structure

```
shesha/
├── shesha/
│   ├── gui/
│   │   ├── __init__.py              # Package init
│   │   ├── main_window.py           # Main GUI window
│   │   ├── supervisor_thread.py     # Worker thread
│   │   ├── display_widgets.py       # Display components
│   │   ├── gui_example.py           # Usage examples
│   │   └── README.md                # Detailed docs
│   └── scripts/
│       └── compass_gui.py           # Launcher script
└── requirements-gui.txt             # GUI dependencies
```

## Advanced Usage

See `shesha/gui/gui_example.py` for examples of:
- Programmatic GUI creation
- Custom button actions
- Auto-initialization
- Custom monitoring and logging

## Support

- Documentation: `shesha/gui/README.md`
- Examples: `shesha/gui/gui_example.py`
- Issues: https://github.com/COSMIC-RTC/compass/issues

## License

LGPL v3+ (same as COMPASS)
Copyright (C) 2011-2024 COSMIC Team
