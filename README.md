# COMPASS

Main status:
[![Main status](https://gitlab.obspm.fr/cosmic-rtc/compass/badges/main/pipeline.svg)](https://gitlab.obspm.fr/cosmic-rtc/compass/commits/main)

Develop status:
[![Develop status](https://gitlab.obspm.fr/cosmic-rtc/compass/badges/develop/pipeline.svg)](https://gitlab.obspm.fr/cosmic-rtc/compass/commits/develop)
[![coverage report](https://gitlab.obspm.fr/cosmic-rtc/compass/badges/develop/coverage.svg)](https://cosmic-rtc.pages.obspm.fr/compass/coverage/index.html)
[![Code Quality](https://gitlab.obspm.fr/cosmic-rtc/compass/badges/develop/code_quality.svg)](https://gitlab.obspm.fr/cosmic-rtc/compass/pipelines)

- [COMPASS](#compass)
  - [Citations](#citations)
  - [Overview](#overview)
    - [Hardware requirements](#hardware-requirements)
    - [Environment requirements](#environment-requirements)
  - [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Step 1: Download COMPASS](#step-1-download-compass)
    - [Step 2: Setup Python Environment](#step-2-setup-python-environment)
    - [Step 3: Activate the Environment](#step-3-activate-the-environment)
    - [Step 4: Install COMPASS CLI](#step-4-install-compass-cli)
    - [Step 5: Initialize Modulefiles](#step-5-initialize-modulefiles)
    - [Step 6: Build COMPASS](#step-6-build-compass)
  - [Usage](#usage)
    - [Quick Start](#quick-start)
    - [Simulation Commands](#simulation-commands)
    - [Other Useful Commands](#other-useful-commands)
    - [Configuration](#configuration)
  - [Contributing](#contributing)
  - [License](#license)

## Citations
If you use COMPASS in your research, please cite one of the following papers:

- [Ferreira, F. et al, “COMPASS: an efficient GPU-based simulation software for adaptive optics systems”, HPCS 2018](https://doi.org/10.1109/HPCS.2018.00043)
- [Ferreira, F. et al., “Real-time end-to-end AO simulations at ELT scale on multiple GPUs with the COMPASS platform”, SPIE 2018](https://doi.org/10.1117/12.2312593)
- [Gratadour, D. et al, "COMPASS: an efficient, scalable and versatile numerical platform for the development of ELT AO systems"](https://doi.org/10.1117/12.2056358)
- [Gratadour, D. et al, “GPUs for adaptive optics: simulations and real-time control”, SPIE, 2012](https://doi.org/10.1117/12.925723)

## Overview

The COMPASS platform is distributed as a single bundle of CArMA and SuTrA C++ / Cuda libraries and their Python extension SHESHA.

### Hardware requirements

The system must contain at least an x86 CPU and a CUDA capable GPU. list of compatible GPUs can be found here <http://www.nvidia.com/object/cuda_gpus.html>. Specific requirements apply to clusters (to be updated).

### Environment requirements

The system must be running a 64 bit distribution of Linux with the latest NVIDIA drivers and CUDA toolkit. The following installation instructions are valid if the default installation paths have been selected for these components.

## Installation

### Prerequisites

Before installing COMPASS, ensure you have:
- **NVIDIA CUDA Toolkit** installed (compatible GPU required)
- **Lmod** (Environment Modules) installed on your system
  - Install via package manager: `sudo apt install lmod` (Debian/Ubuntu) or `sudo yum install Lmod` (RHEL/CentOS)
  - More info: <https://lmod.readthedocs.io/>
- **Git LFS** (Large File Storage) installed
  - Install via package manager: `sudo apt install git-lfs` (Debian/Ubuntu) or `sudo dnf install git-lfs` (RHEL/CentOS)
  - Initialize: `git lfs install`
- **Python 3** 

### Step 1: Download COMPASS

Clone the repository:

```bash
git clone https://gitlab.obspm.fr/cosmic-rtc/compass.git
cd compass
```

### Step 2: Install COMPASS CLI

```bash
pip install -e .
```

### Step 3: Initialize and build COMPASS

```bash
compass init
compass build
```

This will:
- Adds COMPASS modulefiles directory to your `MODULEPATH`
- Updates your `~/.bashrc` with the necessary configuration
- Enables `module load compass/local` for future sessions
- Setup the compass modulefile 
- Build **libcarma** and **libsutra** C++/CUDA libraries
- Build Python extensions (**carma.so**, **sutra.so**)
- Install everything to `local/` directory

You can check the installation status at any time with:
```bash
compass check
```

## Usage

COMPASS provides a powerful command-line interface for running AO simulations. The main simulation commands are available through `compass sim`.

### Quick Start

1. **List available parameter files:**
   ```bash
   compass sim list
   ```
   This shows all parameter directories in `shesha/data/par/`

2. **List parameter files in a specific directory:**
   ```bash
   compass sim list MICADO
   ```
   Shows all `.py` parameter files in the MICADO directory

3. **Run a simulation (interactive CLI mode):**
   ```bash
   compass sim run shesha/data/par/MICADO/micado_full.py
   ```
   Launches an interactive IPython session with your simulation

4. **Run a GUI simulation:**
   ```bash
   compass sim gui shesha/data/par/MICADO/micado_full.py
   ```
   Launches the graphical interface for your simulation

### Simulation Commands

#### `compass sim run` - Interactive CLI Simulations

Run simulations in an interactive IPython session with the default script (default: `closed_loop.py`):

```bash
# Basic usage
compass sim run <parameter_file>

# With options
compass sim run shesha/data/par/MICADO/micado_full.py \
  --iterations 1000 \
  --devices 0,1

# Pass additional arguments to the script
compass sim run parfile.py -- --custom-arg value
```

**Options:**
- `--iterations N` - Set number of iterations
- `--devices 0,1` - Specify GPU devices (comma-separated)
- `-- <args>` - Pass additional arguments to the simulation script

#### `compass sim gui` - GUI Simulations

Run simulations with a graphical interface using the default GUI script (default: `widget_ao.py`):

```bash
# Basic usage
compass sim gui <parameter_file>

# With options
compass sim gui shesha/data/par/MICADO/micado_full.py \
  --iterations 1000 \
  --devices 0
```

#### `compass sim script` - Manage CLI Scripts

Configure which script to use for CLI simulations:

```bash
# Show current default script
compass sim script show

# List available scripts
compass sim script list

# Set default script (by name from shesha/scripts/)
compass sim script set closed_loop.py

# Set custom script (absolute path)
compass sim script set /path/to/my_custom_script.py
```

Available CLI scripts in `shesha/scripts/`:
- `closed_loop.py` - Standard closed-loop AO simulation
- `cosmic_simulator.py` - COSMIC instrument simulator
- `dm_standalone.py` - Deformable mirror standalone test
- `micado_loop.py` - MICADO-specific simulation loop

#### `compass sim gui-script` - Manage GUI Scripts

Configure which widget to use for GUI simulations:

```bash
# Show current GUI script
compass sim gui-script show

# List available GUI widgets
compass sim gui-script list

# Set default GUI script
compass sim gui-script set widget_ao.py
```

Available GUI widgets in `shesha/widgets/`:
- `widget_ao.py` - Standard AO widget (default)
- `widget_bench.py` - Bench testing widget
- `widget_ao_expert.py` - Expert mode AO widget
- `widget_canapass.py` - CANAPASS instrument widget
- `widget_cosmic_simulator.py` - COSMIC simulator widget
- `widget_twoStages.py` - Two-stage AO widget

### Other Useful Commands

```bash
# Check system status (packages, modules, builds)
compass check

# Build specific components
compass build libcarma
compass build libsutra
compass build python_module

# Show current configuration
compass config show

# Display version
compass --version
```

### Configuration

COMPASS stores simulation preferences in `~/.compass/sim_config.json`:
- `default_script` - Default CLI simulation script (absolute path)
- `default_gui_script` - Default GUI simulation script (absolute path)

You can edit this file manually or use the `compass sim script` and `compass sim gui-script` commands.

## Contributing

If you want to contribute to COMPASS, please read the [CONTRIBUTING.md](CONTRIBUTING.md) file.
Contributions must go through merge requests.

## License

COMPASS is distributed under the LGPLv3 license. See [LICENSE](LICENSE) for more information.