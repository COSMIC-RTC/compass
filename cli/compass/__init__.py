"""
COMPASS - COMputing Platform for Adaptive optics SystemS
============================================================================

A modern Python-based deployment and management framework for COMPASS AO systems.

COMPASS is an adaptive optics simulation platform that includes:
- CARMA: CUDA-based AO Real-time Modules and Applications
- SUTRA: Simulation Utilities for Tomographic Reconstruction of Adaptive optics
- Shesha: Python interface and high-level simulation control

Modules
-------
- compass.core: Configuration and environment management
- compass.deployment: System deployment and build management
- compass.simulation: Simulation control and monitoring
"""

__version__ = "6.2.0"
__author__ = "COMPASS Team"
__license__ = "LGPL-3.0-or-later"

from compass.core.config import Config, Environment
from compass.core.logger import setup_logger

__all__ = [
    "__version__",
    "__author__",
    "__license__",
    "Config",
    "Environment",
    "setup_logger",
]
