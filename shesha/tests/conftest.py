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
Top-level pytest configuration for shesha tests.

Provides GPU detection and automatic skipping of tests marked with
@pytest.mark.gpu when no CUDA-capable device is available.

Files listed in _GPU_ONLY_SCRIPTS contain module-level GPU initialisation
and are ignored entirely during collection when no GPU is present.
"""

import pytest

# Files (relative to this conftest) that run GPU code at module level and
# therefore cannot be safely imported during pytest collection without a GPU.
_GPU_ONLY_SCRIPTS = {"test_fp16.py"}


def _gpu_available() -> bool:
    """Return True if at least one CUDA-capable GPU can be initialised."""
    try:
        import carma  # noqa: F401 – the pybind11 extension
        ctx = carma.context.get_instance_1gpu(0)
        del ctx
        return True
    except Exception:
        return False


GPU_AVAILABLE = _gpu_available()


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "gpu: marks tests that require a CUDA-capable GPU (skip with -m 'not gpu')",
    )


def pytest_ignore_collect(collection_path, config):
    """Skip GPU-only scripts when no GPU is available."""
    if not GPU_AVAILABLE and collection_path.name in _GPU_ONLY_SCRIPTS:
        return True
    return None


def pytest_collection_modifyitems(config, items):
    """Automatically skip gpu-marked tests when no GPU is available."""
    if GPU_AVAILABLE:
        return
    skip_no_gpu = pytest.mark.skip(reason="No CUDA-capable GPU available")
    for item in items:
        if item.get_closest_marker("gpu"):
            item.add_marker(skip_no_gpu)
