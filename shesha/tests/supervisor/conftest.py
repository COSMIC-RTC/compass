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
Pytest configuration for supervisor integration tests.

These tests instantiate a full CompassSupervisor which requires a live GPU.
The entire directory is skipped when no CUDA device is available so that
collection does not crash on the module-level supervisor initialisation.
"""

import sys
import pytest


def _gpu_available() -> bool:
    try:
        import carma  # noqa: F401
        ctx = carma.context.get_instance_1gpu(0)
        del ctx
        return True
    except Exception:
        return False


def pytest_ignore_collect(collection_path, config):
    """Skip collecting this directory when no GPU is available."""
    if not _gpu_available():
        return True
    return None
