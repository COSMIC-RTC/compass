"""
Pytest configuration and fixtures for CLI tests.
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

try:
    import pytest
except ImportError:
    pytest = None


if pytest is not None:
    @pytest.fixture
    def temp_dir():
        """Create a temporary directory for tests."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def mock_compass_env(temp_dir, monkeypatch):
        """Mock COMPASS environment variables."""
        compass_root = temp_dir / "compass"
        compass_root.mkdir()
        (compass_root / "local").mkdir()
        (compass_root / "shesha").mkdir()
        (compass_root / "build").mkdir()
        (compass_root / "libcarma").mkdir()
        (compass_root / "libsutra").mkdir()
        (compass_root / "python_module").mkdir()
        
        monkeypatch.setenv("COMPASS_ROOT", str(compass_root))
        monkeypatch.setenv("COMPASS_INSTALL_ROOT", str(compass_root / "local"))
        monkeypatch.setenv("SHESHA_ROOT", str(compass_root / "shesha"))
        
        return {
            "root": compass_root,
            "install": compass_root / "local",
            "shesha": compass_root / "shesha",
        }

    @pytest.fixture
    def mock_home_env(temp_dir, monkeypatch):
        """Mock home directory for config tests."""
        home_dir = temp_dir / "home"
        home_dir.mkdir()
        
        monkeypatch.setenv("HOME", str(home_dir))
        
        return home_dir
