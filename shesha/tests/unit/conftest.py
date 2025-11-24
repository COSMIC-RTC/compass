"""
Pytest configuration for shesha unit tests.
"""

import os
import pytest
from pathlib import Path


@pytest.fixture
def shesha_root():
    """Get SHESHA_ROOT path."""
    if "SHESHA_ROOT" in os.environ:
        return Path(os.environ["SHESHA_ROOT"])
    else:
        # Try to find it relative to this file
        shesha_path = Path(__file__).parent.parent.parent
        if (shesha_path / "data").exists():
            return shesha_path
        return None


@pytest.fixture
def test_data_dir(shesha_root):
    """Get test data directory."""
    if shesha_root:
        return shesha_root / "data"
    return None


@pytest.fixture
def temp_config(tmp_path):
    """Create temporary configuration for tests."""
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    return config_dir
