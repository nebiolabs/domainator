"""Suite-wide pytest configuration.

Markers (registered in pyproject.toml):

- ``requires_binary(*names)``: skip the test unless every named external program is on PATH.
- ``network``: the test downloads from the internet. Deselect with ``-m "not network"``.
"""

import shutil

import pytest


def pytest_runtest_setup(item):
    for marker in item.iter_markers(name="requires_binary"):
        missing = [name for name in marker.args if shutil.which(name) is None]
        if missing:
            pytest.skip(f"requires external program(s) not on PATH: {', '.join(missing)}")
