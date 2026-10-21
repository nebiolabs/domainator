"""Suite-wide pytest configuration.

Markers (registered in pyproject.toml):

- ``requires_binary(*names)``: skip the test unless every named external program is on PATH.
- ``network``: the test downloads from the internet. Deselect with ``-m "not network"``.
"""

import logging
import shutil

import pytest


def pytest_runtest_setup(item):
    for marker in item.iter_markers(name="requires_binary"):
        missing = [name for name in marker.args if shutil.which(name) is None]
        if missing:
            pytest.skip(f"requires external program(s) not on PATH: {', '.join(missing)}")


@pytest.fixture(autouse=True)
def _reset_domainator_logger():
    """Undo the logging setup a tool's main() does, so a test's verbosity can't leak into the next one."""
    yield
    logger = logging.getLogger("domainator")
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
    logger.setLevel(logging.NOTSET)
