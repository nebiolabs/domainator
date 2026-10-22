"""Build the stable-ABI (abi3) wheel, except on free-threaded Python, which it doesn't support.

Everything else is in pyproject.toml. The stable ABI is incompatible with free-threaded builds
before Python 3.15, so those get a wheel for their own version (cp314t, ...). PyO3 builds that
on its own: its abi3 feature (rust/Cargo.toml) has no effect on a free-threaded interpreter.
"""
import sysconfig

from setuptools import setup

if sysconfig.get_config_var("Py_GIL_DISABLED"):
    setup()
else:
    setup(options={"bdist_wheel": {"py_limited_api": "cp312"}})
