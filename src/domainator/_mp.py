"""Multiprocessing with a start method that is safe for Domainator's workers.
"""
import multiprocessing
import os
import warnings


_MP_CONTEXT = None


def get_multiprocessing_context():
    """Return a cached multiprocessing context that avoids unsafe fork start."""

    global _MP_CONTEXT
    if _MP_CONTEXT is not None:
        return _MP_CONTEXT

    preferred_method = os.environ.get("DOMAINATOR_MP_START_METHOD", "spawn")
    try:
        ctx = multiprocessing.get_context(preferred_method)
    except ValueError:
        warnings.warn(
            f"Requested multiprocessing start method '{preferred_method}' is unavailable; "
            "falling back to Python's default start method.",
            RuntimeWarning,
            stacklevel=2,
        )
        ctx = multiprocessing.get_context()

    _MP_CONTEXT = ctx
    return ctx


def make_pool(*args, **kwargs):
    """Create a multiprocessing Pool using the safe start method."""

    return get_multiprocessing_context().Pool(*args, **kwargs)


def make_manager():
    """Create a multiprocessing Manager using the safe start method."""

    return get_multiprocessing_context().Manager()
