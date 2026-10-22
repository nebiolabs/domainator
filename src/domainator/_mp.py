"""Pools of workers: threads when Python runs without the GIL, otherwise processes.

On free-threaded Python with the GIL off (every loaded extension supports running without it),
workers are threads, which share memory. Otherwise they are processes started with spawn, the
one start method that is safe on every platform. ``DOMAINATOR_MP_START_METHOD`` overrides the
choice: ``threads``, or a multiprocessing start method (``spawn``, ``forkserver``, ``fork``).

Pool tasks get large read-only data (references, matrices) through :class:`Shared`, which works
for both: threads read the parent's copy, and each worker process gets its own copy once.
"""
import contextlib
import multiprocessing
import multiprocessing.pool
import os
import queue
import sys
import threading
import uuid
import warnings
from typing import Any, Callable, Iterable, Iterator, Optional, Sequence


def _start_method() -> str:
    """``threads``, or the multiprocessing start method that pools use."""
    method = os.environ.get("DOMAINATOR_MP_START_METHOD")
    if method:
        return method
    gil_enabled = getattr(sys, "_is_gil_enabled", lambda: True)()
    return "spawn" if gil_enabled else "threads"


def _process_context():
    method = _start_method()
    if method == "threads":
        method = "spawn"
    try:
        return multiprocessing.get_context(method)
    except ValueError:
        warnings.warn(
            f"Requested multiprocessing start method '{method}' is unavailable; "
            "falling back to Python's default start method.",
            RuntimeWarning,
            stacklevel=3,
        )
        return multiprocessing.get_context()


class Shared:
    """A read-only value for pool tasks, given to each worker once rather than with every task.

    A pool sends a task's function and arguments to a worker process with every task, so a
    large value held by them would be copied again and again. Hold it in a Shared instead,
    pass that to :func:`make_pool` (``shared=``), and read it in the task with :meth:`get`. A
    Shared pickles as just its key, and the pool gives each worker process the value once.
    Threads read the parent's value directly.

    Make it with the value, or with ``build``, a picklable function of no arguments that makes
    the value. Each worker process then builds its own copy (once) instead of being sent it,
    for values that are expensive to pickle; threads share one copy, built on first use.
    """

    def __init__(self, value: Any = None, *, build: Optional[Callable[[], Any]] = None):
        self._key = uuid.uuid4().hex
        self._value = value
        self._build = build
        self._built = build is None
        self._lock = threading.Lock()

    def get(self) -> Any:
        if not self._built:
            with self._lock:
                if not self._built:
                    self._value = self._build()
                    self._built = True
        return self._value

    def _for_worker(self):
        """What a worker process is sent once: the key, and the value or how to build it."""
        if self._build is not None:
            return self._key, None, self._build
        return self._key, self._value, None

    def __reduce__(self):
        return (_InstalledShared, (self._key,))


# The Shared values of this worker process, by key (filled by make_pool's initializer).
_INSTALLED = {}


class _InstalledShared:
    """A Shared as a task sees it in a worker process."""

    def __init__(self, key):
        self._key = key

    def get(self):
        return _INSTALLED[self._key]


def _install_shared(entries):
    for key, value, build in entries:
        _INSTALLED[key] = build() if build is not None else value


@contextlib.contextmanager
def make_pool(processes: int, *, shared: Sequence[Shared] = ()) -> Iterator[multiprocessing.pool.Pool]:
    """A pool of ``processes`` workers (threads or processes, see the module docstring).

    ``shared`` lists the :class:`Shared` values the tasks read. Leaving the ``with`` block stops
    the workers.
    """
    if _start_method() == "threads":
        pool = multiprocessing.pool.ThreadPool(processes)
    else:
        pool = _process_context().Pool(processes, initializer=_install_shared,
                                       initargs=(tuple(s._for_worker() for s in shared),))
    with pool:
        yield pool


def make_manager():
    """Create a multiprocessing Manager using the safe start method."""

    return _process_context().Manager()


class _Failed:
    def __init__(self, error):
        self.error = error


def imap_unordered_bounded(pool, function: Callable, items: Iterable, max_pending: int) -> Iterator:
    """``pool.imap_unordered(function, items)``, but reading ``items`` only ``max_pending`` tasks ahead.

    ``imap_unordered`` reads all of ``items`` at once, which for a stream of records would hold
    the whole stream in memory. Results come back in the order they finish; the first task to
    fail raises its exception here.
    """
    finished = queue.SimpleQueue()
    pending = 0

    def wait_for_one():
        result = finished.get()
        if isinstance(result, _Failed):
            raise result.error
        return result

    for item in items:
        pool.apply_async(function, (item,), callback=finished.put,
                         error_callback=lambda error: finished.put(_Failed(error)))
        pending += 1
        if pending >= max_pending:
            yield wait_for_one()
            pending -= 1
    while pending:
        yield wait_for_one()
        pending -= 1


def resolve_cpu(cpu: int) -> int:
    """Thread count for a ``cpu`` value: ``0`` means one per physical core available to this process."""
    if cpu > 0:
        return cpu
    if cpu < 0:
        raise ValueError(f"--cpu must be 0 (all cores) or a positive number of threads, not {cpu}.")
    try:
        available = len(os.sched_getaffinity(0))
    except AttributeError:  # not on Linux
        available = os.cpu_count() or 1
    try:
        import psutil
        physical = psutil.cpu_count(logical=False)
    except ImportError:
        physical = None
    return max(1, min(available, physical or available))
