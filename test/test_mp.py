"""The worker pools in domainator._mp, with processes and with threads."""
import functools
import os
import pickle
import threading

import pytest

from domainator import _mp
from domainator._mp import Shared, imap_unordered_bounded, make_pool


@pytest.fixture(params=["spawn", "threads"])
def start_method(request, monkeypatch):
    monkeypatch.setenv("DOMAINATOR_MP_START_METHOD", request.param)
    return request.param


def _sum_shared(shared, offset):
    return sum(shared.get()) + offset


def _build_and_note_pid():
    return os.getpid()


def _shared_value(shared, _):
    return shared.get()


def _fail_on_three(x):
    if x == 3:
        raise ValueError("three")
    return x


def test_shared_values_reach_the_workers(start_method):
    shared = Shared(list(range(10)))
    with make_pool(2, shared=[shared]) as pool:
        assert sorted(pool.imap(functools.partial(_sum_shared, shared), range(4))) == [45, 46, 47, 48]


def test_a_shared_pickles_as_its_key():
    shared = Shared(list(range(100000)))
    assert len(pickle.dumps(shared)) < 200


def test_built_once_per_worker_process(monkeypatch):
    monkeypatch.setenv("DOMAINATOR_MP_START_METHOD", "spawn")
    shared = Shared(build=_build_and_note_pid)
    with make_pool(2, shared=[shared]) as pool:
        # Each worker process builds its own value (its pid), so a task sees its own process.
        builders = set(pool.map(functools.partial(_shared_value, shared), range(20), chunksize=1))
    assert 1 <= len(builders) <= 2
    assert os.getpid() not in builders


def test_threads_build_one_copy(monkeypatch):
    monkeypatch.setenv("DOMAINATOR_MP_START_METHOD", "threads")
    calls = []
    barrier = threading.Barrier(4)

    def build():
        calls.append(None)
        return object()

    shared = Shared(build=build)

    def task(_):
        barrier.wait()
        return id(shared.get())

    with make_pool(4) as pool:
        assert len(set(pool.map(task, range(4)))) == 1
    assert len(calls) == 1


def test_threads_unless_the_gil_is_on(monkeypatch):
    monkeypatch.delenv("DOMAINATOR_MP_START_METHOD", raising=False)
    monkeypatch.setattr(_mp.sys, "_is_gil_enabled", lambda: False, raising=False)
    assert _mp._start_method() == "threads"
    monkeypatch.setattr(_mp.sys, "_is_gil_enabled", lambda: True, raising=False)
    assert _mp._start_method() == "spawn"
    monkeypatch.setenv("DOMAINATOR_MP_START_METHOD", "forkserver")
    assert _mp._start_method() == "forkserver"


def test_bounded_imap_reads_only_ahead_by_max_pending(monkeypatch):
    monkeypatch.setenv("DOMAINATOR_MP_START_METHOD", "threads")
    read = 0

    def items():
        nonlocal read
        for i in range(100):
            read += 1
            yield i

    with make_pool(2) as pool:
        results = imap_unordered_bounded(pool, abs, items(), max_pending=3)
        first = next(results)
        assert read <= 3
        assert sorted([first, *results]) == list(range(100))


def test_bounded_imap_returns_every_result(start_method):
    with make_pool(2) as pool:
        assert sorted(imap_unordered_bounded(pool, abs, range(-20, 0), max_pending=4)) == list(range(1, 21))


def test_bounded_imap_raises_a_task_error(start_method):
    with make_pool(2) as pool:
        with pytest.raises(ValueError, match="three"):
            list(imap_unordered_bounded(pool, _fail_on_three, range(10), max_pending=2))


def test_hmm_database_gives_each_thread_its_own_profiles(shared_datadir):
    # HMMER reconfigures an optimized profile for each target's length while it searches, so
    # threads must not share them.
    from domainator.domainate import HMMDatabase
    from domainator.hmmio import read_hmms
    database = HMMDatabase(next(iter(read_hmms([str(shared_datadir / "CcdB.hmm")]).values())))
    mine = database.optimized_profiles()
    assert database.optimized_profiles() is mine
    other = []
    thread = threading.Thread(target=lambda: other.append(database.optimized_profiles()))
    thread.start()
    thread.join()
    assert other[0] is not mine and len(other[0]) == len(mine)

    copy = pickle.loads(pickle.dumps(database))
    assert dict(copy).keys() == dict(database).keys()
    assert len(copy.optimized_profiles()) == len(mine)
