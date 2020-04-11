"""Test timeout on Lock.acquire()."""
from ase.utils import Lock
from ase.test import must_raise

lock = Lock('lockfile', timeout=0.3)
with lock:
    with must_raise(TimeoutError):
        with lock:
            ...
