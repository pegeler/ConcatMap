"""
Utility functions and classes.
"""
import os
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def chdir_temporarily(path: Path):
    cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)
