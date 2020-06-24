"""
Test the build.py script.

This imports and uses that script like a module of course.
"""
import json
import logging
from pathlib import Path
import sys

# third-party
import pytest

# assume we are at ROOT/tests/something.py and the script is at ROOT/build.py
sys.path.insert(0, str(Path(__file__).parent.parent.absolute()))

print(f"path: {sys.path}")
import build


@pytest.mark.unit
def test_notify():
    build.notify("hello")
    build.notify("hello", level=1)
    build.notify("hello", level=10)


@pytest.mark.unit
def test_script_usage():
    sys.argv = ["build.py", "--usage"]
    exit_code = build.main()
    assert exit_code == 0
