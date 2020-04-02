"""
Test notebooks
"""
# stdlib
import os
import sys

# third-party
import pytest

# import "build" module from top dir
_root = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, _root)
import build


@pytest.mark.integration
def test_build_notebooks():
    os.chdir(_root)
    settings = build.Settings(open("build.yml", "r"))
    nb = build.NotebookBuilder(settings)
    res = nb.build({"execute": True, "rebuild": True, "continue": True,
                    "test_mode": True, "error_file": "__stdout__"})
    assert res.total > 0  # something ran
    assert len(res.failed) == 0  # nothing failed
