"""
Test notebooks
"""
# stdlib
import logging
import os
import subprocess
import sys

# third-party
import nbformat
import pytest

# import "build" module from top dir
_root = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, _root)
import build


@pytest.fixture(scope="module")
def settings_ci():
    os.chdir(_root)
    return build.Settings(open("build-circleci.yml", "r"))


@pytest.mark.component
def test_convert_some_notebooks(settings_ci):
    build._log.setLevel(logging.INFO)  # otherwise DEBUG for some reason
    os.chdir(_root)
    nb = build.NotebookBuilder(settings_ci)
    nb.build({"rebuild": True})
    total, num_failed = nb.report()
    assert total > 0
    assert num_failed == 0


@pytest.mark.unit
def test_parse_notebook(notebook):
    """The parameter 'notebook' is parameterized in `conftest.py`, so that
    this test is called for every Jupyter notebook found under the "src/" dir.
    """
    nbformat.read(notebook, as_version=build.NotebookBuilder.JUPYTER_NB_VERSION)


@pytest.mark.integration
def test_run_all_notebooks():
    os.chdir(_root)
    # make sure we've cleaned up first
    cmd = ["python", "build.py", "--remove"]
    proc = subprocess.Popen(cmd)
    proc.wait()
    assert proc.returncode == 0
    # now run
    cmd = ["python", "build.py", "--test"]
    proc = subprocess.Popen(cmd)
    proc.wait()
    assert proc.returncode == 0

