"""
Test notebooks
"""
# stdlib
import logging
import os
import sys

# third-party
import nbformat
import pytest

# import "build" module from top dir
_root = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, _root)
import build


@pytest.fixture(scope="module")
def settings():
    os.chdir(_root)
    settings = build.Settings(open("build.yml", "r"))


@pytest.mark.integration
def test_build_notebooks(settings):
    os.chdir(_root)
    settings = build.Settings(open("build.yml", "r"))
    nb = build.NotebookBuilder(settings)
    res = nb.build(
        {
            "rebuild": True,
            "continue_on_error": True,
            "test_mode": True,
            "error_file": "__stdout__",
        }
    )
    assert res.n_success > 0  # something ran
    assert res.n_fail == 0  # nothing failed


@pytest.mark.component
def test_convert_notebooks():
    build._log.setLevel(logging.INFO)  # otherwise DEBUG for some reason
    print(f"@@ log level = {build._log.getEffectiveLevel()}")
    os.chdir(_root)
    settings = build.Settings(open("circleci-test.yml", "r"))
    nb = build.NotebookBuilder(settings)
    nb.build({"rebuild": True})
    total, num_failed = nb.report()
    assert total > 0
    assert num_failed == 0


@pytest.mark.unit
def test_parse_notebook(notebook):
    """The parameter 'notebook' is parameterized in conftest.py, so that
    this test is called for every Jupyter notebook found under the "src/" dir.
    """
    nbformat.read(notebook, as_version=build.NotebookBuilder.JUPYTER_NB_VERSION)
