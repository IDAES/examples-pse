"""
Test notebooks
"""
# stdlib
import json
import logging
import os
from pathlib import Path
import re
import subprocess
import sys

# third-party
import nbformat
import pytest
import yaml

# import "build" module from top dir
_root = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, _root)
import build

newline = "\n"  # useful for f-strings


@pytest.fixture(scope="module")
def settings_ci():
    os.chdir(_root)
    return build.Settings(open("build-ci.yml", "r"))


@pytest.mark.unit
def test_parse_notebook(notebook):
    """The parameter 'notebook' is parameterized in `conftest.py`, so that
    this test is called for every Jupyter notebook found under the "src/" dir.
    """
    nbformat.read(notebook, as_version=build.ParallelNotebookWorker.JUPYTER_NB_VERSION)


@pytest.mark.integration
def test_run_all_notebooks():
    os.chdir(_root)
    # make sure we've cleaned up first
    cmd = ["python", "build.py", "--remove"]
    proc = subprocess.Popen(cmd)
    proc.wait()
    assert proc.returncode == 0
    # now run
    cmd = ["python", "build.py", "--config", get_build_config(), "--exec"]
    proc = subprocess.Popen(cmd)
    proc.wait()
    assert proc.returncode == 0
    find_broken_links()


@pytest.mark.component
def test_broken_links():
    find_broken_links()


def find_broken_links(rebuild=True):
    """Run the Sphinx link checker.

    This was created in response to a number of broken links in Jupyter notebook
    cells, but would also find broken links in any documentation pages.
    """
    os.chdir(_root)
    config = get_build_config()
    config_dict = load_build_config(config)
    # Copy notebooks to docs. -S suppresses Sphinx output.
    args = ["python", "build.py", "--config", config, "-Sy"]
    proc = subprocess.Popen(args)
    rc = proc.wait()
    assert rc == 0, "Copying notebooks to docs failed"
    # Run linkchecker (-l). -S suppresses Sphinx output.
    # output will be in dir configured in sphinx.linkcheck_dir (see below)
    proc = subprocess.Popen(["python", "build.py", "--config", config, "-Sl"])
    rc = proc.wait()
    assert rc == 0, "Linkchecker process failed"
    # find links marked [broken], report them
    link_file = Path(".") / config_dict["sphinx"]["linkcheck_dir"] / "output.json"
    assert link_file.exists()
    links = []
    for line in link_file.open(mode="r", encoding="utf-8"):
        obj = json.loads(line)
        if obj["status"] == "broken":
            num = len(links) + 1
            links.append(f"{num}) {obj['filename']}:{obj['lineno']} -> {obj['uri']}")
    # fail if there were any broken links
    assert len(links) == 0, f"{len(links)} broken links:\n" f"{newline.join(links)}"


def test_index_page():
    config = get_build_config()
    config_dict = load_build_config(config)
    args = ["python", "build.py", "--config", config, "--index", "--index-dev"]
    print(f"Build index page (in 'dev' mode) with command: {' '.join(args)}")
    proc = subprocess.Popen(args)
    rc = proc.wait()
    assert rc == 0, "Failed to build Jupyter notebook index page"

# Utility


def get_build_config():
    if os.environ.get("GITHUB_ACTIONS", False):
        return "build-ci.yml"
    return "build.yml"


def load_build_config(config):
    return yaml.safe_load(open(config, "r"))