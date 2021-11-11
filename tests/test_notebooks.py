"""
Test notebooks
"""
# stdlib
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


@pytest.mark.component
def test_convert_some_notebooks(settings_ci):
    build._log.setLevel(logging.DEBUG)  # otherwise DEBUG for some reason
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
    cmd = ["python", "build.py", "--config", get_build_config(), "--test"]
    proc = subprocess.Popen(cmd)
    proc.wait()
    assert proc.returncode == 0
    find_broken_links(permissive=False)


@pytest.mark.component
def test_broken_links():
    find_broken_links(permissive=True)


def find_broken_links(permissive=True):
    """Run the Sphinx link checker.

    This was created in response to a number of broken links in Jupyter notebook
    cells, but would also find broken links in any documentation pages.

    For it to be useful, you need to have run/converted all the notebooks.
    This function is called at the end of the main notebook integration test.
    """
    os.chdir(_root)
    config = get_build_config()
    config_dict = yaml.safe_load(open(config, "r"))
    docs_root = Path(_root) / config_dict["paths"]["output"]
    # verify that notebooks are copied into the docs tree
    empty_dirs, num_subdirs, empty_dir_paths = 0, 0, []
    for subdir in config_dict["notebook"]["directories"]:
        num_subdirs += 1
        subdir_name = subdir["source"]
        #  debug print(f"Look in {str(docs_root / subdir_name)}")
        if len(list((docs_root / subdir_name).glob("*.rst"))) <= 1:
            empty_dirs += 1
            empty_dir_paths.append(str(docs_root / subdir_name))
    # print warnings, but only fail if there are NO notebooks at all
    if empty_dirs > 0:
        lvl = "WARNING" if empty_dirs < num_subdirs else "ERROR"
        print(f"{lvl}: {empty_dirs}/{num_subdirs} directories did not have notebooks")
        print(
            "Perhaps you need to run (in the repo root):\n\n"
            "    python build.py -cd\n\n"
            "This executes the Jupyter Notebooks in 'src'\n"
            "and copies them into the 'docs' directory tree."
        )
    if permissive:
        # continue if there are some non-empty dirs, skip if there are
        # no non-empty dirs
        if empty_dirs == num_subdirs:
            pytest.skip("No notebooks in any directories")
    else:
        assert empty_dirs == 0, (
            f"Notebooks are missing in some directories:\n"
            f"{newline.join(empty_dir_paths)}"
        )
    # run linkchecker, -S means suppress normal Sphinx output.
    # output will be in dir configured in sphinx.linkcheck_dir (see below)
    proc = subprocess.Popen(["python", "build.py", "--config", config, "-Sl"])
    rc = proc.wait()
    assert rc == 0, "Linkchecker process failed"
    # find links marked [broken], report them
    link_file = Path(".") / config_dict["sphinx"]["linkcheck_dir"] / "output.txt"
    assert link_file.exists()
    links = []
    for line in link_file.open(mode="r", encoding="utf-8"):
        m = re.search(r"^([^:]*):(\d+):.*\[broken]\s*(https?://[^:]*)", line)
        if m:
            num = len(links) + 1
            links.append(f"{num}) {m.group(1)}:{m.group(2)} -> {m.group(3)}")
    # fail if there were any broken links
    assert len(links) == 0, f"{len(links)} broken links:\n" f"{newline.join(links)}"


def get_build_config():
    if os.environ.get("GITHUB_ACTIONS", False):
        return "build-ci.yml"
    return "build.yml"
