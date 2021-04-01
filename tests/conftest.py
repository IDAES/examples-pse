"""
This file is automatically run by pytest.

The `pytest_generate_tests` function is a hook called during pytest's collection
process. We use this hook to have the "notebook" parameter include a list of all
the Jupyter notebooks under the top-level "src" directory. Given this,
any test that takes a parameter called "notebook" will actually be run on each of
these notebooks.

See https://docs.pytest.org/en/stable/parametrize.html#pytest-generate-tests
for more details on how this works.
"""
import os
from pathlib import Path

_root = os.path.join(os.path.dirname(__file__), "..")


def _find_notebooks(p: Path, file_list):
    """Find all notebooks under path `p`, and add them to the file_list.
    """
    # the return value of Path.iterdir() should be sorted to ensure consistency across different OSes
    for entry in sorted(p.iterdir()):
        if entry.is_dir():
            _find_notebooks(entry, file_list)
        elif entry.name.endswith(".ipynb"):
            file_list.append(str(entry))


def pytest_generate_tests(metafunc):
    """Parameterize tests in this directory.
    """
    all_notebooks = []
    _find_notebooks(Path(_root) / "src", all_notebooks)
    if "notebook" in metafunc.fixturenames:
        metafunc.parametrize("notebook", all_notebooks)