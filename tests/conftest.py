import os
from pathlib import Path

_root = os.path.join(os.path.dirname(__file__), "..")


def _find_notebooks(p: Path, file_list):
    """Find all notebooks under path `p`, and add them to the file_list.
    """
    for entry in p.iterdir():
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