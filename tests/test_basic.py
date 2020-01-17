"""
Basic sanity-check unit tests
"""
# stdlib
import importlib.util
import os
from pathlib import Path
import sys
import time
from typing import List
# third-party
import pytest

sep = os.path.sep
nthreads = 4

class BadSyntax(Exception):
    pass

# Tests


@pytest.mark.unit
def test_import_syntax():
    srcdir = Path(__file__).parent.parent / "src"
    items = find_python_modules(srcdir)
    errs = []
    for item in items:
        try:
            import_item(item)
        except BadSyntax as err:
            errs.append(str(err))
    if errs:
        errmsg = "\n".join(errs)
        assert False, f"Syntax errors encountered: {errmsg}"


def import_item(item):
    try:
        import_python_file(*item)
    # report syntax errors
    except SyntaxError as err:
        raise BadSyntax(f"[{type(err)}] {item[0]}: {err}")
    # ignore everything else
    except Exception:
        pass

# Helper functions


def find_python_modules(target_dir: Path) -> List:
    """Find all python modules from target_dir, on down, that contain a
    Python module or sub-package.
    """
    tgtdir_len = len(str(target_dir))
    modules = []
    for path in target_dir.rglob("*.py"):
        spath = str(path)
        if spath.endswith("__init__.py"):
            continue
        if ".ipynb_checkpoints" in spath:
            continue
        name = module_name_from_path(path, tgtdir_len)
        modules.append((path, name))
    return modules


def module_name_from_path(p, n):
    # strip leading target dir and separator
    s = str(p)[n + 1:]
    # strip ext
    s = s[:-3]
    # change separator to .
    s = s.replace(sep, ".")
    return s


def import_python_file(file_path, module_name):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
