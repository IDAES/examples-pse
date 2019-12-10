"""
Test notebooks
"""
import os
from pathlib import Path
import subprocess

# assume parent of this test's directory is top
TOP_DIR = Path(__file__).parent.parent


def test_build_notebooks():
    ensure_build_directory()
    cmd = str(TOP_DIR / "build.py")
    cmdargs = [cmd, "--no-sphinx", "--config",
               str(TOP_DIR / "build.yml"), "-vv"]
    proc = subprocess.run(cmdargs)
    assert proc.returncode == 0


def ensure_build_directory():
    """In case docs aren't built, ensure the docs
    'build' directory exists; this is needed because
    the notebook creation tries to copy images in there.
    """
    build_dir = TOP_DIR / "docs" / "build"
    if not build_dir.exists():
        os.mkdir(str(build_dir))
        print("Created empty docs build dir at: {build_dir}")
