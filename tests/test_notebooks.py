"""
Test notebooks
"""
from pathlib import Path
import subprocess

# assume parent of this test's directory is top
TOP_DIR = Path(__file__).parent.parent


def test_build_notebooks():
    cmd = str(TOP_DIR / "build.py")
    cmdargs = [cmd, "--no-sphinx", "--config",
               str(TOP_DIR / "build.yml"), "-vv"]
    proc = subprocess.run(cmdargs)
    assert proc.returncode == 0
