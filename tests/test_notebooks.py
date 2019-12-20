"""
Test notebooks
"""
import os
from pathlib import Path
import subprocess

# assume parent of this test's directory is top
TOP_DIR = Path(__file__).parent.parent

# optionally, control the Jupyter Notebook kernel
_kernel_envvar = "IDAES_KERNEL"
if _kernel_envvar in os.environ:
    g_kernel = os.environ[_kernel_envvar]
else:
    g_kernel = None


def test_build_notebooks():
    ensure_build_directory()
    cmd = str(TOP_DIR / "build.py")
    cmdargs = [cmd, "--no-sphinx", "--config",
               str(TOP_DIR / "build.yml"), "-vv"]
    if g_kernel:
        cmdargs.append(f"--kernel={g_kernel}")
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
