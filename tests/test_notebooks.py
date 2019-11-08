"""
Test notebooks
"""
import pytest
import os
from pathlib import Path
from subprocess import Popen, PIPE

# assume parent of this test's directory is top
TOP_DIR = Path(__file__).parent.parent


class Runner:
    def __init__(self, command=None, args=None):
        self.command = [command] + args
        self.stderr = None

    def run(self, srcdir, match=None):
        command = self.command.copy()
        if match:
            command.extend(["--match", match])
        top_dir = TOP_DIR / "src" / Path(srcdir)
        command.append(top_dir)
        command.append("/")

        proc = Popen(command)
        proc.wait()
        return proc.returncode


@pytest.fixture()
def runner():
    kernel = os.environ.get("CONDA_DEFAULT_ENV", "python")
    cmd_path = Path(__file__)
    rnr = Runner(
        command=TOP_DIR / "docs" / "build_notebooks.py",
        args=["--test", "-v", f"--kernel={kernel}"],
    )
    return rnr


def test_workshops(runner):
    assert 0 == runner.run("workshops", match="Solution")


def test_surrogate(runner):
    assert 0 == runner.run("surrogate")


def test_dmf(runner):
    assert 0 == runner.run(Path("properties") / "Workshop_Module_2_DMF")
