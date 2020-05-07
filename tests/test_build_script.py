"""
Test the build.py script.

This imports and uses that script like a module of course.
"""
import json
import logging
from pathlib import Path
import sys

# third-party
import pytest

# assume we are at ROOT/tests/something.py and the script is at ROOT/build.py
sys.path.insert(0, str(Path(__file__).parent.parent.absolute()))

print(f"path: {sys.path}")
import build


@pytest.fixture
def hello_world_nb():
    return {
        "cells": [
            {
                "cell_type": "code",
                "execution_count": 2,
                "metadata": {},
                "outputs": [
                    {
                        "name": "stdout",
                        "output_type": "stream",
                        "text": ["Hello, world\n"],
                    }
                ],
                "source": ['print("Hello, world")'],
            },
            {
                "cell_type": "code",
                "execution_count": 3,
                "metadata": {"tags": ["test", "remove_cell"]},
                "outputs": [],
                "source": ["assert 2 + 2 == 4"],
            },
        ],
        "metadata": {
            "celltoolbar": "Tags",
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "codemirror_mode": {"name": "ipython", "version": 3},
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "nbconvert_exporter": "python",
                "pygments_lexer": "ipython3",
                "version": "3.8.1",
            },
        },
        "nbformat": 4,
        "nbformat_minor": 2,
    }


@pytest.mark.unit
def test_notify():
    build.notify("hello")
    build.notify("hello", level=1)
    build.notify("hello", level=10)


def _setup_builder(tmp_path, hello_world_nb):
    json.dump(hello_world_nb, (tmp_path / "hello_world.ipynb").open("w"))
    (tmp_path / "settings.yml").open("w").write(json.dumps(build.Settings.DEFAULTS))
    settings = build.Settings((tmp_path / "settings.yml").open())
    nb_builder = build.NotebookBuilder(settings)
    options = {
        "paths.source": str(tmp_path),
        "paths.output": str(tmp_path),
        "paths.html": str(tmp_path),
        "notebook.directories": [{
            "source": ".",
            "output": "."
        }],
    }
    # create dummy template file
    (tmp_path / Path(settings.get("notebook.template")).name).open("w")
    # strip the notebooks
    build._log.setLevel(logging.DEBUG)
    return nb_builder, options


@pytest.mark.component
def test_strip_tagged_cells(tmp_path, hello_world_nb):
    nb_builder, options = _setup_builder(tmp_path, hello_world_nb)
    options["notebook.rebuild"] = True
    nb_builder.build(options)


