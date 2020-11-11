# IDAES Examples

This information is for people creating or modifying the 
content of this repo -- i.e. adding or modifying examples.

Please see the file [README.md](README.md), in this directory, if you are looking for instructions on running
and using the examples.

## Adding an example notebook

To add a new example, follow these instructions

1. Put the notebook in an existing directory, or create a new one, under `src`. Use underscores in
filenames instead of spaces or dashes, e.g. "my_new_notebook.ipynb" and NOT "my-new-notebook.ipynb"
or "My New Notebook.ipynb". 

2. If your notebook will be using tagged cells, please add (before the .ipynb extension)
these suffixes:

   - for notebooks with "testing" cells, add "_testing"
   - for notebooks with "exercise" / "solution" (and potentially testing) cells, use "_solution_testing"

3. Update the `build.yml` file, make sure your notebook will be found under one
 of the notebook directories. Take note of the corresponding output directory in `docs`

4. Add a corresponding entries in the `docs/<dir>/index.rst`, where "<dir>" is the
output directory from `build.yml`. If this is a new directory, add the directory to the table of contents 
in the `docs/index.rst` top-level index, and add some text there to describe the
contents of the directory.

5. Test your new notebook by running `python build.py -cd` (see Build Script below) and seeing that it
gets run and converted. *The first run may take a while, as it needs to execute
the notebooks to create their output. After that, only modified notebooks will be re-run.*

6. Open the HTML documentation generated under
`docs/_build/html` in yur browser and make sure the notebook is there and
linked properly. You can also run the test suite with `pytest`, though this
should usually be unaffected.

7. If you have any Python scripts included, write tests for them and put that
under the `tests` directory.

8. Add the notebook, supporting scripts, tests, and images in `src`, and the `index.rst` files
that you added/changed in `docs`, to Git. Commit and push your result into a 
pull request to this repo. Make sure tests pass on CircleCI.

## Build Script

The documentation is built by running the `build.py` script in this directory.
This script uses the `build.yml` file as its configuration file. You can specify an
alternate configuration file for testing, etc.

Run `python build.py --help` to see basic help information and `python build.py --usage`
to see a more detailed explanation of usage.

If you add new directories with Jupyter notebooks, you will need to let the
script know about them by adding the directories (and the desired target directory
for the rendered documentation) to the `build.yml` file in this directory.
The comments in that file should explain how to add the new information.

You can run the `build.py` script in testing mode (see `-h` option for details) in order
to test the notebooks.

A more limited set of notebooks to examine is configured in the
`build-ci.yml` file, which you can pass to the `--config` option of the build
script. This will emulate how the code is tested on CircleCI during a pull request.

