# IDAES Examples

This information is for people creating or modifying the 
example notebooks, scripts, etc.
Please see the main README.md for user instructions.

## Building

The documentation is built by running the `build.py` script in this directory.

Run `python build.py --help` to see basic help information and `python build.py --usage`
to see a more detailed explanation of usage.

If you add new directories with Jupyter notebooks, you will need to let the
script know about them by adding the directories (and the desired target directory
for the rendered documentation) to the `build.yml` file in this directory.
The comments in that file should explain how to add the new information.

## Testing

You can run the `build.py` script in testing mode (see `-h` option for details) in order
to test the notebooks.

A more limited set of notebooks to examine is configured in the
`build-circleci.yml` file, which you can pass to the `--config` option of the build
script. This will emulate how the code is tested on CircleCI during a pull request.

## Adding new examples

To add a new example, add a directory (or use an existing one) under `src`. Place all
code, data files, Jupyter Notebooks, etc. in that directory. Also add a corresponding
directory under `docs`, with an `index.rst` file (see other directories for examples
of how to do this) and add this directory to the table of contents in the `index.rst`
file at top level of `docs`. Finally, if this is a new directory, add an entry in the
"directories" section of the `build.yml` file. Follow the pattern of existing entries,
and use documentation 
provided in comments in that file.

Add all these files to Git, as these are all source material. Files generated with
subsequent steps should *not* be added to Git.

To test your new examples, simply run `build.py` with appropriate arguments and look
at the output, including opening the built documentation with a browser to make sure
nothing is missing or malformed. 

Note: The first run may take a while, as it needs to execute
the notebooks to create their output. After that, only modified notebooks will be re-run.