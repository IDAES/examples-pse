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
