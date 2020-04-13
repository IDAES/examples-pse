# IDAES Examples

See the main README.md for user instructions.

## Developer information

### Rebuild the documentation

To rebuild the docs, the easiest and most useful thing to do is  
run in "_dev_" mode with no options:

    ./build.py dev

This will re-run any notebooks that have changed since the last
time, and generate versions of the notebooks needed for the
documentation.

If you really don't want to re-run the notebooks, even ones that
have changed since you last did this, then add the _--no-notebooks_
option:

    # only build docs
    ./build.py dev --no-notebooks

### Just test notebooks

If you want to test an individual notebook, you can of course
just open it and run it in the Jupyter Notebook environment.
You can also run the `nbconvert` tool from the command-line, which
has an option to execute the notebook and "convert" it to a
Jupyter notebook (which you can throw away or keep, as you wish).

To test all the notebooks at once, you can use the "_test_" mode of the
build script:

    ./build.py test

This will execute all the notebooks, but not modify anything under the
`docs` directory.

### Build for a release

If you are releasing the code, first you must be on a branch
that starts with the word "release". This is to avoid accidentally
doing release actions in a local developer checkout. The release
build will, in addition to running notebooks and generating
documentation, remove any "test" cells in the notebooks in the `src`
directory, adding a notebook with a suffix "_test.ipynb" for every
notebook that had any test cells in it. The reason for doing this is
to make the notebook with the shorter name test-free for the user.
The release mode may also do some miscellaneous cleanup. Like the other  
build commands, it will create a bunch of generated docs.  
To run the release-mode:

    ./build.py release

When actually doing a release, you should rebuild the "release"
branch of the _examples-dev_ repository, tag and push it, and then
use this branch to populate the master branch of the _examples-pse_ 
repository.