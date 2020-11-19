# IDAES Examples

This information is for people creating or modifying the 
content of this repo -- i.e. adding or modifying examples.

Please see the file [README.md](README.md), in this directory, if you are looking for instructions on running
and using the examples.

## Adding an example notebook

To add a new example, follow these instructions

1. Decide where the notebook will go. The first level of folders has two choices: Examples or Tutorials.
If your notebook as "exercise" cells and "solution" cells in it, then it is a tutorial and it belongs under
Tutorials. Otherwise, it is an example and belongs under Examples. The next levels of folders determine the
category and sub-category to which the example belongs. Please read the "notebook_index.yml" file to see
the full descriptions of each category (given to the user) and, if in doubt, ask. In rare cases, you will need
to add a new sub-folder. In this case, use the "camel-case" convention to name it, e.g., "SurrMod" rather
than "Surr Mod" or "surr_mod", and try to stick with short names (please no "AdvancedMethodsToDoSomethingInteresting"
types of names).

2. Name the notebook. The convention is to Use underscores in filenames instead of spaces or dashes, e.g.
"my_new_notebook.ipynb" and NOT "my-new-notebook.ipynb" or "My New Notebook.ipynb". Do not capitalize any 
letters in the notebook name, except acronyms. Don't use the word "example" or "tutorial" in the notebook title.

3. If your notebook uses tagged cells (see below for explanation), please add (before the .ipynb extension)
these suffixes:

   - for notebooks with "testing" cells, add "_testing"
   - for notebooks with "exercise" / "solution" (and potentially testing) cells, use "_solution_testing",
     even if there aren't any testing cells

4. Once you have picked a location and named your notebook, you can add the notebook file and any supporting files
in the appropriate folder under `src/` with that name. For example, if you are adding a notebook named
"pacman_game", with testing cells, under the "Examples/Tools" subfolder, then you would put the notebook under
`src/Examples/Tools/pacman_game_testing.ipynb`.

5. Add your notebook to `notebook_index.yml`, which is in the root directory of this repository. This is a YAML
   file, so be careful to follow the exact format used for the other entries. If you added a subfolder, you will
   need to add a new entry under a "subfolders:" heading, with a name, title, description, and list of notebooks
   that are being added under that subfolder. For each notebook, add an entry with the name of the notebook,
   exactly matching the filename, a colon, and the description. Use existing entries as a template.
   
6. If you added a subfolder, add this to the `build.yml` file, so the build process will see it. Otherwise,
   you should not need to do anything with this file.

7. Add the notebook to the documentation section. This means that you need to add a line in a file called `index.rst`
   in a directory under `docs/` matching the directory for the notebook under `src`. So, for the 
   notebook `src/Examples/Tools/pacman_game_testing.ipynb`, you would add an entry into `docs/Examples/Tools/index.rst`.
   The entry goes under the ".. toctree::" directive in that file. You can imitate the form of the existing
   entries, which is basically: "Description of notebook <name_of_notebook_doc>". For tutorial notebooks, the suffixes
   you added in Step 3 are replaced with "_solution" and for all other notebooks the sufixes are removed.
   Returning to our fake example, you would add an entry like: "Use IDAES to play Pac-Man <pacman_game_doc>".
   
8. Now everything should be ready to test! There are three phases to this.

9. Re-generate the notebook index (src/notebook_index.ipynb) by running:
 
        python build.py -x
        
    This will create a new version of the index. Errors here are usually due to indentation 
    or other formatting issues with the YAML in `notebook_index.yml`.

10. Test the notebook by running:

        python build.py -t
        
    This will run all notebooks that have changed since the last time you tested them. If
    you don't see your notebook running, it's possible there is a problem with the `build.py`
    configuration that is causing it not to be seen. If your notebook runs and fails, you
    should see that clearly in the output at the end.
    
11. Build your notebook as a static HTML page in the docs by running:

        python build.py -cd
    
    The "c" is convert and the "d" is for docs, so this does two things: converts the notebook
    to forms used by the docs (reStructuredText and HTML) and rebuilds the docs with Sphinx.
    If you see Sphinx errors about not finding the notebook, that probably means there is a
    problem with the entry you made in `index.rst`.
    
12. You're done! Check that the HTML looks good by opening the page
  `docs/_build/html/index.html` in your browser and make sure the notebook is there and
  linked properly. Remember to check in to Git the notebook, supporting files, and updated
  `notebook_index.yml` and `src/notebook_index.ipynb`, as well as any supporting files
  you need.

## Notebook tags

Coming soon!

## Build Script

The documentation is built by running the `build.py` script in this directory.
This script uses the `build.yml` file as its configuration file. You can specify an
alternate configuration file for testing, etc.

Run `python build.py --help` to see basic help information and `python build.py --usage`
to see a more detailed explanation of usage.

You can run the `build.py` script in testing mode (see `-h` option for details) in order
to test the notebooks.

A more limited set of notebooks to examine is configured in the
`build-ci.yml` file, which you can pass to the `--config` option of the build
script. This will emulate how the code is tested on CircleCI during a pull request.

