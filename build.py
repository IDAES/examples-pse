#! /usr/bin/env python
"""
Build Jupyter notebooks and Sphinx docs.

This encapsulates generating Sphinx-ready versions of the Jupyter Notebooks and
calling 'sphinx-build' to build the HTML docs from source. You can jump to the
bottom of this text to see some sample command-lines.

This program uses a configuration file, by default "build.yml" in the same
directory, to set some options and to tell it which directories to process.

The code is structured into "Builder"s, which take the Settings object in
their constructor and actually do something when you call "build()". The
default values for settings are defined in the class.

Most of the code is devoted to running and exporting the Jupyter notebooks,
in the NotebookBuilder class.

There are some quirks to the NotebookBuilder that require further explanation.
First, you will notice some postprocess functions for both the RST and HTML output
files. These deal with some problems with the nbconvert output.

For RST, the image files are all numbered "image0", which causes Sphinx to complain
that they are duplicate references. The code walks through the RST file and numbers the
images sequentially.

For HTML, a whole different kind of problem with the images is that
the Sphinx builder will copy all the images from each subdirectory into a single
"_images" directory in the HTML build area. This will break the image references,
which are simple filenames for the current directory. The Sphinx build step now
copies anything that looks like an image file from the source directory to the
HTML build directory.

Finally, you can see that there is a "wrapper" file generated in ReStructuredText,
i.e. for Sphinx, that uses a template file in `docs/jupyter_notebook_sphinx.tpl`.
This wrapper file makes the notebooks usable as pages in the generated docs by
giving them a title and then including, depending on the documentation mode,
either the raw HTML of the exported notebook (HTML pages), or the RST version (for
LaTeX, PDF, text). The wrapper file will have the same name as the notebook with
"_doc" at the end. For example, if your Jupyter Notebook was called "foo.ipynb",
the whole procedure would generate these files:

  * foo.html    - Notebook as HTML (what you want to see on a web page)
  * foo.rst     - Notebook as RST (other formats)
  * foo_doc.rst - Wrapper doc

Any images exported by running the notebook would also be in the directory.
At any rate, this means that references to "the notebook" in the Sphinx documentation
should be to the "_doc" version of the notebook.

The SphinxBuilder pretty much just runs Sphinx in HTML and PDF modes,
dumping errors to, by default, "sphinx-errors.txt" as well as logging them.
As mentioned above, it copies images and also the raw notebooks (.ipynb) themselves
into the appropriate HTML directory. This will make the notebooks render
properly and also makes it easy to see the "source" of the notebook.

For example command lines see the README.md in this directory.
"""
from abc import ABC, abstractmethod
import argparse
import logging
import os
from pathlib import Path
import shutil
import subprocess
import re
from string import Template
import sys
import tempfile
import time
from typing import TextIO, List
import yaml

#
import nbconvert
from nbconvert.exporters import HTMLExporter, RSTExporter, NotebookExporter
from nbconvert.writers import FilesWriter
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
import nbformat
from traitlets.config import Config

_log = logging.getLogger("build_notebooks")
_hnd = logging.StreamHandler()
_hnd.setLevel(logging.NOTSET)
_hnd.setFormatter(logging.Formatter("%(asctime)s %(levelname)s - %(message)s"))
_log.addHandler(_hnd)
_log.propagate = False
_log.setLevel(logging.INFO)


class NotebookError(Exception):
    pass


class NotebookExecError(NotebookError):
    pass


class NotebookFormatError(NotebookError):
    pass


class SphinxError(Exception):
    pass


class SphinxCommandError(Exception):
    def __init__(self, cmdline, errmsg, details):
        msg = (
            f"Sphinx error while running '{cmdline}': {errmsg}. " f"Details: {details}"
        )
        super().__init__(msg)


# Images to copy to the build directory
IMAGE_SUFFIXES = ".jpg", ".jpeg", ".png", ".gif", ".svg", ".pdf"
# These should catch all data files in the same directory as the notebook, needed for execution.
DATA_SUFFIXES = ".csv", ".json", ".svg", ".xls", ".xlsx", ".txt", ".zip", ".pdf"
CODE_SUFFIXES = (".py",)
NOTEBOOK_SUFFIX = ".ipynb"
TERM_WIDTH = 60  # for notification message underlines


class Settings:
    class ConfigError(Exception):
        def __init__(self, f, msg):
            super().__init__(f"in file '{f.name}': {msg}")

    DEFAULTS = {
        "paths": {"source": "src", "output": "docs", "html": "_build/html"},
        "sphinx": {
            "args": "-b html -T docs {output}/{html}",
            "error_file": "sphinx-errors.txt",
        },
        "notebook": {
            "error_file": "notebook-errors.txt",
            "continue_on_error": True,
            "template": "{output}/jupyter_notebook_sphinx.tpl",
            "directories": [],
            "kernel": None,
            "test_mode": False,
        },
    }

    _NULL = "__null__"

    def __init__(self, f: TextIO):
        self._dsec = None
        try:
            self.d = yaml.safe_load(f)
        except yaml.error.YAMLError as err:
            raise Settings.ConfigError(f, f"YAML format error: {err}")
        missing = [k for k in self.DEFAULTS.keys() if k not in self.d]
        if missing:
            raise Settings.ConfigError(f, f"missing sections: {', '.join(missing)}.")

    def __str__(self):
        return str(self.d)

    def set_default_section(self, section: str):
        """Set default section for get/set.
        """
        self._dsec = section

    def get(self, key: str, default=_NULL):
        """Get a value by key.

        Args:
            key: Key for argument. If it contains ".", taken as "section.name"
                 otherwise use current default section (otherwise ValueError)
            default: Default value. A little different than dict.get(), since
                 providing no value will cause a missing key to raise an error,
                 you have to explicitly pass default=<value> to get the default
                 value to be used.

        Raises:
            ValueError: if the key is malformed
            KeyError: if there is no such value in the section, and no default is provided.

        Returns:
            setting value (dict, list, or str)
        """
        section, val = self._split_key(key)
        values = self.d.get(section, self.DEFAULTS[section])
        if val not in values:
            if default is self._NULL:
                raise KeyError(f"get() failed for value '{val}' in section '{section}'")
            result = self._subst_paths(default)
        else:
            result = self._subst_paths(values[val])
        return result

    def set(self, key: str, value):
        section, val = self._split_key(key)
        self.d[section][val] = value

    def _split_key(self, key):
        parts = key.split(".", 1)
        if len(parts) == 1:
            if self._dsec is None:
                raise ValueError(f"expected <section>.<val>, got '{key}'")
            section, val = self._dsec, parts[0]
        else:
            section, val = parts
        if section not in self.DEFAULTS:
            raise KeyError(
                f"get() failed for section '{section}',"
                f" not in: {', '.join(self.DEFAULTS)}"
            )
        return section, val

    def _subst_paths(self, value):
        if hasattr(value, "format"):
            return value.format(**self.d["paths"])
        return value


class Builder(ABC):
    """Abstract base class for notebook and sphinx builders.
    """

    _rpath = None

    def __init__(self, settings):
        self.s = settings

    @abstractmethod
    def build(self, options):
        pass

    def _merge_options(self, options):
        for key, value in options.items():
            self.s.set(key, value)

    @classmethod
    def root_path(cls):
        if cls._rpath is None:
            cls._rpath = Path(__file__).parent.absolute()
        return cls._rpath


class NotebookBuilder(Builder):
    """Run Jupyter notebooks and render them for viewing.
    """

    # Map format to file extension
    FORMATS = {"html": ".html", "rst": ".rst"}

    HTML_IMAGE_DIR = "_images"
    TEST_SUFFIXES = ("_test", "_testing")  # test notebook suffixes
    # Mapping from {meaning: tag_name}
    _cell_tags = {
        "remove": "remove_cell",
        "exercise": "exercise",
        "testing": "testing",
        "solution": "solution",
    }
    JUPYTER_NB_VERSION = 4  # for parsing

    #: CSS stylesheet for notebooks
    STYLESHEET = """
    #notebook-container {
      box-shadow: none;  /* get rid of cheesy shadow */
      -webkit-box-shadow: none;
      padding 5px;
      width: auto !important; /* expand to container */
      min-width: 400px !important;
    }
    .wy-nav-content {
      max-width: 90%; /* top container for notebook */
    }
    .run_this_cell {
        padding-left: 0px;
        padding-right: 0px;
    }
    /* Hide section name (at top) */
    div.section > h1 {
        display: none;
    }
    """

    class Results:
        """Stores results from build().
        """

        def __init__(self):
            self.failed, self.cached = [], []
            self.dirs_processed = []
            self.duration = -1.0
            self.n_fail, self.n_success = 0, 0
            self._t = None

        def start(self):
            self._t = time.time()

        def stop(self):
            self.duration = time.time() - self._t

        def failures(self) -> str:
            return ", ".join(self.failed)

    def __init__(self, *args):
        super().__init__(*args)
        self._ep = None  # nbconvert ExecutePreprocessor instance (re-used)
        self._imgdir = None  # cached location of images in HTML tree
        self._nb_error_file = None  # error file, for notebook execution failures
        self._results = None  # record results here
        self._nb_remove_config = None  # traitlets.Config for removing test cells

    def build(self, options):
        self.s.set_default_section("notebook")
        self._merge_options(options)
        self._test_mode = self.s.get("test_mode")
        self._open_error_file()
        self._ep = self._create_preprocessor()
        self._imgdir = (
            self.root_path()
            / self.s.get("paths.output")
            / self.s.get("paths.html")
            / self.HTML_IMAGE_DIR
        )
        if _log.isEnabledFor(logging.DEBUG):
            _log.debug(f"settings for build: {self.s}")
        self._read_template()
        # Set up configuration for removing specially tagged cells
        c = Config()
        # Configure for exporters
        c.NotebookExporter.preprocessors = [
            "nbconvert.preprocessors.TagRemovePreprocessor"
        ]  # for some reason, this only works with the full module path
        self._nb_remove_config = c
        nb_dirs = self.s.get("directories")
        self._results = self.Results()
        self._results.start()
        for item in nb_dirs:
            self._build_tree(item)
        self._results.stop()
        return self._results

    def report(self) -> (int, int):
        """Print some messages to the user.

        Returns:
            (total, num_failed) number of notebooks processed and failed.
        """
        r = self._results  # alias
        total = r.n_fail + r.n_success
        _dirs = "y" if len(r.dirs_processed) == 1 else "ies"
        notify(
            f"Processed {len(r.dirs_processed)} director{_dirs} in {r.duration:.3f}s",
            level=1,
        )
        if total == 0:
            notify(f"No notebooks converted", level=2)
        else:
            if len(r.cached) > 0:
                s_verb = " was" if len(r.cached) == 1 else "s were"
                notify(f"{len(r.cached)} notebook{s_verb} cached", level=2)
            if r.n_success == total:
                s_verb = " was" if total == 1 else "s were"
                notify(f"{total} notebook{s_verb} converted successfully", level=2)
            else:
                notify(
                    f"  Out of {total} notebook(s), {r.n_success} converted and {r.n_fail} failed",
                    level=2,
                )
                notify(
                    f"The following {len(r.failed)} notebooks had failures: ", level=1
                )
                for failure in r.failed:
                    notify(f"-  {failure}", level=2)
                if self._nb_error_file is not sys.stderr:
                    notify(
                        f"Notebook errors are in '{self._nb_error_file.name}'", level=2
                    )
        return total, r.n_fail

    def _open_error_file(self):
        """Open error file from value in settings.
        If the filename is '__stdout__' or '__stderr__', use the corresponding stream.
        If opening the file fails, use stderr.
        """
        efile = self.s.get("error_file")
        if efile == "__stdout__":
            self._nb_error_file = sys.stdout
        elif efile == "__stderr__":
            self.nb_error_file = sys.stderr
        else:
            try:
                self._nb_error_file = open(efile, "w")
            except IOError as err:
                _log.error(
                    f"Could not open error file '{efile}': {err}. "
                    f"Writing to standard error on console"
                )
                self._nb_error_file = sys.stderr

    def _create_preprocessor(self):
        kwargs, kernel = {}, self.s.get("kernel", None)
        if kernel:
            kwargs["kernel_name"] = kernel
        return ExecutePreprocessor(timeout=600, **kwargs)

    def _read_template(self):
        nb_template_path = self.root_path() / self.s.get("template")
        try:
            with nb_template_path.open("r") as f:
                nb_template = Template(f.read())
        except IOError as err:
            raise NotebookError(f"cannot read template file {nb_template_path}: {err}")
        self.s.set("Template", nb_template)

    def _build_tree(self, info: dict):
        """"""
        try:
            source, output = info["source"], info["output"]
            match = info.get("match", None)
        except KeyError:
            raise NotebookError(
                f"notebook directory requires values for 'source', "
                f"'output', but got: {info}"
            )
        sroot, oroot = (
            self.root_path() / self.s.get("paths.source"),
            self.root_path() / self.s.get("paths.output"),
        )
        srcdir, outdir = sroot / source, oroot / output
        verb = "Running" if self._test_mode else "Converting"
        notify(f"{verb} notebooks in '{srcdir}'", level=1)
        self._match_expr = re.compile(match) if match else None
        # build, starting at this directory
        self._build_subtree(srcdir, outdir, depth=1)
        self._results.dirs_processed.append(srcdir)

    def _build_subtree(self, srcdir: Path, outdir: Path, depth: int):
        """Build all notebooks in a given directory, exporting into the given output
        format.
        """
        _log.debug(f"build.begin subtree={srcdir}")
        success, failed = 0, 0
        notebooks_to_convert = []
        data_files = []
        for entry in srcdir.iterdir():
            filename = entry.parts[-1]
            if filename.startswith(".") or filename.startswith("__"):
                _log.debug(f"skip special file '{entry}'")
                continue  # e.g. .ipynb_checkpoints
            if entry.is_dir():
                # build sub-directory (filename is really directory name)
                self._build_subtree(entry, outdir / filename, depth + 1)
            elif entry.suffix == NOTEBOOK_SUFFIX:
                _log.debug(f"found notebook at: {entry}")
                # check against match expression
                if self._match_expr and not self._match_expr.search(filename):
                    _log.debug(f"Skip: {entry} does not match {self._match_expr}")
                else:
                    notebooks_to_convert.append(entry)
            else:  # these may all be true
                if entry.suffix in IMAGE_SUFFIXES and not self._test_mode:
                    os.makedirs(self._imgdir, exist_ok=True)
                    shutil.copy(str(entry), str(self._imgdir))
                    _log.debug(f"copied image {entry} to {self._imgdir}/")
                # notebooks might try to open any of these files, sneaky as they are
                if entry.suffix in DATA_SUFFIXES or entry.suffix in CODE_SUFFIXES or entry.suffix in IMAGE_SUFFIXES:
                    data_files.append(entry)

        # convert notebooks last, allowing for discovery of all the data files
        if notebooks_to_convert:
            tmpdir = Path(tempfile.mkdtemp())
            if data_files:
                _log.debug(
                    f"copy {len(data_files)} data file(s) into temp dir '{tmpdir}'"
                )
                for fp in data_files:
                    # print(f"@@ copy {fp} -> {tmpdir}")
                    shutil.copy(fp, tmpdir)
            for entry in notebooks_to_convert:
                if not self._test_mode:
                    os.makedirs(outdir, exist_ok=True)
                # build, if the output file is missing/stale
                verb = "Running" if self._test_mode else "Converting"
                notify(f"{verb}: {entry.name}", level=2)
                continue_on_err = self.s.get("continue_on_error", None)
                try:
                    self._convert(tmpdir, entry, outdir, depth)
                    success += 1
                except NotebookExecError as err:
                    failed += 1
                    self._write_notebook_error(entry, err)
                    self._write_failed_marker(entry, outdir)
                    if continue_on_err:
                        _log.warning(
                            f"Execution failed: generating partial output for '{entry}'"
                        )
                        continue
                    raise
                except NotebookError as err:
                    failed += 1
                    if continue_on_err:
                        _log.warning(f"failed to convert {entry}: {err}")
                        continue
                    raise  # abort all processing
                # remove failed marker, if there was one from a previous execution
                failed_marker = self._get_failed_marker(entry, outdir)
                if failed_marker:
                    _log.info(f"Remove stale marker of failed execution: {failed_marker}")
                    failed_marker.unlink()
                # quick cleanup of contents of temporary dir
                for f in tmpdir.iterdir():
                    if f.suffix not in DATA_SUFFIXES and f.suffix not in CODE_SUFFIXES:
                        try:  # may fail for things like __pycache__
                            f.unlink()
                        except Exception as err:
                            _log.debug(
                                f"Could not unlink file in temp execution dir: {err}"
                            )

            _log.debug(f"build.end subtree={srcdir} {success}/{success + failed}")
            self._results.n_fail += failed
            self._results.n_success += success
            # clean up any temporary directories
            _log.debug(f"remove temporary directory at '{tmpdir.name}'")
            try:
                shutil.rmtree(str(tmpdir))
            except Exception as err:
                _log.error(f"could not remove temporary directory '{tmpdir}': {err}")

    def _write_notebook_error(self, entry, error):
        filename = str(entry)
        self._nb_error_file.write(f"\n====> File: {filename}\n")
        self._nb_error_file.write(str(error))
        self._nb_error_file.flush()  # in case someone is tailing the file
        self._results.failed.append(filename)

    @staticmethod
    def _write_failed_marker(entry, outdir):
        """Put a marker into the output directory for the failed notebook, so we
        can tell whether we need to bother trying to re-run it later.
        """
        marker = outdir / (entry.stem + ".failed")
        _log.debug(f"write failed marker '{marker}' for entry={entry} outdir={outdir}")
        marker.open("w").write("This file is a marker for avoiding re-runs of failed notebooks that haven't changed")

    @staticmethod
    def _get_failed_marker(entry, outdir) -> Path:
        marker = outdir / (entry.stem + ".failed")
        if marker.exists():
            return marker
        return None

    def _convert(self, tmpdir: Path, entry: Path, outdir: Path, depth: int):
        """Convert a notebook.

        Args:
            tmpdir: Temporary working directory
            entry: notebook to convert
            outdir: output directory for .html and .rst files
            depth: depth below root, for fixing image paths
        """
        test_mode = self.s.get("test_mode")
        # strip special cells.
        if self._has_tagged_cells(entry, set(self._cell_tags.values())):
            _log.debug(f"notebook '{entry.name}' has test cell(s)")
            orig_entry, entry = entry, self._strip_tagged_cells(
                tmpdir, entry, ("remove", "exercise"), "testing")
            notify(f"Stripped tags from: {orig_entry.name}", 3)
        else:
            # copy to temporary directory just to protect from output cruft
            tmp_entry = tmpdir / entry.name
            shutil.copy(entry, tmp_entry)
            orig_entry, entry = entry, tmp_entry

        # convert all tag-stripped versions of the notebook.
        # before running, check if converted result is newer than source file
        if self._already_converted(orig_entry, entry, outdir):
            notify(f"Skip notebook conversion, output is newer, for: {entry.name}", 3)
            self._results.cached.append(entry)
            return
        notify(f"Running notebook: {entry.name}", 3)
        nb = self._parse_and_execute(entry)
        if test_mode:  # don't do conversion in test mode
            return
        notify(f"Exporting notebook '{entry.name}' to directory {outdir}", 3)
        wrt = FilesWriter()
        # export each notebook into multiple target formats
        created_wrapper = False
        for (exp, postprocess_func, pp_args) in (
            (RSTExporter(), self._postprocess_rst, ()),
            (HTMLExporter(), self._postprocess_html, (depth,)),
        ):
            _log.debug(f"export '{orig_entry}' with {exp} to notebook '{entry}'")
            (body, resources) = exp.from_notebook_node(nb)
            body = postprocess_func(body, *pp_args)
            wrt.build_directory = str(outdir)
            wrt.write(body, resources, notebook_name=entry.stem)
            # create a 'wrapper' page
            if not created_wrapper:
                _log.debug(f"create wrapper page for '{entry.name}' in '{outdir}'")
                self._create_notebook_wrapper_page(entry.stem, outdir)
                created_wrapper = True
            # move notebooks into docs directory
            _log.debug(f"move notebook '{entry} to output directory: {outdir}")
            shutil.copy(entry, outdir / entry.name)

    def _has_tagged_cells(self, entry: Path, tags: set) -> bool:
        """Quickly check whether this notebook has any cells with the given tag(s).

        Returns:
            True = yes, it does; False = no specially tagged cells
        Raises:
            NotebookFormatError, if notebook at 'entry' is not parseable
        """
        # parse the notebook (assuming this is fast; otherwise should cache it)
        try:
            nb = nbformat.read(str(entry), as_version=self.JUPYTER_NB_VERSION)
        except nbformat.reader.NotJSONError:
            raise NotebookFormatError(f"'{entry}' is not JSON")
        # look for tagged cells; return immediately if one is found
        for i, c in enumerate(nb.cells):
            if "tags" in c.metadata:
                for tag in tags:
                    if tag in c.metadata.tags:
                        _log.debug(f"Found tag '{tag}' in cell {i}")
                        return True  # can stop now, one is enough
        # no tagged cells
        return False

    def _strip_tagged_cells(self, tmpdir, entry, tags, remove_name: str):
        """Strip specially tagged cells from a notebook.

        Copy notebook to a temporary location, and generate a stripped
        version there. No files are modified in the original directory.

        Args:
            tmpdir: directory to copy notebook into
            entry: original notebook
            tags: List of tags (strings) to strip
            remove_name: Remove this component from the name

        Returns:
            stripped-entry, original-entry - both in the temporary directory
        """
        _log.debug(f"run notebook in temporary directory: {tmpdir}")
        # Copy notebook to temporary directory
        tmp_entry = tmpdir / entry.name
        shutil.copy(entry, tmp_entry)
        # Remove the given tags
        # Configure tag removal
        tag_names = [self._cell_tags[t] for t in tags]
        self._nb_remove_config.TagRemovePreprocessor.remove_cell_tags = tag_names
        _log.debug(
            f"removing tag(s) <{', '.join(tag_names)}'> from notebook: {entry.name}"
        )
        (body, resources) = NotebookExporter(
            config=self._nb_remove_config
        ).from_filename(str(tmp_entry))
        # Determine output notebook name:
        # remove suffixes, either "Capitalized" or "lowercase" version
        nmc, nml = remove_name.capitalize(), remove_name.lower()
        nb_name = re.sub(f'(_{nmc}|_{nml}|_{nmc}_|_{nml}_|{nmc}_|{nml}_)', '', entry.stem)
        # Create the new notebook
        wrt = nbconvert.writers.FilesWriter()
        wrt.build_directory = str(tmpdir)
        _log.debug(f"writing stripped notebook: {nb_name}")
        wrt.write(body, resources, notebook_name=nb_name)
        # Return both notebook names, and temporary directory (for cleanup)
        stripped_entry = tmpdir / f"{nb_name}.ipynb"
        return stripped_entry

    def _parse_and_execute(self, entry):

        # parse
        _log.debug(f"parsing '{entry}'")
        try:
            nb = nbformat.read(str(entry), as_version=self.JUPYTER_NB_VERSION)
        except nbformat.reader.NotJSONError:
            raise NotebookFormatError(f"'{entry}' is not JSON")
        except AttributeError:
            raise NotebookFormatError(f"'{entry}' has invalid format")

        # execute
        _log.debug(f"executing '{entry}'")
        try:
            self._ep.preprocess(nb, {"metadata": {"path": str(entry.parent)}})
        except (CellExecutionError, NameError) as err:
            raise NotebookExecError(f"execution error for '{entry}': {err}")
        return nb

    def _create_notebook_wrapper_page(self, nb_file: str, output_dir: Path):
        """Generate a Sphinx documentation page for the Module.
        """
        # interpret some characters in filename differently for title
        title = nb_file.replace("_", " ").title()
        title_under = "=" * len(title)
        # create document from template
        doc = self.s.get("Template").substitute(
            title=title, notebook_name=nb_file, title_underline=title_under
        )
        # write out the new doc
        doc_rst = output_dir / (nb_file + "_doc.rst")
        with doc_rst.open("w") as f:
            _log.info(f"generate Sphinx doc wrapper for {nb_file} => {doc_rst}")
            f.write(doc)

    def _already_converted(self, orig: Path, dest: Path, outdir: Path) -> bool:
        """Check if a any of the output files are either missing or older than
        the input file ('entry').

        Returns:
            True if the output is newer than the input, otherwise False.
        """
        source_time = orig.stat().st_mtime

        # First see if there is a 'failed' marker. If so,
        # compare timestamps: if marker is newer than file, it's converted
        failed_file = outdir / (orig.stem + ".failed")
        if failed_file.exists():
            failed_time = failed_file.stat().st_ctime
            ac = failed_time > source_time
            if ac:
                notify(f"Notebook '{orig.stem}.ipynb' unchanged since previous failed conversion", 3)
            return ac

        # Otherwise look at all the output files and see if any one of them is
        # older than the source file (in which case it's NOT converted)
        for fmt, ext in self.FORMATS.items():
            output_file = outdir / f"{dest.stem}{ext}"
            _log.debug(f"checking if cached: {output_file} src={orig}")
            if not output_file.exists():
                return False
            if source_time >= output_file.stat().st_mtime:
                return False

        return True

    def _postprocess_rst(self, body):
        body = self._replace_image_refs(body)
        return body

    # support up to 35 different images
    IMAGE_NUMBERS = " 123456789abcdefghijklmnopqrstuvwxyz"

    def _replace_image_refs(self, body):
        """Replace |image0| references with successive numbers, instead of
            having them all be "image0". The order is "|image0|" then
            ".. |image0|".
        """
        chars = list(body)  # easy to manipulate this way
        pos, n = 0, 0
        while True:
            # find next image reference
            pos = body.find("|image0|\n", pos)
            if pos == -1:
                break  # no more references; stop
            pos2 = body.find(".. |image0|", pos + 8)
            if pos2 == -1:
                raise ValueError(f"Couldn't find image matching ref at {pos}")
            if n > 0:  # don't do anything for the first one
                c = self.IMAGE_NUMBERS[n]
                # replace '0' with another thing in ref
                chars[pos + 6] = c
                # replace '0' with same other thing in image
                chars[pos2 + 9] = c
            pos = pos2 + 10  # skip past the one we just completed
            n += 1  # increment image number
        return "".join(chars)

    def _postprocess_html(self, body, depth):
        """Change path on image refs to point into HTML build dir.
        """
        # create prefix for <img src='..'> attribute values, which is a relative path
        # to the (single) images directory, from within the HTML build tree
        prefix = Path("")
        for i in range(depth):
            prefix = prefix / ".."
        prefix = prefix / self.HTML_IMAGE_DIR
        # split up content around <img> tags
        splits = re.split(r'<img src="(.*\.png"[^>]*)>', body)
        # replace grouping in odd-numbered splits with modified <img>
        for i in range(1, len(splits), 2):
            orig = splits[i]
            splits[i] = f'<img src="{prefix}/{splits[i]}>'
            # print(f"@@ depth={depth}; rewrote image link '{orig}' as '{splits[i]}'")
        # rejoin splits, to make the modified body text
        body = "".join(splits)
        # hack in some CSS, replacing useless link
        custom = re.search(r'<\s*link.*href="custom.css"\s*>', body)
        if custom:
            p1, p2 = custom.span()
            body = "".join(
                (
                    body[:p1],
                    '<style type="text/css">',
                    self.STYLESHEET,
                    "</style>",
                    body[p2:],
                )
            )
        else:
            _log.warning("Could not insert stylesheet in HTML")
        # done
        return body


#


class SphinxBuilder(Builder):
    """Run Sphinx documentation build command.
    """

    def build(self, options):
        self.s.set_default_section("sphinx")
        self._merge_options(options)
        raw_args = self.s.get("args")
        _log.debug(f"Sphinx args from settings: {raw_args}")
        args = raw_args.split()
        # prepend root path to last source and output directories
        args[-2] = str(self.root_path() / args[-2])
        args[-1] = str(self.root_path() / args[-1])
        html_dir = (
            self.root_path() / self.s.get("paths.output") / self.s.get("paths.html")
        )
        doc_dir = self.root_path() / self.s.get("paths.output")
        src_dir = self.root_path() / self.s.get("paths.source")

        # Catch some problems that may cause an ugly stack trace from Sphinx
        if not doc_dir.exists():
            raise SphinxError(f"Input directory does not exist: {doc_dir}")
        if not (doc_dir / "conf.py").exists():
            raise SphinxError(
                f"Input directory does not have 'conf.py' file: {doc_dir}"
            )

        if not os.path.exists(html_dir):
            _log.warning(f"Target HTML directory {html_dir} does not exist: creating")
            html_dir.mkdir(parents=True)

        # Run Sphinx command
        errfile = self.s.get("error_file")
        cmdargs = ["sphinx-build", "-a", "-w", errfile] + args
        cmdline = " ".join(cmdargs)
        notify(f"Running Sphinx command: {cmdline}", level=1)
        proc = subprocess.Popen(cmdargs)
        proc.wait()
        status = proc.returncode
        if status != 0:
            log_error = self._extract_sphinx_error(errfile)
            raise SphinxCommandError(cmdline, f"return code = {status}", log_error)

        # copy notebooks from doc & src directories into html directory
        notify(f"Copying notebooks from '{doc_dir}' -> '{html_dir}'", 1)
        for nb_dir in self.s.get("notebook.directories"):
            nb_output_dir = doc_dir / nb_dir["output"]
            nb_source_dir = src_dir / nb_dir["source"]
            _log.debug(f"find notebooks in path: {nb_output_dir}")
            for nb_path in Path(nb_output_dir).glob("**/*.ipynb"):
                nb_dest = html_dir / nb_path.relative_to(doc_dir)
                notify(f"Copy notebook {nb_path.name} -> {nb_dest.parent}", 2)
                shutil.copy(nb_path, nb_dest)
            _log.info(f"find supporting files in path: {nb_output_dir}")
            for ext in IMAGE_SUFFIXES:
                pattern = f"**/*{ext}"
                for nb_path in Path(nb_source_dir).glob(pattern):
                    nb_dest = html_dir / nb_path.relative_to(src_dir)
                    notify(f"Copy supporting file {nb_path.name} -> {nb_dest.parent}", 3)
                    shutil.copy(nb_path, nb_dest)

    @staticmethod
    def _extract_sphinx_error(errfile: str):
        """Extract error message from a Sphinx output file.
        """
        collect_errors, lines = None, []
        with open(errfile) as f:
            for line in f:
                s = line.strip()
                if collect_errors:
                    lines.append(s)
                elif s.startswith("Sphinx error:"):
                    collect_errors = True
        return "\n".join(lines)


class Cleaner(Builder):
    """Clean up all the pre-generated files, mostly in the docs.
    """

    def build(self, options):
        root = self.root_path()
        docs_path = root / self.s.get("paths.output")
        # clean in each path
        did_remove = False
        did_remove |= self._clean_html(docs_path / self.s.get("paths.html"))
        did_remove |= self._clean_docs(docs_path)
        if not did_remove:
            notify("No files removed", 1)

    @staticmethod
    def _clean_html(html_path):
        if html_path.exists():
            notify(f"remove directory: {html_path}", 1)
            shutil.rmtree(html_path)
            return True
        return False

    def _clean_docs(self, docs_path):
        removed_any = False
        stop_list = {"index.rst"}
        file_types = ["rst", "html", "ipynb", "failed"] + list(IMAGE_SUFFIXES)
        for notebook_dirs in self.s.get("notebook.directories"):
            docs_output = docs_path / notebook_dirs["output"]
            for file_type in file_types:
                for f in docs_output.glob(f"**/*.{file_type}"):
                    if f.name not in stop_list:
                        notify(f"remove: {f}", 1)
                        f.unlink()
                        removed_any = True
        return removed_any


class Color:
    BLACK = "\033[30m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"
    UNDERLINE = "\033[4m"
    RESET = "\033[0m"


def notify(message, level=0):
    """Multicolored, indented, messages to the user.
    """
    c = [Color.MAGENTA, Color.GREEN, Color.CYAN][min(level, 2)]
    indent = "  " * level
    if level == 0:
        print()
    print(f"{indent}{c}{message}{Color.RESET}")


def get_git_branch():
    _log.debug("Run command: git branch --list")
    proc = subprocess.run(
        ["git", "branch", "--list"], encoding="utf-8", capture_output=True
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr)
    branch = None
    for line in proc.stdout.split("\n"):
        if len(line) > 1 and line[0] == "*":  # current branch
            branch = line[1:].strip()
    if branch is None:
        raise RuntimeError("cannot find current branch in output")
    return branch


def print_usage():
    """Print a detailed usage message.
    """
    command = "python build.py"
    message = (
        "\n"
        "# tl;dr To convert notebooks and build docs, use this command:\n"
        "{command} -crd\n"
        "\n"
        "The build.py command is used to create the documentation from\n"
        "the Jupyter Notebooks and hand-written '.rst' files in this\n"
        "repository. It is also used by tests, programmatically, to run\n"
        "the notebooks. Some sample command-lines, with comments as to \n"
        "what they will do, follow. The operation of this script is also\n"
        "controlled by the configuration file, 'build.yml' by default.\n"
        "\n"
        "# Read the default configuration file, but do nothing else\n"
        "{command}\n"
        "\n"
        "# Build the Sphinx documentation, only\n"
        "{command} --docs\n"
        "{command} -d  # <-- short option\n"
        "\n"
        "# Remove all built documentation files\n"
        "{command} --remove\n"
        "{command} -r  # <-- short option\n"
        "\n"
        "# Convert Jupyter notebooks. Only those notebooks\n"
        "# that have not changed since the last time this was run will\n"
        "# be re-executed. Converted notebooks are stored in the 'docs'\n"
        "# directory, in locations configured in the 'build.yml'\n"
        "# configuration file.\n"
        "{command} --convert\n"
        "{command} -c  # <-- short option\n"
        "\n"
        "# Convert Jupyter notebooks, as in previous command,\n"
        "# then build Sphinx documentation.\n"
        "# This can be combined with -r/--remove to convert all notebooks.\n"
        "{command} -cd\n"
        "\n"
        "# Run notebooks, but do not convert them into any other form.\n"
        "# This can be combined with -r/--remove to run all notebooks.\n"
        "{command} --test\n"
        "{command} -t  # <-- short option\n"
        "\n"
        "# Run with <options> at different levels of verbosity\n"
        "{command} <options>      # Show warning, error, fatal messages\n"
        "{command} -v <options>   # Add informational (info) messages\n"
        "{command} -vv <options>  # Add debug messages\n"
        "\n"
    )
    print(message.format(command=command))


def main():
    ap = argparse.ArgumentParser(
        description="Build documentation and/or Jupyter notebooks"
    )
    ap.add_argument(
        "--config",
        "-C",
        default="build.yml",
        metavar="FILE",
        help="Location of configuration file (default=./build.yml)",
    )
    ap.add_argument(
        "--usage", "-U", action="store_true", help="Print a more detailed usage message"
    )
    ap.add_argument("--remove", "-r", action="store_true", help="Remove generated files")
    ap.add_argument("--docs", "-d", action="store_true", help="Build documentation")
    ap.add_argument(
        "--convert", "-c", action="store_true", help="Convert Jupyter notebooks",
    )
    ap.add_argument(
        "--test",
        "-t",
        dest="test_mode",
        action="store_true",
        help="Run notebooks but do not convert them.",
    )
    ap.add_argument(
        "-v",
        "--verbose",
        action="count",
        dest="vb",
        help="Increase verbosity",
        default=0,
    )
    args = ap.parse_args()

    # Check for usage message flag
    if args.usage:
        print_usage()
        return 0

    # Check for confusing option combinations
    if args.convert and args.test_mode:
        ap.error(
            "-t/--test conflicts with notebook conversion -c/--convert; pick one"
        )
    if args.docs and args.test_mode:
        ap.error(
            "-t/--test should not be used with -d/--docs, as it does not convert any notebooks"
        )

    # If building docs, check for working Sphinx command
    if args.docs:
        build_command = "sphinx-build"
        try:
            subprocess.check_call([build_command, "--version"])
        except (FileNotFoundError, subprocess.CalledProcessError) as err:
            ap.error(f"Cannot run doc build command: {err}")

    # Set verbosity
    if args.vb > 1:
        _log.setLevel(logging.DEBUG)
    elif args.vb > 0:
        _log.setLevel(logging.INFO)
    else:
        _log.setLevel(logging.WARNING)

    # Read and parse configuration file
    _log.debug(f"reading settings file '{args.config}'")
    try:
        conf_file = open(args.config, "r")
        settings = Settings(conf_file)
    except IOError as err:
        _log.fatal(f"Cannot open settings file '{args.config}': {err}")
        return 1
    except Settings.ConfigError as err:
        _log.fatal(f"Cannot read settings from '{args.config}': {err}")
        return 1

    # set local variables from command-line arguments
    run_notebooks = args.convert
    build_docs = args.docs
    clean_files = args.remove
    test_mode = args.test_mode

    # Clean first, if requested
    if clean_files:
        notify("Clean all built files")
        cleaner = Cleaner(settings)
        cleaner.build({})

    status_code = 0  # start with success

    nbb = None
    if run_notebooks or test_mode:
        verb = "Run" if test_mode else "Convert"
        notify(f"{verb} Jupyter notebooks")
        nbb = NotebookBuilder(settings)
        try:
            nbb.build(
                {"test_mode": test_mode}
            )
        except NotebookError as err:
            _log.fatal(f"Could not build notebooks: {err}")
            return -1
        run_total, run_failed = nbb.report()
        if run_failed > 0:
            status_code = 1

    if build_docs:
        notify("Build documentation with Sphinx")
        spb = SphinxBuilder(settings)
        try:
            spb.build({})
        except SphinxError as err:
            _log.fatal(f"Could not build Sphinx docs: {err}")
            return -1

    return status_code


if __name__ == "__main__":
    sys.exit(main())
