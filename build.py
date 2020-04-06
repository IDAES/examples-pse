#!/usr/bin/env python
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
in the NotebookBuilder class. This can operate in two basic modes, controlled
by the "notebook.test_mode" setting:

    * Documentation mode (the default), will create ReStructuredText and HTML
      versions of the notebooks, after executing them. Execution of the notebooks
      is optional, but needed for full results. By default, the timestamps of the
      generated files and notebook source are compared, and execution is skipped if
      the generated files are newer.

    * Test mode will run the notebooks, but not try and generate any documentation.
      This is useful, obviously, from within Python tests.

There are some quirks to the NotebookBuilder that require further explanation.
First, you will notice some postprocess functions for both the RST and HTML output
files. These deal with some problems with the nbconvert output.

For RST, the image files are all numbered "image0", which causes Sphinx to complain
that they are duplicate references. The code walks through the RST file and numbers the
images sequentially.

For HTML, a whole different kind of problem with the images is that
the Sphinx builder will copy all the images from each subdirectory into a single
"_images" directory in the HTML build area. This will break the image references,
which are simple filenames for the current directory. The preprocessor goes through
and replaces all <img> tags, with ".png" file source, with a path to the
images directory. This has two corollaries: no two images anywhere in the tree
should have the same name, and all images should be PNG format.

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
should be to the "_doc" version of the notebook, e.g. from the `tutorials/index.rst`
table of contents::

    Welcome tutorial (short) <Module_0_Welcome/Module_0_Welcome_Short_Solution_doc>

The SphinxBuilder pretty much just runs Sphinx in HTML mode (PDF coming R.S.N.),
dumping errors to, by default, "sphinx-errors.txt" as well as logging them.
It doesn't do any other accommodations, as all the work of dealing what Sphinx wants
to do by default is handled in the NotebookBuilder.

For example command lines see the README.md in this directory.
"""
from abc import ABC, abstractmethod
import argparse
from dataclasses import dataclass
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

_log = logging.getLogger("build_notebook_html")
_hnd = logging.StreamHandler()
_hnd.setFormatter(logging.Formatter("%(asctime)s %(levelname)s - %(message)s"))
_log.addHandler(_hnd)


class NotebookError(Exception):
    pass


class NotebookExecError(NotebookError):
    pass


class NotebookFormatError(NotebookError):
    pass


class SphinxError(Exception):
    def __init__(self, cmdline, errmsg, details):
        msg = (
            f"Sphinx error while running '{cmdline}': {errmsg}. " f"Details: {details}"
        )
        super().__init__(msg)


IMAGE_SUFFIXES = ".jpg", ".jpeg", ".png", ".gif", ".svg", ".pdf"
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
            "clean": False,
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
            KeyError: if there is no such value in the section, and no default
                      is provided.

        Returns:
            setting value (dict, list, or str)
        """
        section, val = self._split_key(key)
        values = self.d.get(section, self.DEFAULTS[section])
        if val not in values:
            if default is self._NULL:
                raise KeyError(f"get() failed for value '{val}' in section '{section}'")
            return self._subst_paths(default)
        return self._subst_paths(values[val])

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
    TEST_SUFFIX = "_test"  # for notebooks *with* tests
    CLEAN_SUFFIX = "_clean"  # for notebooks without tests
    REMOVE_CELL_TAG = "remove_cell"

    JUPYTER_NB_VERSION = 4  # for parsing

    @dataclass
    class Results:
        """Stores results from build().
        """

        _t = None
        failed: List[str]
        dirs_processed: List[str]
        duration: float = -1.0
        n_fail: int = 0
        n_success: int = 0

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
        # Configure tag removal
        c.TagRemovePreprocessor.remove_cell_tags = (self.REMOVE_CELL_TAG,)
        # Configure for exporters
        c.NotebookExporter.preprocessors = [
            "nbconvert.preprocessors.TagRemovePreprocessor"
        ]  # for some reason, this only works with the full module path
        self._nb_remove_config = c
        nb_dirs = self.s.get("directories")
        self._results = self.Results(failed=[], dirs_processed=[])
        self._results.start()
        for item in nb_dirs:
            self._build_tree(item)
        self._results.stop()
        return self._results

    def report(self):
        """Print some messages to the user.
        """
        r = self._results  # alias
        total = r.n_fail + r.n_success
        notify(
            f"Built {r.n_success}/{total} notebooks in {len(r.dirs_processed)} "
            f"directories in {r.duration:.3f} seconds"
        )
        notify("Notebook build summary", level=1)
        if total == 0:
            notify(f"No notebooks converted", level=2)
        elif r.n_fail == 0:
            notify(f"All {r.total} notebooks executed successfully", level=2)
        else:
            notify(
                f"The following {len(r.failed)} notebooks had failures: "
                f"{r.failures()}",
                level=2,
            )
            if self._nb_error_file is not sys.stderr:
                notify(f"Notebook errors are in '{self._nb_error_file.name}'", level=2)

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
        notify(f"Converting notebooks in '{srcdir}'", level=1)
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
        for entry in srcdir.iterdir():
            filename = entry.parts[-1]
            if filename.startswith(".") or filename.startswith("__"):
                _log.debug(f"skip special file '{entry}'")
                continue  # e.g. .ipynb_checkpoints
            if entry.is_dir():
                # build sub-directory (filename is really directory name)
                self._build_subtree(entry, outdir / filename, depth + 1)
            elif entry.suffix in IMAGE_SUFFIXES:
                shutil.copy(str(entry), str(self._imgdir))
                _log.debug(f"copied image {entry} to {self._imgdir}/")
            elif entry.suffix == NOTEBOOK_SUFFIX:
                _log.debug(f"found notebook at: {entry}")
                # check against match expression
                if self._match_expr and not self._match_expr.search(filename):
                    _log.debug(f"does not match {self._match_expr}, skip {entry}")
                    continue
                self._ensure_target_dir(outdir)
                # build, if 'rebuild'-mode or the output file is missing/stale
                if not self.s.get("rebuild") and self._is_cached(entry, outdir):
                    _log.info(f"skip converting notebook '{entry}': not changed")
                else:
                    notify(f"{entry}", level=2)
                    continue_on_err = self.s.get("continue_on_error", None)
                    try:
                        self._convert(entry, outdir, depth)
                        success += 1
                    except NotebookExecError as err:
                        failed += 1
                        self._notebook_error(entry, err)
                        if continue_on_err:
                            _log.warning(f"generating partial output for '{entry}'")
                            continue
                        raise
                    except NotebookError as err:
                        failed += 1
                        if continue_on_err:
                            _log.warning(f"failed to convert {entry}: {err}")
                            continue
                        raise  # abort all processing
        _log.debug(f"build.end subtree={srcdir} {success}/{success + failed}")
        self._results.n_fail += failed
        self._results.n_success += success

    def _notebook_error(self, entry, error):
        filename = str(entry)
        self._nb_error_file.write(f"\n====> File: {filename}\n")
        self._nb_error_file.write(str(error))
        self._nb_error_file.flush()  # in case someone is tailing the file
        self._results.failed.append(filename)

    def _convert(self, entry: Path, outdir: Path, depth: int):
        """Convert a notebook.

        Args:
            entry: notebook to convert
            outdir: output directory for .html and .rst files
            depth: depth below root, for fixing image paths
        """
        if self.s.get("test_mode"):
            self._parse_and_execute(entry)
            return  # skip export
        if self._has_tagged_cells(entry):
            _log.debug(f"notebook '{entry}' has test cell(s)")
            entries, tmpdir = self._strip_tagged_cells(entry)
        else:
            entries, tmpdir = [entry], None
        # main loop
        try:
            for e in entries:  # one or two notebooks to export
                _log.debug(f"exporting '{e}' to directory {outdir}")
                nb = self._parse_and_execute(e)
                wrt = FilesWriter()
                # export each notebook into multiple target formats
                for (exp, postprocess_func, pp_args) in (
                    (RSTExporter(), self._postprocess_rst, ()),
                    (HTMLExporter(), self._postprocess_html, (depth,)),
                ):
                    _log.debug(f"export '{e}' with {exp}")
                    (body, resources) = exp.from_notebook_node(nb)
                    body = postprocess_func(body, *pp_args)
                    wrt.build_directory = str(outdir)
                    wrt.write(body, resources, notebook_name=e.stem)
            # create a 'wrapper' page for the main (first) entry
            self._create_notebook_wrapper_page(entries[0].stem, outdir)
        finally:
            # clean up any temporary directories
            if tmpdir is not None:
                _log.debug(f"remove temporary directory at '{tmpdir.name}'")
                shutil.rmtree(str(tmpdir))

    def _has_tagged_cells(self, entry: Path) -> bool:
        """Quickly check whether this notebook has any cells with the "special" tag.

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
            if "tags" in c.metadata and self.REMOVE_CELL_TAG in c.metadata.tags:
                _log.debug(f"Found {self.REMOVE_CELL_TAG} tag in cell {i}")
                return True  # can stop now, one is enough
        # no tagged cells
        return False

    def _strip_tagged_cells(self, entry):
        """Split a notebook into one with and without the specially tagged cells.
        """
        if self.s.get("clean"):
            # generate a "clean" notebook in place of original
            # copy 'raw' original to <name>_test.ipynb
            entry_copy = Path(entry.parent) / f"{entry.stem}{self.TEST_SUFFIX}.ipynb"
            _log.debug(f"copy raw notebook '{entry}' -> '{entry_copy}'")
            shutil.copy(entry, entry_copy)
            # replace 'raw' notebook with stripped one
            (body, resources) = NotebookExporter(
                config=self._nb_remove_config  # config removes cells with special tags
            ).from_filename(str(entry))
            wrt = nbconvert.writers.FilesWriter()
            wrt.build_directory = str(entry.parent)
            wrt.write(body, resources, notebook_name=entry.stem)
            _log.debug(f"stripped tags from '{entry}'")
            return [entry, entry_copy], None
        else:
            # copy notebook to a temporary location, and generate a stripped
            # version there. No files are modified in the original directory.
            tmpdir = Path(tempfile.mkdtemp())  # already checked this
            _log.debug(f"run notebook in temporary directory: {tmpdir}")
            raw_entry = tmpdir / entry.name
            shutil.copy(entry, raw_entry)  # copy into temp dir
            (body, resources) = NotebookExporter(
                config=self._nb_remove_config  # config removes cells with special tags
            ).from_filename(str(raw_entry))
            wrt = nbconvert.writers.FilesWriter()
            wrt.build_directory = str(tmpdir)
            cleaned_entry = tmpdir / f"{entry.stem}{self.CLEAN_SUFFIX}.ipynb"
            wrt.write(body, resources, notebook_name=cleaned_entry.stem)
            _log.debug(f"stripped tags from '{raw_entry}' -> '{cleaned_entry}'")
            return [raw_entry, cleaned_entry], tmpdir

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
        title = nb_file.replace("_", " ")
        title_under = "=" * len(title)
        # create document from template
        doc = self.s.get("Template").substitute(
            title=title, notebook_name=nb_file, title_underline=title_under
        )
        # write out the new doc
        doc_rst = output_dir / (nb_file + "_doc.rst")
        with doc_rst.open("w") as f:
            _log.debug(f"generate Sphinx doc wrapper for {nb_file} => {doc_rst}")
            f.write(doc)

    @staticmethod
    def _ensure_target_dir(p: Path):
        if not p.exists():
            _log.debug(f"Creating output directory {p}")
            p.mkdir()

    def _is_cached(self, entry: Path, outdir: Path) -> bool:
        """Check if a any of the output files are either missing or older than
        the input file ('entry').

        Returns:
            True if the output is newer than the input, otherwise False.
        """
        for fmt, ext in self.FORMATS.items():
            output_file = outdir / f"{entry.stem}{ext}"
            if not output_file.exists():
                return False
            if entry.stat().st_mtime >= output_file.stat().st_mtime:
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
        """Change path on `.png` image refs to point into HTML build dir.
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
            splits[i] = f'<img src="{prefix}/{splits[i]}>'
        # rejoin splits, to make the modified body text
        return "".join(splits)


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
        if not os.path.exists(html_dir):
            _log.warning("Target HTML directory {html_dir} does not exist: creating")
            html_dir.mkdir(parents=True)
        errfile = self.s.get("error_file")
        cmdargs = ["sphinx-build", "-a", "-w", errfile] + args
        cmdline = " ".join(cmdargs)
        notify(f"Running Sphinx command: {cmdline}", level=1)
        proc = subprocess.Popen(cmdargs)
        proc.wait()
        status = proc.returncode
        if status != 0:
            log_error = self._extract_sphinx_error(errfile)
            raise SphinxError(cmdline, f"return code = {status}", log_error)

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


def notify(message, level=0):
    if level == 0:
        print("*" * TERM_WIDTH)
        print(f"* {message}")
        print("*" * TERM_WIDTH)
    elif level == 1:
        print(f"===> {message}")
    elif level == 2:
        print(f".... {message}")


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


# Run modes
MODE_DEV, MODE_REL, MODE_TEST = "dev", "release", "test"


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Build documentation and Jupyter notebooks",
        epilog=f"Modes of operation\n"
        f"------------------\n"
        f"* {MODE_DEV:7s}  Build docs for local development\n"
        f"* {MODE_TEST:7s}  Run all notebooks for local testing. "
        f"* {MODE_REL:7s}  Build notebooks and docs for release\n"
        f"Implies '--no-docs' option.\n",
    )
    ap.add_argument(
        "mode",
        help="Mode of operation (see below)",
        choices=[MODE_DEV, MODE_REL, MODE_TEST],
    )
    ap.add_argument(
        "-c",
        "--config",
        default="build.yml",
        help="Location of configuration file (default=./build.yml)",
    )
    ap.add_argument(
        "--no-notebooks", action="store_true", help="Do not run any Jupyter notebooks"
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

    # Need temporary directories in "dev" mode
    if args.mode == MODE_DEV and tempfile.gettempdir() == os.getcwd():
        ap.error(
            "ABORT: temporary directory is going to be in "
            "the current directory. Please set one of the "
            "following environment variables to a suitable "
            "directory: TMPDIR, TEMP, or TMP"
        )

    # Cannot be in master branch for "release" mode
    if args.mode == MODE_REL:
        branch = None
        try:
            branch = get_git_branch()
        except RuntimeError as err:
            ap.error(f"Cannot run 'git' to get branch: {err}")
        if not branch.lower().startswith("release"):
            ap.error(
                f"Current git branch is '{branch}', but in '{MODE_REL}' mode,\n"
                f"the branch name must start with 'release' (case-insensitive)"
            )

    # Set verbosity
    if args.vb > 1:
        _log.setLevel(logging.DEBUG)
    elif args.vb > 0:
        _log.setLevel(logging.INFO)

    # Read and parse configuration file
    _log.debug(f"reading settings file '{args.config}'")
    try:
        conf_file = open(args.config, "r")
        settings = Settings(conf_file)
    except IOError as err:
        _log.fatal(f"Abort: cannot open settings file '{args.config}': {err}")
        return 1
    except Settings.ConfigError as err:
        _log.fatal(f"Abort: cannot read settings from '{args.config}': {err}")
        return 1

    # Process Jupyter Notebooks
    if (args.mode == MODE_REL) or (not args.no_notebooks):
        if args.mode == MODE_REL and args.no_notebooks:
            notify(
                "Warning: 'release' mode does not honor the --no-notebooks flag. "
                "Building notebooks anyways, sorry.",
                1,
            )
        try:
            notify("Convert Jupyter notebooks")
            nbb = NotebookBuilder(settings)
            # set options based on mode
            clean, rebuild, test_mode, cont = None, None, None, None
            if args.mode == MODE_REL:
                clean, rebuild, test_mode, cont = True, True, False, False
            elif args.mode == MODE_DEV:
                clean, rebuild, test_mode, cont = False, False, False, True
            elif args.mode == MODE_TEST:
                clean, rebuild, test_mode, cont = False, True, True, True
            nbb.build(
                dict(
                    rebuild=rebuild,
                    test_mode=test_mode,
                    clean=clean,
                    continue_on_error=cont,
                )
            )
        except NotebookError as err:
            _log.fatal(f"Could not build notebooks: {err}")
            return -1
        nbb.report()

    # Build Sphinx documentation
    if args.mode != MODE_TEST:
        try:
            notify("Build documentation with Sphinx")
            spb = SphinxBuilder(settings)
            spb.build({})
        except SphinxError as err:
            _log.fatal(f"Could not build Sphinx docs: {err}")
            return -1

    return 0


if __name__ == "__main__":
    sys.exit(main())
