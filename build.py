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

Some example command-lines::

    # build the docs, only re-running notebooks whose converted output is "stale",
    # i.e. the notebook changed more recently
    ./build.py

    # just re-run the Sphinx documentation generator
    ./build.py -N

    # rebuild the notebooks, though this will also convert them to HTML
    # output. Do not run Sphinx. Rebuild *all* the notebooks, whether they
    # are "stale" or not.
    ./build.py -SR

    # run with a custom configuration file
    ./build.py -c my-config.yaml

    # don't really do anything except read/parse configuration file; so this
    # can serve as a test of that file's syntax
    ./build.py -c my-config.yaml -NS

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
import time
from typing import TextIO, List
import yaml

#
import nbformat
from nbconvert.exporters import HTMLExporter, RSTExporter
from nbconvert.writers import FilesWriter
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError

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
            "error-file": "notebook-errors.txt",
            "continue": True,
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
                 you have to explicitly pass default=None (or a value of your
                 choosing) to get the default value to be used.

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
                raise KeyError(f"unknown value '{val}' in section '{section}'")
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
                f"unknown section '{section}' not in: " f"{', '.join(self.DEFAULTS)}"
            )
        return section, val

    def _subst_paths(self, value):
        if hasattr(value, "format"):
            return value.format(**self.d["paths"])
        return value


class Builder(ABC):
    """Abstract base class for notebook and sphinx builders.
    """

    def __init__(self, settings):
        self.s = settings

    @abstractmethod
    def build(self, options):
        pass

    def _merge_options(self, options):
        for key, value in options.items():
            self.s.set(key, value)


class NotebookBuilder(Builder):
    """Run Jupyter notebooks and render them for viewing.
    """

    # Map format to file extension
    FORMATS = {"html": ".html", "rst": ".rst"}
    HTML_IMAGE_DIR = "_images"

    @dataclass
    class Results:
        """Stores results from build().
        """
        _t = None
        failed: List[str]
        dirs_processed: List[str]
        total: int = 0
        duration: float = -1.0

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

    def build(self, options):
        self.s.set_default_section("notebook")
        self._merge_options(options)
        self._open_error_file()
        self._ep = self._create_preprocessor()
        self._imgdir = (
            Path(self.s.get("paths.output"))
            / self.s.get("paths.html")
            / self.HTML_IMAGE_DIR
        )
        if _log.isEnabledFor(logging.DEBUG):
            _log.debug("settings for build: {self.s}")
        self._read_template()
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

        notify(
            f"{r.total} notebooks in {len(r.dirs_processed)} "
            f"directories in {r.duration:.3f} seconds"
        )
        notify("Notebook build summary", level=1)
        if r.total == 0 and len(r.failed) == 0:
            notify(f"No notebooks converted", level=2)
        elif len(r.failed) == 0:
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
        nb_template_path = Path(self.s.get("template"))
        try:
            with nb_template_path.open("r") as f:
                nb_template = Template(f.read())
        except IOError as err:
            raise NotebookError(f"cannot read template file {nb_template_path}: {err}")
        self.s.set("Template", nb_template)

    def _build_tree(self, info: dict):
        """"""
        ex_msg = "" if self.s.get("execute", True) else " (without executing)"
        try:
            source, output = info["source"], info["output"]
            match = info.get("match", None)
        except KeyError:
            raise NotebookError(
                f"notebook directory requires values for 'source', "
                f"'output', but got: {info}"
            )
        sroot, oroot = (
            Path(self.s.get("paths.source")),
            Path(self.s.get("paths.output")),
        )
        srcdir, outdir = sroot / source, oroot / output
        notify(f"Converting notebooks{ex_msg} in '{srcdir}'", level=1)
        self._match_expr = re.compile(match) if match else None
        # build, starting at this directory
        n = self._build_subtree(srcdir, outdir, depth=1)
        self._results.total += n
        self._results.dirs_processed.append(srcdir)

    def _build_subtree(self, srcdir: Path, outdir: Path, depth: int) -> int:
        """Build all notebooks in a given directory, exporting into the given output
        format.
        """
        _log.debug(f"build.begin subtree={srcdir}")
        total = 0
        for entry in srcdir.iterdir():
            filename = entry.parts[-1]
            if filename.startswith(".") or filename.startswith("__"):
                _log.debug(f"skip special file '{entry}'")
                continue  # e.g. .ipynb_checkpoints
            if entry.is_dir():
                # build sub-directory (filename is really directory name)
                total += self._build_subtree(entry, outdir / filename, depth + 1)
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
                    _log.info(f"converting notebook '{entry}'")
                    continue_on_err = self.s.get("continue", None)
                    try:
                        self._convert(entry, outdir, depth)
                        total += 1
                    except NotebookExecError as err:
                        self._notebook_error(entry, err)
                        if continue_on_err:
                            _log.warning(f"generating partial output for '{entry}'")
                            continue
                        raise
                    except NotebookError as err:
                        if continue_on_err:
                            _log.warning(f"failed to convert {entry}: {err}")
                            continue
                        raise  # abort all processing
        _log.debug(f"build.end subtree={srcdir} n={total}")
        return total

    def _notebook_error(self, entry, error):
        filename = str(entry)
        self._nb_error_file.write(f"\n====> File: {filename}\n")
        self._nb_error_file.write(str(error))
        self._nb_error_file.flush()  # in case someone is tailing the file
        self._results.failed.append(filename)

    def _convert(self, entry: Path, outdir: Path, depth: int):
        """Convert a notebook.

        Parameters
        ----------
        entry: notebook to convert
        outdir: output directory for .html and .rst files
        depth: depth below root, for fixing image paths
        """
        # parse
        _log.debug(f"parsing '{entry}'")
        try:
            nb = nbformat.read(str(entry), as_version=4)
        except nbformat.reader.NotJSONError:
            raise NotebookFormatError(f"'{entry}' is not JSON")
        except AttributeError:
            raise NotebookFormatError(f"'{entry}' has invalid format")
        # execute
        if self.s.get("execute"):
            _log.debug(f"executing '{entry}'")
            try:
                self._ep.preprocess(nb, {"metadata": {"path": str(entry.parent)}})
            except (CellExecutionError, NameError) as err:
                raise NotebookExecError(f"execution error for '{entry}': {err}")
        # export
        if self.s.get("test_mode"):
            return  # skip export
        _log.debug(f"exporting '{entry}'")
        wrt = FilesWriter()
        for (exp, postprocess_func, pp_args) in (
            (RSTExporter(), self._postprocess_rst, ()),
            (HTMLExporter(), self._postprocess_html, (depth,)),
        ):
            (body, resources) = exp.from_notebook_node(nb)
            body = postprocess_func(body, *pp_args)
            wrt.build_directory = str(outdir)  # write converted NB to output dir
            wrt.write(body, resources, notebook_name=entry.stem)
        # create a 'wrapper' page
        self._create_notebook_wrapper_page(entry.stem, outdir)

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
        args = self.s.get("args").split()
        html_dir = Path(self.s.get("paths.output")) / self.s.get("paths.html")
        if not os.path.exists(html_dir):
            _log.warning("Target HTML directory {html_dir} does not exist: creating")
            html_dir.mkdir(parents=True)
        errfile = self.s.get("error_file")
        cmdargs = ["sphinx-build", "-a", "-w", errfile] + args
        cmdline = " ".join(cmdargs)
        _log.info(f"Running Sphinx command: {cmdline}")
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-N", "--no-notebooks", action="store_true")
    ap.add_argument("-E", "--no-exec", action="store_true")
    ap.add_argument("-S", "--no-sphinx", action="store_true")
    ap.add_argument("-c", "--config", default="build.yml")
    ap.add_argument("-k", "--kernel", default=None, dest="kernel")
    ap.add_argument(
        "-R",
        "--rebuild-all",
        default=False,
        action="store_true",
        help="Rebuild all notebooks",
    )
    ap.add_argument("-v", "--verbose", action="count", dest="vb", default=0)
    args = ap.parse_args()

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
    if not args.no_notebooks:
        try:
            notify("Convert Jupyter notebooks")
            nbb = NotebookBuilder(settings)
            nbb.build(
                dict(
                    kernel=args.kernel,
                    execute=not args.no_exec,
                    rebuild=args.rebuild_all,
                )
            )
        except NotebookError as err:
            _log.fatal(f"Could not build notebooks: {err}")
            return -1
        nbb.report()

    # Build Sphinx documentation
    if not args.no_sphinx:
        try:
            notify("Sphinx documentation")
            spb = SphinxBuilder(settings)
            spb.build({})
        except SphinxError as err:
            _log.fatal(f"Could not build Sphinx docs: {err}")
            return -1

    return 0


if __name__ == "__main__":
    sys.exit(main())
