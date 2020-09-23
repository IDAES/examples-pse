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
# stdlib
from abc import ABC, abstractmethod
import argparse
from collections import namedtuple
from datetime import datetime
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
from typing import List, TextIO, Tuple, Optional
import yaml

# third-party
import nbconvert
from nbconvert.exporters import HTMLExporter, RSTExporter, NotebookExporter
from nbconvert.writers import FilesWriter
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
import nbformat
from traitlets.config import Config

# from build.py dir
_script_dir = os.path.abspath(os.path.dirname(__file__))
if _script_dir not in sys.path:
    print("Add script's directory to sys.path")
    sys.path.insert(0, _script_dir)
from build_util import bossy

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

class NotebookPreviouslyFailedError(NotebookError):
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
            "num_workers": 1,
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
        self._path = Path(f.name)

    @property
    def directory(self):
        return self._path.parent.absolute()

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

    def __init__(self, settings):
        self.s = settings
        self.root_path = settings.directory

    @abstractmethod
    def build(self, options):
        pass

    def _merge_options(self, options):
        for key, value in options.items():
            self.s.set(key, value)


Job = namedtuple("Job", ["nb", "tmpdir", "outdir", "depth"])


class NotebookBuilder(Builder):
    """Run Jupyter notebooks and render them for viewing.
    """

    TEST_SUFFIXES = ("_test", "_testing")  # test notebook suffixes
    HTML_IMAGE_DIR = "_images"

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
        self._test_mode, self._num_workers = None, 1
        self._match_expr = None
        # Lists of entries (or just one, for outdir) for each subdirectory
        self.notebooks_to_convert, self.data_files, self.outdir, self.depth = (
            {},
            {},
            {},
            {},
        )

    def build(self, options):
        self.s.set_default_section("notebook")
        self._num_workers = self.s.get("num_workers", default=4)
        self._merge_options(options)
        self._test_mode = self.s.get("test_mode")
        self._open_error_file()
        self._ep = self._create_preprocessor()
        self._imgdir = (
            self.root_path
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
            self.discover_tree(item)
        self.convert_discovered_notebooks()
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
        nb_template_path = self.root_path / self.s.get("template")
        try:
            with nb_template_path.open("r") as f:
                nb_template = Template(f.read())
        except IOError as err:
            raise NotebookError(f"cannot read template file {nb_template_path}: {err}")
        self.s.set("Template", nb_template)

    def discover_tree(self, info: dict):
        """Discover all the notebooks, recursively, in directories below `info['source']`
        and convert and build their output into `info['output']`.
        """
        try:
            source, output = info["source"], info["output"]
            match = info.get("match", None)
        except KeyError:
            raise NotebookError(
                f"notebook directory requires values for 'source', "
                f"'output', but got: {info}"
            )
        sroot, oroot = (
            self.root_path / self.s.get("paths.source"),
            self.root_path / self.s.get("paths.output"),
        )
        srcdir, outdir = sroot / source, oroot / output
        notify(f"Looking for notebooks in '{srcdir}'", level=1)
        self._match_expr = re.compile(match) if match else None
        # build, starting at this directory
        self.discover_subtree(srcdir, outdir, depth=1)
        self._results.dirs_processed.append(srcdir)

    def discover_subtree(self, srcdir: Path, outdir: Path, depth: int):
        """Discover all notebooks in a given directory.
        """
        _log.debug(f"Discover.begin subtree={srcdir}")

        # Iterate through directory and get list of notebooks to convert (and data files)
        notebooks_to_convert, data_files = [], []
        for entry in srcdir.iterdir():
            filename = entry.parts[-1]
            if filename.startswith(".") or filename.startswith("__"):
                _log.debug(f"skip special file '{entry}'")
                continue  # e.g. .ipynb_checkpoints
            if entry.is_dir():
                # build sub-directory (filename is really directory name)
                self.discover_subtree(entry, outdir / filename, depth + 1)
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
                if (
                    entry.suffix in DATA_SUFFIXES
                    or entry.suffix in CODE_SUFFIXES
                    or entry.suffix in IMAGE_SUFFIXES
                ):
                    data_files.append(entry)
        self.notebooks_to_convert[srcdir] = notebooks_to_convert
        self.data_files[srcdir] = data_files
        self.depth[srcdir], self.outdir[srcdir] = depth, outdir

    def convert_discovered_notebooks(self):
        # create temporary directories
        temporary_dirs = {
            srcdir: Path(tempfile.mkdtemp()) for srcdir in self.notebooks_to_convert
        }
        # copy datafiles (into temporary directories)
        for srcdir, dfs in self.data_files.items():
            tmpdir = temporary_dirs[srcdir]
            _log.debug(f"copy {len(dfs)} data file(s) into temp dir '{tmpdir}'")
            for fp in dfs:
                _log.debug(f"copy data file: {fp} -> {tmpdir}")
                shutil.copy(fp, tmpdir)
        # create list of jobs
        jobs = []
        for srcdir, nb_list in self.notebooks_to_convert.items():
            # template job
            tj = Job(
                None, temporary_dirs[srcdir], self.outdir[srcdir], self.depth[srcdir]
            )
            # per-notebook real jobs
            for nb in nb_list:
                jobs.append(Job(nb, tj.tmpdir, tj.outdir, tj.depth))
        # process list of jobs, in parallel
        _log.info(f"Process {len(jobs)} notebooks")
        num_workers = min(self._num_workers, len(jobs))
        # Run conversion in parallel
        _log.info(f"Convert notebooks with {num_workers} worker(s)")
        worker = ParallelNotebookWorker(
            processor=self._ep,
            wrapper_template=self.s.get("Template"),
            remove_config=self._nb_remove_config,
            test_mode=self._test_mode,
        )
        b = bossy.Bossy(
            jobs,
            num_workers=num_workers,
            worker_function=worker.convert,
            output_log=_log,
        )
        results = b.run()
        # Report results, which returns summary
        successes, failures = self._report_results(results)
        # Clean up any temporary directories
        for tmpdir in temporary_dirs.values():
            _log.debug(f"remove temporary directory at '{tmpdir.name}'")
            try:
                shutil.rmtree(str(tmpdir))
            except Exception as err:
                _log.error(f"could not remove temporary directory '{tmpdir}': {err}")
        # Record summary of success/fail and return
        _log.debug(f"Convert.end {successes}/{successes + failures}")
        self._results.n_fail += failures
        self._results.n_success += successes

    def _report_results(
        self, result_list: List["ParallelNotebookWorker.ConversionResult"]
    ) -> Tuple[int, int]:
        s, f = 0, 0
        for worker_id, result in result_list:
            if result.ok:
                s += 1
            else:
                f += 1
                filename = str(result.entry)
                self._write_notebook_error(self._nb_error_file, filename, result.why)
                self._results.failed.append(filename)
        return s, f

    @staticmethod
    def _write_notebook_error(error_file, nb_filename, error):
        error_file.write(f"\n====> File: {nb_filename}\n")
        error_file.write(str(error))
        error_file.flush()  # in case someone is tailing the file


class ParallelNotebookWorker:
    """For converting notebooks.

    State is avoided where possible to ensure that the ForkingPickler succeeds.

     Main method is `convert`.
     """

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

    # Map format to file extension
    FORMATS = {"html": ".html", "rst": ".rst"}

    # Mapping from {meaning: tag_name}
    CELL_TAGS = {
        "remove": "remove_cell",
        "exercise": "exercise",
        "testing": "testing",
        "solution": "solution",
    }
    JUPYTER_NB_VERSION = 4  # for parsing

    class ConversionResult:
        def __init__(self, id_, ok, converted, why, entry):
            self.id_, self.ok, self.converted, self.why, self.entry = (
                id_,
                ok,
                converted,
                why,
                entry,
            )

        def __str__(self):
            success = "SUCCESS" if self.ok else "FAILURE"
            cached = " (cached)" if not self.converted else ""
            common_prefix = f"{success} for '{self.entry}'{cached}"
            if not self.ok:
                summary = f"{common_prefix}: {self.why}"
            else:
                summary = common_prefix
            return summary

    def __init__(
        self,
        processor: ExecutePreprocessor = None,
        wrapper_template: Optional[Template] = None,
        remove_config: Optional[Config] = None,
        test_mode: bool = False,
    ):
        self.processor = processor
        self.template, self.rm_config = (
            wrapper_template,
            remove_config,
        )
        self.test_mode = test_mode
        self.log_q, self.id_ = None, 0

    # Logging utility functions

    def log(self, level, msg):
        self.log_q.put((level, f"[Worker {self.id_}] {msg}"))

    def log_error(self, msg):
        return self.log(logging.ERROR, msg)

    def log_warning(self, msg):
        return self.log(logging.WARNING, msg)

    def log_info(self, msg):
        return self.log(logging.INFO, msg)

    def log_debug(self, msg):
        return self.log(logging.DEBUG, msg)

    # Main function

    def convert(self, id_, job, log_q) -> ConversionResult:
        """Parallel 'worker' to convert a single notebook.
        """
        self.log_q, self.id_ = log_q, id_

        self.log_info(f"Convert notebook name={job.nb}: begin")

        ok, why = True, ""
        if not self.test_mode:
            if not job.outdir.exists():
                job.outdir.mkdir(parents=True)
        # build, if the output file is missing/stale
        verb = "Running" if self.test_mode else "Converting"
        self.log_info(f"{verb}: {job.nb.name}")
        # continue_on_err = self.s.get("continue_on_error", None)
        converted = False
        try:
            converted = self._convert(job)
        except NotebookPreviouslyFailedError as err:
            ok, why = False, f"Previously failed {err}"
        except NotebookExecError as err:
            ok, why = False, err
            self._write_failed_marker(job)
            self.log_error(
                f"Execution failed: generating partial output for '{job.nb}'"
            )
        except NotebookError as err:
            ok, why = False, f"NotebookError: {err}"
            self.log_error(f"Failed to convert {job.nb}: {err}")
        except Exception as err:
            ok, why = False, f"Unknown error: {err}"
            self.log_error(f"Failed due to unknown error: {err}")

        if converted:
            # remove failed marker, if there was one from a previous execution
            failed_marker = self._get_failed_marker(job)
            if failed_marker:
                self.log_info(
                    f"Remove stale marker of failed execution: {failed_marker}"
                )
                failed_marker.unlink()

        self.log_info(f"Convert notebook name={job.nb}: end, ok={ok}")

        return self.ConversionResult(id_, ok, converted, why, job.nb)

    def _convert(self, job) -> bool:
        """Convert a notebook.

        Returns:
            True if conversion was performed, False if no conversion was needed
        """
        info, dbg = logging.INFO, logging.DEBUG  # aliases
        # strip special cells.
        if self._has_tagged_cells(job.nb, set(self.CELL_TAGS.values())):
            self.log_debug(f"notebook '{job.nb.name}' has test cell(s)")
            entry = self._strip_tagged_cells(job, ("remove", "exercise"), "testing")
            self.log_info(f"Stripped tags from: {job.nb.name}")
        else:
            # copy to temporary directory just to protect from output cruft
            entry = job.tmpdir / job.nb.name
            shutil.copy(job.nb, entry)

        # Convert all tag-stripped versions of the notebook.
        # Stop if failure marker is newer than source file.
        failed_time = self._previously_failed(job)
        if failed_time is not None:
            self.log_info(
                f"Skip notebook conversion, failure marker is newer, for: {entry.name}"
            )
            failed_datetime = datetime.fromtimestamp(failed_time)
            raise NotebookPreviouslyFailedError(f"at {failed_datetime}")
        # Stop if converted result is newer than source file.
        if self._previously_converted(job, entry):
            self.log_info(
                f"Skip notebook conversion, output is newer, for: {entry.name}"
            )
            return False
        self.log_info(f"Running notebook: {entry.name}")
        try:
            nb = self._parse_and_execute(entry)
        except NotebookExecError as err:
            self.log_error(f"Notebook execution failed: {err}")
            raise
        if self.test_mode:  # don't do export in test mode
            return True

        self.log_info(f"Exporting notebook '{entry.name}' to directory {job.outdir}")
        wrt = FilesWriter()
        # export each notebook into multiple target formats
        created_wrapper = False
        for (exp, post_process_func, pp_args) in (
            (RSTExporter(), self._postprocess_rst, ()),
            (HTMLExporter(), self._postprocess_html, (job.depth,)),
        ):
            self.log_debug(f"export '{job.nb}' with {exp} to notebook '{entry}'")
            (body, resources) = exp.from_notebook_node(nb)
            body = post_process_func(body, *pp_args)
            wrt.build_directory = str(job.outdir)
            wrt.write(body, resources, notebook_name=entry.stem)
            # create a 'wrapper' page
            if not created_wrapper:
                self.log_debug(
                    f"create wrapper page for '{entry.name}' in '{job.outdir}'"
                )
                self._create_notebook_wrapper_page(job, entry.stem)
                created_wrapper = True
            # move notebooks into docs directory
            self.log_debug(f"move notebook '{entry} to output directory: {job.outdir}")
            shutil.copy(entry, job.outdir / entry.name)
        return True

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
                        self.log_debug(f"Found tag '{tag}' in cell {i}")
                        return True  # can stop now, one is enough
        # no tagged cells
        return False

    def _previously_failed(self, job: Job) -> Optional[float]:
        orig = job.nb

        source_time = orig.stat().st_mtime

        # First see if there is a 'failed' marker. If so,
        # compare timestamps: if marker is newer than file, it's converted
        failed_file = job.outdir / (orig.stem + ".failed")
        if failed_file.exists():
            failed_time = failed_file.stat().st_mtime
            ac = failed_time > source_time
            if ac:
                self.log_info(
                    f"Notebook '{orig.stem}.ipynb' unchanged since previous failed conversion"
                )
            return failed_time
        return None

    def _previously_converted(self, job: Job, dest: Path) -> bool:
        """Check if any of the output files are either missing or older than
        the input file ('entry').

        Returns:
            True if the output is newer than the input, otherwise False.
        """
        orig = job.nb
        source_time = orig.stat().st_mtime

        # First see if there is a 'failed' marker. If so,
        # compare timestamps: if marker is newer than file, it's converted
        failed_file = job.outdir / (orig.stem + ".failed")
        if failed_file.exists():
            failed_time = failed_file.stat().st_ctime
            ac = failed_time > source_time
            if ac:
                self.log_info(
                    f"Notebook '{orig.stem}.ipynb' unchanged since previous failed conversion"
                )
            return ac

        # Otherwise look at all the output files and see if any one of them is
        # older than the source file (in which case it's NOT converted)
        for fmt, ext in self.FORMATS.items():
            output_file = job.outdir / f"{dest.stem}{ext}"
            self.log_debug(f"checking if cached: {output_file} src={orig}")
            if not output_file.exists():
                return False
            if source_time >= output_file.stat().st_mtime:
                return False

        return True

    def _write_failed_marker(self, job: Job):
        """Put a marker into the output directory for the failed notebook, so we
        can tell whether we need to bother trying to re-run it later.
        """
        if not job.outdir.exists():
            try:
                job.outdir.mkdir(parents=True)
            except Exception as err:
                self.log_error(
                    f"Could not write failed marker '{marker}' for entry={job.nb} "
                    f"outdir={job.outdir}: {err}"
                )
                return  # oh, well
        marker = job.outdir / (job.nb.stem + ".failed")
        self.log_debug(
            f"write failed marker '{marker}' for entry={job.nb} outdir={job.outdir}"
        )
        marker.open("w").write(
            "This file is a marker for avoiding re-runs of failed notebooks that haven't changed"
        )

    def _get_failed_marker(self, job: Job) -> Optional[Path]:
        marker = job.outdir / (job.nb.stem + ".failed")
        if marker.exists():
            self.log_debug(f"Found 'failed' marker: {marker}")
            return marker
        self.log_debug(
            f"No 'failed' marker for notebook '{job.nb}' in directory '{job.outdir}'"
        )
        return None

    def _strip_tagged_cells(self, job: Job, tags, remove_name: str):
        """Strip specially tagged cells from a notebook.

        Copy notebook to a temporary location, and generate a stripped
        version there. No files are modified in the original directory.

        Args:
            job: Notebook and parameters
            tags: List of tags (strings) to strip
            remove_name: Remove this component from the name

        Returns:
            stripped-entry, original-entry - both in the temporary directory
        """
        self.log_debug(f"run notebook in temporary directory: {job.tmpdir}")
        # Copy notebook to temporary directory
        tmp_nb = job.tmpdir / job.nb.name
        shutil.copy(job.nb, tmp_nb)
        # Remove the given tags
        # Configure tag removal
        tag_names = [self.CELL_TAGS[t] for t in tags]
        self.rm_config.TagRemovePreprocessor.remove_cell_tags = tag_names
        self.log_debug(
            f"removing tag(s) <{', '.join(tag_names)}'> from notebook: {job.nb.name}"
        )
        (body, resources) = NotebookExporter(config=self.rm_config).from_filename(
            str(tmp_nb)
        )
        # Determine output notebook name:
        # remove suffixes, either "Capitalized" or "lowercase" version,
        # with underscores before, after, or on both sides.
        suffixes = []
        for name in (remove_name.capitalize(), remove_name.lower()):
            suffixes.extend(["_" + name, name + "_", "_" + name + "_"])
        suffixes_disjunction = "|".join(suffixes)
        nb_name = re.sub(f"({suffixes_disjunction})", "", job.nb.stem)
        # Create the new notebook
        wrt = nbconvert.writers.FilesWriter()
        wrt.build_directory = str(job.tmpdir)
        self.log_debug(f"writing stripped notebook: {nb_name}")
        wrt.write(body, resources, notebook_name=nb_name)
        # Return both notebook names, and temporary directory (for cleanup)
        stripped_entry = job.tmpdir / f"{nb_name}.ipynb"
        return stripped_entry

    def _parse_and_execute(self, entry):
        # parse
        self.log_debug(f"parsing '{entry}'")
        try:
            nb = nbformat.read(str(entry), as_version=self.JUPYTER_NB_VERSION)
        except nbformat.reader.NotJSONError:
            raise NotebookFormatError(f"'{entry}' is not JSON")
        except AttributeError:
            raise NotebookFormatError(f"'{entry}' has invalid format")

        # execute
        self.log_debug(f"executing '{entry}'")
        try:
            self.processor.preprocess(nb, {"metadata": {"path": str(entry.parent)}})
        except (CellExecutionError, NameError) as err:
            raise NotebookExecError(f"execution error for '{entry}': {err}")
        return nb

    def _create_notebook_wrapper_page(self, job: Job, nb_file: str):
        """Generate a Sphinx documentation page for the Module.
        """
        # interpret some characters in filename differently for title
        title = nb_file.replace("_", " ").title()
        title_under = "=" * len(title)
        # create document from template
        doc = self.template.substitute(
            title=title, notebook_name=nb_file, title_underline=title_under
        )
        # write out the new doc
        doc_rst = job.outdir / (nb_file + "_doc.rst")
        with doc_rst.open("w") as f:
            self.log_info(f"generate Sphinx doc wrapper for {nb_file} => {doc_rst}")
            f.write(doc)

    def _postprocess_rst(self, body):
        return self._replace_image_refs(body)

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
        prefix = prefix / NotebookBuilder.HTML_IMAGE_DIR
        # split up content around <img> tags
        splits = re.split(r'<img src="(.*\.png"[^>]*)>', body)
        # replace grouping in odd-numbered splits with modified <img>
        for i in range(1, len(splits), 2):
            # orig = splits[i]
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
            self.log_warning("Could not insert stylesheet in HTML")
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
        args[-2] = str(self.root_path / args[-2])
        args[-1] = str(self.root_path / args[-1])
        html_dir = (
            self.root_path / self.s.get("paths.output") / self.s.get("paths.html")
        )
        doc_dir = self.root_path / self.s.get("paths.output")
        src_dir = self.root_path / self.s.get("paths.source")

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
                if not nb_dest.parent.exists():
                    nb_dest.parent.mkdir(parents=True)
                notify(f"Copy notebook {nb_path.name} -> {nb_dest.parent}", 2)
                shutil.copy(nb_path, nb_dest)
            _log.info(f"find supporting files in path: {nb_output_dir}")
            for ext in IMAGE_SUFFIXES:
                pattern = f"**/*{ext}"
                for nb_path in Path(nb_source_dir).glob(pattern):
                    nb_dest = html_dir / nb_path.relative_to(src_dir)
                    notify(
                        f"Copy supporting file {nb_path.name} -> {nb_dest.parent}", 3
                    )
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
        root = self.root_path
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
    ap.add_argument(
        "--remove", "-r", action="store_true", help="Remove generated files"
    )
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
    ap.add_argument(
        "--workers",
        "-w",
        dest="np",
        default=None,
        type=int,
        help="Number of parallel processes to run. Overrides `notebook.num_workers` in settings",
    )
    args = ap.parse_args()

    # Check for usage message flag
    if args.usage:
        print_usage()
        return 0

    # Check for confusing option combinations
    if args.convert and args.test_mode:
        ap.error("-t/--test conflicts with notebook conversion -c/--convert; pick one")
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

    # Update settings from command-line arguments
    if args.np:
        settings.set("notebook.num_workers", args.np)

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

    if run_notebooks or test_mode:
        verb = "Run" if test_mode else "Convert"
        notify(f"{verb} Jupyter notebooks")
        nbb = NotebookBuilder(settings)
        try:
            nbb.build({"test_mode": test_mode})
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
