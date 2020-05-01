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
            "strip": False,
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
    TEST_SUFFIXES = ("_test", "_testing")  # test notebook suffixes
    REMOVE_CELL_TAG = "remove_cell"

    JUPYTER_NB_VERSION = 4  # for parsing

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
        self._cleaned = []  # remember generated "clean" entries

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
        self._results = self.Results()
        self._results.start()
        for item in nb_dirs:
            self._build_tree(item)
        self._results.stop()
        return self._results

    def remove_generated_files(self):
        if not self._cleaned:
            return 0
        n = 0
        for c in self._cleaned:
            _log.debug(f"remove: {c}")
            c.unlink()
            n += 1
        return n

    def report(self):
        """Print some messages to the user.
        """
        r = self._results  # alias
        total = r.n_fail + r.n_success
        notify(
            f"Processed {len(r.dirs_processed)} directories: "
            f"cached={len(r.cached)}, "
            f"converted={r.n_success}/{total}, "
            f"duration={r.duration:.3f}s",
            level=1,
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
                os.makedirs(self._imgdir, exist_ok=True)
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
                    self._results.cached.append(entry)
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
        # optionally just execute notebook and stop
        if self.s.get("test_mode"):
            self._parse_and_execute(entry)
            return
        # optionally strip special cells
        if self.s.get("strip") and self._has_tagged_cells(entry):
            _log.debug(f"notebook '{entry}' has test cell(s)")
            entries, tmpdir = self._strip_tagged_cells(entry)
        else:
            entries, tmpdir = [entry], None
        # main loop
        try:
            for e in entries:  # notebooks to export
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
            self._create_notebook_wrapper_page(entries[0].stem, entry.stem, outdir)
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
        """Strip specially tagged cells from a notebook.
        Copy notebook to a temporary location, and generate a stripped
        version there. No files are modified in the original directory.
        """
        tmpdir = Path(tempfile.mkdtemp())  # already checked this
        _log.debug(f"run notebook in temporary directory: {tmpdir}")
        # Copy notebook to temporary directory
        raw_entry = tmpdir / entry.name
        shutil.copy(entry, raw_entry)
        # Remove the special tags
        (body, resources) = NotebookExporter(
            config=self._nb_remove_config
        ).from_filename(str(raw_entry))
        # Determine outbook notebook name:
        # either strip test suffix, or append "clean" suffix
        cleaned_name = None
        for suffix in self.TEST_SUFFIXES:
            if entry.stem.endswith(suffix):
                cleaned_name = entry.stem[: entry.stem.rfind("_")]
                break
        if cleaned_name is None:
            cleaned_name = f"{entry.stem}{self.CLEAN_SUFFIX}"
        # Create the new notebook
        wrt = nbconvert.writers.FilesWriter()
        wrt.build_directory = str(entry.parent)
        wrt.write(body, resources, notebook_name=cleaned_name)
        _log.debug(
            f"stripped tags from '{raw_entry}' -> '{cleaned_name}' in "
            f"{wrt.build_directory}"
        )
        # Return both notebook names, and temporary directory (for cleanup)
        cleaned_entry = entry.parent / f"{cleaned_name}.ipynb"
        self._cleaned.append(cleaned_entry)
        return [cleaned_entry, raw_entry], tmpdir

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

    def _create_notebook_wrapper_page(
        self, nb_file: str, nb_base_file: str, output_dir: Path
    ):
        """Generate a Sphinx documentation page for the Module.
        """
        # interpret some characters in filename differently for title
        title = nb_base_file.replace("_", " ")
        title_under = "=" * len(title)
        # create document from template
        doc = self.s.get("Template").substitute(
            title=title, notebook_name=nb_file, title_underline=title_under
        )
        # write out the new doc
        doc_rst = output_dir / (nb_base_file + "_doc.rst")
        with doc_rst.open("w") as f:
            _log.info(f"generate Sphinx doc wrapper for {nb_file} => {doc_rst}")
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
        doc_dir = self.root_path() / self.s.get("paths.output")
        if not os.path.exists(html_dir):
            _log.warning(f"Target HTML directory {html_dir} does not exist: creating")
            html_dir.mkdir(parents=True)
        # copy images
        notify(
            f"Copying image files from '{self.s.get('paths.source')}' -> " f"'{doc_dir}'",
            1,
        )
        self._copy_aux(doc_dir, ["png", "jpg", "jpeg", "pdf"])
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
        # copy notebooks
        notify(f"Copying notebooks from '{self.s.get('paths.source')}' -> '{html_dir}'", 1)
        self._copy_aux(html_dir, ["ipynb"])

    def _copy_aux(self, dest, ext_list):
        """Copy auxiliary files in 'src' into a built directory.
        """
        root = self.root_path()
        for nbdir in self.s.get("notebook.directories"):
            source, output = nbdir["source"], nbdir["output"]
            src_dir = root / self.s.get("paths.source") / source
            files = []
            for ext in ext_list:
                files.extend(list(Path(src_dir).glob(f"**/*.{ext}")))
            for aux_file in files:
                if any((part.startswith(".") for part in aux_file.parts)):
                    continue
                copy_to = dest / output / aux_file.relative_to(src_dir)
                _log.info(f"copy: {aux_file} -> {copy_to}")
                try:
                    shutil.copy(str(aux_file), str(copy_to))
                except IOError as err:
                    _log.warning(f"Copy failed: {err}")

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

    CLEAN_SUFFIX = NotebookBuilder.CLEAN_SUFFIX

    def build(self, options):
        root = self.root_path()
        docs_path = root / self.s.get("paths.output")
        # clean in each path
        did_remove = False
        did_remove |= self._clean_html(docs_path / self.s.get("paths.html"))
        did_remove |= self._clean_docs(docs_path)
        did_remove |= self._clean_src(root / self.s.get("paths.source"))
        if not did_remove:
            notify("No files removed", 1)

    def _clean_html(self, html_path):
        if html_path.exists():
            notify(f"remove directory: {html_path}", 1)
            shutil.rmtree(html_path)
            return True
        return False

    def _clean_docs(self, docs_path):
        removed_any = False
        stop_list = set(["index.rst"])
        for notebook_dirs in self.s.get("notebook.directories"):
            docs_output = docs_path / notebook_dirs["output"]
            for file_type in "rst", "html", "png", "jpg":
                for f in docs_output.glob(f"**/*.{file_type}"):
                    if not f.name in stop_list:
                        notify(f"remove: {f}", 1)
                        f.unlink()
                        removed_any = True
        return removed_any

    def _clean_src(self, src_path):
        removed_any = False
        for notebook_dirs in self.s.get("notebook.directories"):
            src_input = src_path / notebook_dirs["source"]
            for f in src_input.glob(f"**/*{self.CLEAN_SUFFIX}.ipynb"):
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
    c = [Color.MAGENTA, Color.GREEN, Color.WHITE][min(level, 2)]
    if level == 0:
        print()
    print(f"{c}{message}{Color.RESET}")


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
    message = (
        "\n"
        "# tl;dr - To build the documentation, running Jupyter Notebooks\n"
        "# as needed, etc., use this command: ./build.py -drs\n"
        "\n"
        "The build.py command is used to create the documentation from\n"
        "the Jupyter Notebooks and hand-written '.rst' files in this\n"
        "repository. It is also used by tests, programmatically, to run\n"
        "the notebooks. Some sample command-lines, with comments as to \n"
        "what they will do, follow. The operation of this script is also\n"
        "controlled by the configuration file, 'build.yml' by default.\n"
        "Options specified on the command line will override options of\n"
        "the same name in the configuration file.\n"
        "\n"
        "\n"
        "# Read the default configuration file, but do nothing else\n"
        "./build.py\n"
        "\n"
        "# Build the Sphinx documentation, only\n"
        "./build.py --docs\n"
        "./build.py -d  # <-- short option\n"
        "\n"
        "# Run and convert Jupyter notebooks. Only those notebooks\n"
        "# that have not changed since the last time this was run will\n"
        "# be re-executed. Converted notebooks are stored in the 'docs'\n"
        "# directory, in locations configured in the 'build.yml'\n"
        "# configuration file.\n"
        "./build.py --notebooks\n"
        "./build.py -r  # <-- short option\n"
        "\n"
        "# Run and convert Jupyter notebooks, as in previous command,\n"
        "# but also strip any cells marked as 'test' cells, and run those \n"
        "# notebooks as well.\n"
        "./build.py -rs\n"
        "\n"
        "# Run and convert Jupyter notebooks, as in previous command,\n"
        "# but also force a rebuild of all notebooks whether or not they\n"
        "# were changed since the last time they were converted.\n"
        "./build.py -Rrs\n"
        "\n"
        "# Build documentation, run and convert notebooks, both testing and\n"
        "# test-cell-stripped versions.\n"
        "./build.py -drs\n"
        "\n"
        "# Remove all built documentation files\n"
        "./build.py --clean\n"
        "./build.py -c  # <-- short option\n"
        "\n"
        "# Run with <options> at different levels of verbosity\n"
        "./build.py <options>      # Show warning, error, fatal messages\n"
        "./build.py -v <options>   # Add informational (info) messages\n"
        "./build.py -vv <options>  # Add debug messages\n"
        "\n"
    )
    print(message)


def main():
    ap = argparse.ArgumentParser(
        description="Build documentation and/or Jupyter notebooks"
    )
    ap.add_argument(
        "--config",
        default="build.yml",
        help="Location of configuration file (default=./build.yml)",
    )
    ap.add_argument(
        "--usage", "-U", action="store_true", help="Print a more detailed usage message"
    )
    ap.add_argument("--clean", "-c", action="store_true", help="Clean built objects")
    ap.add_argument("--docs", "-d", action="store_true", help="Build documentation")
    ap.add_argument("--exit", "-x", action="store_true", help="Exit on first error")
    ap.add_argument(
        "--notebooks",
        "-r",
        action="store_true",
        help="Run/convert" " Jupyter notebooks",
    )
    ap.add_argument(
        "--rebuild",
        "-R",
        action="store_true",
        help="For Jupyter notebooks, re-run even if no change",
    )
    ap.add_argument(
        "--strip",
        "-s",
        action="store_true",
        help="For Jupyter notebooks, strip cells marked for tests",
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
    exit_on_error = args.exit
    run_notebooks = args.notebooks
    rebuild_notebooks = args.rebuild
    strip_notebooks = args.strip
    build_docs = args.docs
    clean_files = args.clean

    if clean_files:
        notify("Clean all built files")
        cleaner = Cleaner(settings)
        cleaner.build({})

    nbb = None
    if run_notebooks:
        notify("Convert Jupyter notebooks")
        nbb = NotebookBuilder(settings)
        try:
            nbb.build(
                {
                    "rebuild": rebuild_notebooks,
                    "continue_on_error": not exit_on_error,
                    "strip": strip_notebooks,
                    "test_mode": False,
                }
            )
        except NotebookError as err:
            _log.fatal(f"Could not build notebooks: {err}")
            return -1
        nbb.report()

    if build_docs:
        notify("Build documentation with Sphinx")
        spb = SphinxBuilder(settings)
        try:
            spb.build({})
        except SphinxError as err:
            _log.fatal(f"Could not build Sphinx docs: {err}")
            return -1

    # cleanup, if the notebooks were run
    if nbb:
        notify("Removing any generated notebook files")
        n = nbb.remove_generated_files()
        notify(f"Removed {n} files", 1)

    return 0


if __name__ == "__main__":
    sys.exit(main())
