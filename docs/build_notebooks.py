#!/usr/bin/env python
import argparse
import logging
import pathlib
import re
import shutil
import sys

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


IMAGE_SUFFIXES = ".jpg", ".jpeg", ".png", ".gif", ".svg", ".pdf"


def convert(
    srcdir: pathlib.Path = None,
    outdir: pathlib.Path = None,
    htmldir: pathlib.Path = None,
    test: bool = False,
    options: dict = None,
):
    """Convert notebooks under `srcdir`, placing output in `outdir`.

    Args:
        srcdir: Input directory
        outdir: Output directory
        htmldir: Where HTML files will end up (Sphinx build directory)
        test: Run in 'test' mode
        options: Options controlling the conversion:
            * "continue": if true, continue if there is an error executing the
              notebook; if false, raise NotebookError
            * "kernel": name of kernel for notebook execution, ie conda env name
            * "format": output format, either "html" or "rst"
            * "pat": if not None, a `re` expression that must match notebook
               filenames, or they will be skipped
    Raises:
        NotebookError: error executing or parsing the notebook
    Returns:
        None
    """
    ep_kw = {}
    if options["kernel"]:
        ep_kw["kernel_name"] = options["kernel"]
    ep = ExecutePreprocessor(timeout=600, **ep_kw)
    if options["format"] == "html":
        exp = HTMLExporter()
    elif options["format"] == "rst":
        exp = RSTExporter()
    else:
        raise ValueError(f"Invalid output format: {options['fmt']}")
    wrt = FilesWriter()

    test_errors = [] if test else None
    _convert(srcdir, outdir, htmldir, wrt, exp, ep, options, test_errors)

    return test_errors


def _convert(srcdir, outdir, htmldir, wrt, exp, ep, options, test_errors):
    """Recurse through directory, converting notebooks.
    """
    _log.debug(f"looking for notebooks in {srcdir}")

    if test_errors is None:  # conversion mode (not test-only)
        wrt.build_directory = str(outdir)

    for entry in srcdir.iterdir():
        filename = entry.parts[-1]
        if filename.startswith("."):
            _log.debug(f"skip: {entry}")
            continue  # e.g. .ipynb_checkpoints
        if entry.is_dir():
            new_outdir = outdir / entry.parts[-1]
            new_htmldir = htmldir / entry.parts[-1]
            _convert(entry, new_outdir, new_htmldir, wrt, exp, ep, options, test_errors)
        elif entry.suffix == ".ipynb":
            if options["pat"] and not options["pat"].search(filename):
                _log.debug(f"does not match {options['pat']}, skip {entry}")
                continue

            if test_errors is not None:  # testing mode
                _log.info(f"testing '{entry}'")
                try:
                    nb = nbformat.read(str(entry), as_version=4)
                    try:
                        ep.preprocess(nb, {"metadata": {"path": str(entry.parent)}})
                    except (CellExecutionError, NameError) as err:
                        test_errors.append(f"'{entry}' execution error: {err}")
                except (nbformat.reader.NotJSONError, AttributeError) as err:
                    test_errors.append(f"'{entry}' parse error: {err}")
                continue

            _log.info(f"converting '{entry}' to {options['format']}")

            # parse
            _log.debug(f"parsing '{entry}'")
            nb = None
            try:
                nb = nbformat.read(str(entry), as_version=4)
            except nbformat.reader.NotJSONError:
                if options["continue"]:
                    _log.error(f"parsing of '{entry}' failed: invalid JSON")
                else:
                    raise NotebookFormatError(f"'{entry}' is not JSON")
            except AttributeError as err:
                if options["continue"]:
                    _log.error(f"parsing of '{entry}' failed: format error, {err}")
                else:
                    raise NotebookFormatError(f"'{entry}' format error: {err}")

            if nb is None:
                _log.warning(f"skip conversion of '{entry}'")
            else:
                # execute
                _log.debug(f"executing '{entry}'")
                try:
                    ep.preprocess(nb, {"metadata": {"path": str(entry.parent)}})
                except (CellExecutionError, NameError) as err:
                    if options["continue"]:
                        _log.error(f"execution of '{entry}' failed: {err}")
                        _log.warning(f"generating partial output for '{entry}'")
                    else:
                        raise NotebookExecError(f"execution error for '{entry}': {err}")

                # export
                _log.debug(f"exporting '{entry}'")
                (body, resources) = exp.from_notebook_node(nb)
                wrt.write(body, resources, notebook_name=entry.stem)
        elif test_errors is None and entry.suffix in IMAGE_SUFFIXES:
            _log.debug(f"copying image '{entry}' to html output dir '{htmldir}'")
            shutil.copy(entry, htmldir)

    _log.debug(f"leaving {srcdir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source_dir")
    ap.add_argument("output_dir")
    ap.add_argument("--html-dir")
    ap.add_argument("--continue", dest="cont", action="store_true", default=False)
    ap.add_argument("--format", dest="fmt", choices=["html", "rst"], default="html")
    ap.add_argument("--kernel", dest="kernel", default="")
    ap.add_argument("--match", dest="match", default=None)
    ap.add_argument("--test", dest="test", action="store_true", default=False)
    ap.add_argument("-v", "--verbose", action="count", dest="vb", default=0)
    args = ap.parse_args()

    if args.vb > 1:
        _log.setLevel(logging.DEBUG)
    elif args.vb > 0:
        _log.setLevel(logging.INFO)

    srcdir = pathlib.Path(args.source_dir)
    if not srcdir.exists():
        _log.fatal(f"source directory does not exist: {srcdir}")
        return -1

    outdir = pathlib.Path(args.output_dir)
    if not outdir.exists():
        _log.warning(f"output directory does not exist: {outdir}")

    if args.html_dir is None:
        htdir = pathlib.Path("build") / outdir
    else:
        htdir = pathlib.Path(args.html_dir)
    if not htdir.exists():
        _log.warning(f"HTML (build) directory does not exist: {htdir}")

    options = {
        "continue": args.cont,
        "kernel": args.kernel,
        "format": args.fmt,
        "pat": re.compile(args.match) if args.match else None,
    }
    test_mode = args.test

    status_code, te = 0, None
    try:
        te = convert(
            srcdir=srcdir, outdir=outdir, htmldir=htdir, options=options, test=test_mode
        )
    except ValueError as err:
        _log.fatal(f"error converting notebooks: {err}")
        status_code = 1
    except NotebookError as err:
        _log.fatal(f"error running notebook: {err}")
        status_code = 2

    if test_mode:
        if te:
            errlist = "\n".join([f"[ERROR {i + 1}] {s}" for i, s in enumerate(te)])
            _log.error(f"{len(te):d} test error(s):\n{errlist}")
            status_code = len(te)
        else:
            _log.info(f"all tests passed")
            status_code = 0

    return status_code


if __name__ == "__main__":
    sys.exit(main())
