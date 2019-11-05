#!/usr/bin/env python
import argparse
import logging
import pathlib
import sys
#
import nbformat
from nbconvert.exporters import HTMLExporter
from nbconvert.writers import FilesWriter
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError

_log = logging.getLogger("build_notebook_html")
_hnd = logging.StreamHandler()
_hnd.setFormatter(logging.Formatter("%(asctime)s %(levelname)s - %(message)s"))
_log.addHandler(_hnd)


class NotebookError(Exception):
    pass


def notebooks_to_html(srcdir=None, outdir=None, options=None):
    _log.debug(f"looking for notebooks in {srcdir}")

    hexp = HTMLExporter()
    wrt = FilesWriter()
    wrt.build_directory = str(outdir)
    ep_kw = {}
    if options["kernel"]:
        ep_kw["kernel_name"] = options["kernel"]
    ep = ExecutePreprocessor(timeout=600, **ep_kw)

    for entry in srcdir.iterdir():
        if entry.parts[-1].startswith("."):
            _log.debug(f"skip: {entry}")
            continue  # e.g. .ipynb_checkpoints
        if entry.is_dir():
            new_outdir = outdir / entry.parts[-1]
            notebooks_to_html(srcdir=entry, outdir=new_outdir, options=options)
        elif entry.suffix == ".ipynb":
            html_target = str(outdir / entry.stem) + ".html"
            _log.info(f"converting '{entry}' to '{html_target}'")
            nb = nbformat.read(str(entry), as_version=4)
            try:
                ep.preprocess(nb, {'metadata': {'path': str(entry.parent)}})
            except (CellExecutionError, NameError) as err:
                if options["continue"]:
                    _log.error(f"execution of '{entry}' failed: {err}")
                    _log.warning(f"continuing to generate partial HTML for '{entry}' "
                                 f"in '{html_target}'")
                else:
                    raise NotebookError(f"execution error for '{entry}': {err}")
            (body, resources) = hexp.from_notebook_node(nb)
            wrt.write(body, resources, notebook_name=entry.stem)

    _log.debug(f"leaving {srcdir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source_dir")
    ap.add_argument("output_dir")
    ap.add_argument("-v", "--verbose", action="count", dest="vb", default=0)
    ap.add_argument("--continue", dest="cont", action="store_true", default=False)
    ap.add_argument("--kernel", dest="kernel", default="")
    args = ap.parse_args()

    if args.vb > 1:
        _log.setLevel(logging.DEBUG)
    elif args.vb > 0:
        _log.setLevel(logging.INFO)
    srcdir = pathlib.Path(args.source_dir)
    if not srcdir.exists():
        _log.fatal(f"source directory does not exist: {srcdir}")
    outdir = pathlib.Path(args.output_dir)
    if not outdir.exists():
        _log.fatal(f"output directory does not exist: {outdir}")

    options = {"continue": args.cont, "kernel": args.kernel}

    status_code = 0
    try:
        notebooks_to_html(srcdir=srcdir, outdir=outdir, options=options)
    except ValueError as err:
        _log.fatal(f"error converting notebooks: {err}")
        status_code = 1
    except NotebookError as err:
        _log.fatal(f"error running notebook: {err}")
        status_code = 2
    return status_code


if __name__ == "__main__":
    sys.exit(main())
