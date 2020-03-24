#!/usr/bin/env python
"""
Build docs from Python.

This encapsulates generating Sphinx-ready versions of the Jupyter Notebooks and
calling 'sphinx-build' to build the HTML docs from source. No changes in the
Sphinx default build procedure should be required, which means this should work
on Windows, Mac OSX, and Linux equally well.

This program uses a configuration file, by default "build.yml" in the same
directory, to set some options and tell it which directories to process.
"""
import argparse
import logging
import os
import pathlib
import shutil
import subprocess
import re
from string import Template
import sys
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

g_doc_template = None
g_force_rebuild = False


def convert(
    srcdir: pathlib.Path = None,
    outdir: pathlib.Path = None,
    htmldir: pathlib.Path = None,
    options: dict = None,
):
    """Convert notebooks under `srcdir`, placing output in `outdir`.

    Args:
        srcdir: Input directory
        outdir: Output directory
        htmldir: Where HTML files will end up (Sphinx build directory)
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

    _convert(srcdir, outdir, htmldir, wrt, exp, ep, options)


def _convert(srcdir, outdir, htmldir, wrt, exp, ep, options):
    """Recurse through directory, converting notebooks.
    """
    _log.debug(f"looking for notebooks in {srcdir}")

    output_type = options["format"]  # rst or html
    for entry in srcdir.iterdir():
        filename = entry.parts[-1]
        if filename.startswith("."):
            _log.debug(f"skip: {entry}")
            continue  # e.g. .ipynb_checkpoints
        if entry.is_dir():
            new_outdir = outdir / entry.parts[-1]
            new_htmldir = htmldir / entry.parts[-1]
            _convert(entry, new_outdir, new_htmldir, wrt, exp, ep, options)
        elif entry.suffix == NOTEBOOK_SUFFIX:
            _log.debug(f"found notebook at: {entry}")
            nb = None
            if options["pat"] and not options["pat"].search(filename):
                _log.debug(f"does not match {options['pat']}, skip {entry}")
                continue
            # always guarantee target HTML directory exists
            if not os.path.exists(htmldir):
                _log.info(f"Creating output directory {htmldir} for notebook")
                os.makedirs(htmldir)
            # don't rebuild the notebook if a converted
            # one already exists in `outdir` that is as recent, UNLESS
            # the "g_force_rebuild" flag is turned on
            if not g_force_rebuild and converted_notebook_is_current(
                entry, outdir, output_type
            ):
                _log.info(f"Using previously-converted notebook in {outdir}")
                nb = True  # change from None so post-processing steps occur
            else:
                _log.info(f"converting '{entry}' to {options['format']}")

                # parse
                _log.debug(f"parsing '{entry}'")
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
                            raise NotebookExecError(
                                f"execution error for '{entry}': {err}"
                            )

                    # export
                    _log.debug(f"exporting '{entry}'")
                    (body, resources) = exp.from_notebook_node(nb)
                    if output_type == "rst":
                        body = replace_image_refs(body)
                    wrt.build_directory = str(
                        outdir
                    )  # write converted NB to output dir
                    wrt.write(body, resources, notebook_name=entry.stem)
            if nb is not None:
                # Do these regardless of whether we built the notebook or not.
                # This will restore files in case the HTML build directory was
                # removed/cleaned.
                copy_to_html_dir(entry, htmldir)
                # create Sphinx .rst doc page wrapping the notebook
                generate_doc_page(entry.stem, outdir)

        elif entry.suffix in IMAGE_SUFFIXES:
            copy_to_html_dir(entry, htmldir)

    _log.debug(f"leaving {srcdir}")


def converted_notebook_is_current(
    nb_file: pathlib.Path, dst_dir: pathlib.Path, otype: str
):
    output_nb_name = str(nb_file.stem) + "." + otype
    output_file = dst_dir / output_nb_name
    if not output_file.exists():
        return False
    nb_mod_ts = nb_file.stat().st_mtime
    output_mod_ts = output_file.stat().st_mtime
    return output_mod_ts >= nb_mod_ts


def generate_doc_page(nb_file: str, output_dir: pathlib.Path):
    """Generate a Sphinx documentation page for the Module.
    """
    # interpret some characters in filename differently for title
    title = nb_file.replace("__", ": ")
    title = title.replace("_", " ")
    # replacements in template
    doc_kw = {
        "title": title,
        "notebook_name": nb_file,
        "title_underline": "=" * len(title),
    }
    # fill in the template (use simple "string.substitute")
    doc = g_doc_template.substitute(**doc_kw)
    # write out the new doc
    doc_rst = output_dir / (nb_file + "_doc.rst")
    _log.debug(f"generate Sphinx doc wrapper for {nb_file} => {doc_rst}")
    with doc_rst.open("w") as f:
        _log.debug(f"Writing document to: {doc_rst}")
        f.write(doc)


def copy_to_html_dir(entry: str, htmldir: pathlib.Path):
    """Copy file in 'entry' over to HTML dir.

    Converted notebooks expect their image files in the same directory as the .html
    file, but Sphinx will copy them all into some _images/ sub-directory instead.
    """
    target_file = htmldir / entry
    if not target_file.exists():
        _log.debug(f"copying file '{entry}' to html output dir '{htmldir}'")
        shutil.copy(entry, htmldir)


def replace_image_refs(body):
    """Replace |image0| references with successive numbers, instead
        having them all be "image0". The order is "|image0|" then
        ".. |image0|".
    """
    chars = list(body)
    pos, imgn = 0, 0
    while True:
        pos = body.find("|image0|\n", pos)
        if pos == -1:
            break
        pos2 = body.find(".. |image0|", pos + 8)
        if pos2 == -1:
            raise ValueError(f"Couldn't find image matching ref at {pos}")
        if imgn > 0:
            # support up to 35 different images
            c = " 123456789abcdefghijklmnopqrstuvwxyz"[imgn]
            # replace '0' with another thing in ref
            chars[pos + 6] = c
            # replace '0' with same other thing in image
            chars[pos2 + 9] = c
        pos = pos2 + 10
        imgn += 1
    return "".join(chars)


def build_notebooks(config, **kwargs):
    """Build RST and HTML versions of Jupyter Notebooks, after
    executing them to ensure the output is populated.
    """
    global g_doc_template
    nb = config["notebook"]
    # pull options out of the config file
    kwargs.update({"continue": nb.get("continue", False)})
    source_base = nb.get("source_base", "src")
    output_base = nb.get("output_base", "docs")
    # read sphinx notebook page template into a global var
    if "template" not in nb:
        raise NotebookError("In configuration, `notebook.template` is missing")
    try:
        f = open(nb["template"], "r")
    except IOError as err:
        raise NotebookError(f"Cannot open notebook Sphinx template: {err}")
    g_doc_template = Template(f.read())
    # build all input/output directory pairs
    num_dirs = len(nb["directories"])
    for i, item in enumerate(nb["directories"]):
        srcdir = pathlib.Path(source_base) / item["source"]
        notify(f"Converting notebooks in '{srcdir}'", level=1)
        outdir = pathlib.Path(output_base) / item["output"]
        if "match" in item:
            kwargs["pat"] = re.compile(item["match"])
        else:
            kwargs["pat"] = None
        # build HTML and RST versions of the notebooks
        for ofmt in "html", "rst":
            htmldir = (
                pathlib.Path(output_base) / nb.get("html_dir", "build") / item["output"]
            )
            kwargs["format"] = ofmt
            try:
                convert(srcdir=srcdir, outdir=outdir, options=kwargs, htmldir=htmldir)
            except NotebookError as err:
                _log.error(
                    f"Failed converting notebooks in directory "
                    f"{i + 1} out of {num_dirs}. Abort."
                )
                raise
    _log.info(f"Converted {num_dirs} notebook directories")


def build_sphinx(config):
    """Run 'sphinx-build' command.
    """
    spx = config["sphinx"]
    args = spx["args"].split()
    html_dir = args[-1]
    if not os.path.exists(html_dir):
        _log.warning("Target HTML directory {html_dir} does not exist: creating")
        os.makedirs(html_dir)
    errfile = spx.get("error_file", "sphinx-errors.txt")
    cmdargs = ["sphinx-build", "-a", "-w", errfile] + args
    cmdline = " ".join(cmdargs)
    _log.info(f"Running Sphinx command: {cmdline}")
    proc = subprocess.Popen(cmdargs)
    proc.wait()
    status = proc.returncode
    if status != 0:
        log_error = extract_sphinx_error(errfile)
        raise SphinxError(cmdline, f"return code = {status}", log_error)


def extract_sphinx_error(errfile: str):
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


TERM_WIDTH = 60  # TODO: figure this out


def notify(message, level=0):
    if level == 0:
        print("*" * TERM_WIDTH)
        print(f"* {message}")
        print("*" * TERM_WIDTH)
    else:
        print(f"===> {message}")


def main():
    global g_force_rebuild

    ap = argparse.ArgumentParser()
    ap.add_argument("-N", "--no-notebooks", action="store_true")
    ap.add_argument("-S", "--no-sphinx", action="store_true")
    ap.add_argument("-c", "--config", default="build.yml")
    ap.add_argument("-k", "--kernel", default=None, dest="kernel")
    ap.add_argument(
        "-R",
        "--rebuild-all",
        default=False,
        action="store_true",
        dest="force_rebuild",
        help="Rebuild all notebooks",
    )
    ap.add_argument("-v", "--verbose", action="count", dest="vb", default=0)
    args = ap.parse_args()

    # further arg processing
    if args.vb > 1:
        _log.setLevel(logging.DEBUG)
    elif args.vb > 0:
        _log.setLevel(logging.INFO)
    if args.kernel is None:
        _log.warning("using default kernel for notebook execution")
        kernel = None
    else:
        _log.info(f"using kernel {args.kernel} for notebook execution")
        kernel = args.kernel

    g_force_rebuild = args.force_rebuild

    try:
        f = open(args.config)
        config = yaml.load(f, Loader=yaml.SafeLoader)
    except IOError as err:
        _log.error(f"Configuration file is required: {err}")
        return 1
    except yaml.error.YAMLError as err:
        _log.error(f"Cannot parse configuration file, expected YAML: {err}")
        return 1

    if not args.no_notebooks:
        try:
            notify("Convert Jupyter notebooks")
            build_notebooks(config, kernel=kernel)
        except NotebookError as err:
            _log.fatal(f"Could not build notebooks: {err}")
            return -1

    if not args.no_sphinx:
        try:
            notify("Sphinx documentation")
            build_sphinx(config)
        except SphinxError as err:
            _log.fatal(f"Could not build Sphinx docs: {err}")
            return -1

    return 0


if __name__ == "__main__":
    sys.exit(main())
