#!/usr/bin/env python3
"""
apply_report_rebuild.py

Wire the rebuilt report tabs into the existing dashboard builder.

What it changes
---------------
bin/dashboard_builder/build.py
    - imports the three new parsers (translocations, ichor, qc)
    - populates ctx["translocations"], ctx["ichor"], ctx["qc"]
    The existing parser calls are left untouched; this patch is additive on
    the Python side, so a rollback is a template-only concern.

bin/dashboard_builder/templates/sample_report.html.j2
    - adds a Translocations tab
    - relabels the CNV tab to "Copy number" and replaces its pane with the
      ichorCNA include (the segment call table goes)
    - replaces the QC pane with the adaptive-sampling QC include
    - removes the FLT3 tab and pane, which belong to the reference pipeline
      this builder was ported from and have no input here

Safety
------
- Writes a timestamped .bak of every file it edits.
- Refuses to apply twice: each edit carries a sentinel that is checked first.
- Validates every anchor before writing anything. If any anchor is missing the
  script reports which one and exits without touching a file, so a template
  that has drifted fails loudly instead of being half-patched.
- --dry-run prints the planned edits, including the full text of any block it
  would delete, and writes nothing.

Usage
-----
    python3 bin/apply_report_rebuild.py --dry-run
    python3 bin/apply_report_rebuild.py
    python3 bin/apply_report_rebuild.py --repo /goast/mm-awgs-nextflow
"""

import argparse
import os
import re
import shutil
import sys
from datetime import datetime


SENTINEL_BUILD_IMPORT = "from parsers import translocations as p_translocations"
SENTINEL_BUILD_CTX = 'ctx["translocations"]'
SENTINEL_TPL_TX = "translocations_tab.html.j2"
SENTINEL_TPL_QC = "qc_tab.html.j2"
SENTINEL_TPL_ICHOR = "ichor_tab.html.j2"
SENTINEL_TPL_OVERVIEW = "overview_tab.html.j2"

NEW_IMPORTS = (
    "from parsers import translocations as p_translocations\n"
    "from parsers import ichor as p_ichor\n"
    "from parsers import qc as p_qc\n"
)


class PatchError(Exception):
    """Raised when an anchor cannot be located; aborts before any write."""


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def backup_and_write(path, text, dry_run):
    if dry_run:
        return
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.copy2(path, "%s.bak_report_rebuild_%s" % (path, stamp))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def find_line(text, needle, description):
    """Return (start, end) of the first line containing needle."""
    index = text.find(needle)
    if index == -1:
        raise PatchError("anchor not found: %s (%r)" % (description, needle))
    start = text.rfind("\n", 0, index) + 1
    end = text.find("\n", index)
    if end == -1:
        end = len(text)
    return start, end + 1


def indent_of(line):
    return line[: len(line) - len(line.lstrip())]


def find_element(text, anchor, tag, description):
    """Return (start, end) of the non-nesting element containing an anchor.

    Used for navigation list items, which never nest. Operating on the element
    rather than on the whole line means an edit stays surgical even if the
    template puts several items on one line.
    """
    index = text.find(anchor)
    if index == -1:
        raise PatchError("anchor not found: %s (%r)" % (description, anchor))
    start = text.rfind("<%s" % tag, 0, index)
    if start == -1:
        raise PatchError("no opening <%s> before anchor: %s" % (tag, description))
    close = "</%s>" % tag
    end = text.find(close, index)
    if end == -1:
        raise PatchError("no closing %s after anchor: %s" % (close, description))
    return start, end + len(close)


def find_div_block(text, anchor, description):
    """Return (start, end) of a <div ...> block located by a unique anchor.

    Scans forward from the anchor's opening tag, counting nested <div and
    </div> until the depth returns to zero. Bootstrap tab panes are well
    formed, so this is exact; if the depth never closes the caller is told
    rather than the file being truncated.
    """
    index = text.find(anchor)
    if index == -1:
        raise PatchError("anchor not found: %s (%r)" % (description, anchor))
    start = text.rfind("<div", 0, index)
    if start == -1:
        raise PatchError("no opening <div> before anchor: %s" % description)
    # Extend back over the opening tag's own indentation so the replacement
    # controls the whole line rather than inheriting stale leading whitespace.
    line_start = text.rfind("\n", 0, start) + 1
    if not text[line_start:start].strip():
        start = line_start

    depth = 0
    position = start
    pattern = re.compile(r"<div\b|</div\s*>", re.IGNORECASE)
    while True:
        match = pattern.search(text, position)
        if not match:
            raise PatchError("unbalanced <div> after anchor: %s" % description)
        if match.group(0).lower().startswith("<div"):
            depth += 1
        else:
            depth -= 1
            if depth == 0:
                end = match.end()
                # Consume the trailing newline so removal does not leave a gap.
                if text[end : end + 1] == "\n":
                    end += 1
                return start, end
        position = match.end()


def find_function_return(text, def_line, description="function return"):
    """Return (offset, indent) of a function's top-level return statement.

    Scans the body of the named function for a ``return`` at exactly one
    indentation level deeper than the ``def``, which is the statement that
    hands the context dictionary back. Insertion at that point is guaranteed
    to run on every call, unlike insertion next to a parser call that happens
    to sit inside a conditional.
    """
    index = text.find(def_line)
    if index == -1:
        raise PatchError("function not found: %s" % def_line)

    def_start = text.rfind("\n", 0, index) + 1
    def_indent = indent_of(text[def_start:])
    body_indent = def_indent + "    "

    # Walk the file line by line from the def, tracking byte offsets. The
    # function ends at the first non-blank line indented no deeper than the
    # def itself; bounding the search this way stops a return belonging to a
    # later function from being selected.
    offset = def_start
    in_body = False
    last_return = None
    for line in text[def_start:].splitlines(keepends=True):
        stripped = line.strip()
        if in_body and stripped and len(indent_of(line)) <= len(def_indent):
            break
        if in_body and stripped.startswith("return") and (
            indent_of(line) == body_indent
        ):
            last_return = offset
        if not in_body and stripped.startswith("def "):
            in_body = True
        offset += len(line)

    if last_return is None:
        raise PatchError(
            "no top-level return found in %s; cannot place the context "
            "assignments safely" % def_line
        )
    return last_return, body_indent


# ---------------------------------------------------------------------------
# build.py
# ---------------------------------------------------------------------------

def patch_build(path, dry_run):
    text = read(path)
    actions = []

    if SENTINEL_BUILD_IMPORT in text:
        actions.append("SKIP  build.py imports: already applied")
    else:
        start, end = find_line(
            text, "from parsers import baf_loh", "last parser import"
        )
        text = text[:end] + NEW_IMPORTS + text[end:]
        actions.append("EDIT  build.py: add translocations/ichor/qc imports")

    if SENTINEL_BUILD_CTX in text:
        actions.append("SKIP  build.py context: already applied")
    else:
        # The existing parser calls sit inside conditionals, so anchoring on
        # one of them would nest these three inside that same guard and make
        # the translocation, copy-number and QC tabs depend on an unrelated
        # input being present. Insert at the function's return instead, which
        # is unconditional and at a known indentation.
        start, pad = find_function_return(
            text, "def collect_sample_context"
        )
        block = (
            "%s# Structural variants, copy number and adaptive-sampling QC.\n"
            "%s# Each parser returns None when its inputs are absent, and the\n"
            "%s# matching template renders an empty state rather than failing.\n"
            "%sctx[\"translocations\"] = p_translocations.parse(effective_dir, sample)\n"
            "%sctx[\"ichor\"] = p_ichor.parse(effective_dir, sample)\n"
            "%sctx[\"qc\"] = p_qc.parse(effective_dir, sample)\n"
            % (pad, pad, pad, pad, pad, pad)
        )
        text = text[:start] + block + text[start:]
        actions.append(
            "EDIT  build.py: add ctx translocations/ichor/qc before the "
            "return of collect_sample_context (indent %d)" % len(pad)
        )

    return text, actions


# ---------------------------------------------------------------------------
# sample_report.html.j2
# ---------------------------------------------------------------------------

def patch_template(path, dry_run):
    text = read(path)
    actions = []
    removals = []

    # --- navigation ---------------------------------------------------------
    if 'data-bs-target="#tab-translocations"' in text:
        actions.append("SKIP  template nav: translocations item already present")
    else:
        start, end = find_element(
            text, 'data-bs-target="#tab-cnv"', "li", "CNV nav item"
        )
        cnv_item = text[start:end]
        line_start = text.rfind("\n", 0, start) + 1
        pad = text[line_start:start] if not text[line_start:start].strip() else ""
        nav = (
            '<li class="nav-item"><button class="nav-link"        '
            'data-bs-toggle="pill" data-bs-target="#tab-translocations" '
            'type="button" role="tab">Translocations</button></li>\n%s' % pad
        )
        relabelled = cnv_item.replace(">CNV<", ">Copy number<")
        text = text[:start] + nav + relabelled + text[end:]
        actions.append(
            "EDIT  template nav: add Translocations, relabel CNV to Copy number"
        )

    # FLT3 belongs to the pipeline this builder was ported from.
    if 'data-bs-target="#tab-flt3"' in text:
        start, end = find_element(
            text, 'data-bs-target="#tab-flt3"', "li", "FLT3 nav item"
        )
        removals.append(("FLT3 nav item", text[start:end]))
        # Swallow the item's own indentation and trailing newline so the
        # removal does not leave a blank line behind.
        line_start = text.rfind("\n", 0, start) + 1
        if not text[line_start:start].strip():
            start = line_start
        if text[end : end + 1] == "\n":
            end += 1
        text = text[:start] + text[end:]
        actions.append("EDIT  template nav: remove FLT3")
    else:
        actions.append("SKIP  template nav: FLT3 already removed")

    # --- panes --------------------------------------------------------------
    if SENTINEL_TPL_OVERVIEW in text:
        actions.append("SKIP  template Overview pane: already replaced")
    else:
        start, end = find_div_block(text, 'id="tab-overview"', "Overview pane")
        removals.append(
            ("Overview pane (%d chars)" % (end - start), text[start:end])
        )
        text = (
            text[:start]
            + "      {% include 'overview_tab.html.j2' %}\n"
            + text[end:]
        )
        actions.append(
            "EDIT  template: replace Overview pane with overview_tab include"
        )

    if SENTINEL_TPL_QC in text:
        actions.append("SKIP  template QC pane: already replaced")
    else:
        start, end = find_div_block(text, 'id="tab-qc"', "QC pane")
        removals.append(("QC pane (%d chars)" % (end - start), text[start:end]))
        text = (
            text[:start]
            + "      {% include 'qc_tab.html.j2' %}\n"
            + text[end:]
        )
        actions.append("EDIT  template: replace QC pane with qc_tab include")

    if 'id="tab-flt3"' in text:
        start, end = find_div_block(text, 'id="tab-flt3"', "FLT3 pane")
        removals.append(("FLT3 pane (%d chars)" % (end - start), text[start:end]))
        text = text[:start] + text[end:]
        actions.append("EDIT  template: remove FLT3 pane")
    else:
        actions.append("SKIP  template: FLT3 pane already removed")

    if SENTINEL_TPL_ICHOR in text:
        actions.append("SKIP  template CNV pane: already replaced")
    else:
        start, end = find_div_block(text, 'id="tab-cnv"', "CNV pane")
        removals.append(("CNV pane (%d chars)" % (end - start), text[start:end]))
        replacement = ""
        if SENTINEL_TPL_TX not in text:
            replacement += "      {% include 'translocations_tab.html.j2' %}\n"
            actions.append("EDIT  template: add translocations_tab include")
        replacement += "      {% include 'ichor_tab.html.j2' %}\n"
        text = text[:start] + replacement + text[end:]
        actions.append(
            "EDIT  template: replace CNV pane with ichor_tab include"
        )

    return text, actions, removals


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        default=".",
        help="pipeline repository root (default: current directory)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="report planned edits and removed blocks; write nothing",
    )
    args = parser.parse_args(argv)

    builder = os.path.join(args.repo, "bin", "dashboard_builder")
    build_py = os.path.join(builder, "build.py")
    template = os.path.join(builder, "templates", "sample_report.html.j2")

    for path in (build_py, template):
        if not os.path.isfile(path):
            print("ERROR: not found: %s" % path, file=sys.stderr)
            return 2

    required = [
        os.path.join(builder, "parsers", "translocations.py"),
        os.path.join(builder, "parsers", "ichor.py"),
        os.path.join(builder, "parsers", "qc.py"),
        os.path.join(builder, "templates", "translocations_tab.html.j2"),
        os.path.join(builder, "templates", "ichor_tab.html.j2"),
        os.path.join(builder, "templates", "qc_tab.html.j2"),
        os.path.join(builder, "templates", "overview_tab.html.j2"),
    ]
    missing = [p for p in required if not os.path.isfile(p)]
    if missing:
        print(
            "ERROR: install the new parsers and templates first. Missing:",
            file=sys.stderr,
        )
        for path in missing:
            print("  " + path, file=sys.stderr)
        return 2

    try:
        build_text, build_actions = patch_build(build_py, args.dry_run)
        tpl_text, tpl_actions, removals = patch_template(template, args.dry_run)
    except PatchError as error:
        print("ERROR: %s" % error, file=sys.stderr)
        print(
            "Nothing was written. The template has drifted from the expected "
            "structure; inspect it and adjust the anchor.",
            file=sys.stderr,
        )
        return 2

    for action in build_actions + tpl_actions:
        print(action)

    if removals:
        print("\nBlocks that will be removed:")
        for label, block in removals:
            print("\n--- %s ---" % label)
            lines = block.splitlines()
            preview = lines if args.dry_run else lines[:6]
            for line in preview:
                print("  " + line)
            if not args.dry_run and len(lines) > 6:
                print("  ... (%d more lines)" % (len(lines) - 6))

    edits = [a for a in build_actions + tpl_actions if a.startswith("EDIT")]
    if not edits:
        print("\nNothing to do; all edits already applied.")
        return 0

    if args.dry_run:
        print("\nDry run: %d edit(s) planned, nothing written." % len(edits))
        return 0

    backup_and_write(build_py, build_text, args.dry_run)
    backup_and_write(template, tpl_text, args.dry_run)

    print("\nApplied %d edit(s). Backups written alongside each file." % len(edits))
    print("\nNext:")
    print("  python3 bin/dashboard_builder/build.py <bundle_dir>")
    return 0


if __name__ == "__main__":
    sys.exit(main())
