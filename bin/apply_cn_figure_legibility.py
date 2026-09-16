#!/usr/bin/env python3
"""
apply_cn_figure_legibility.py

Follow-up to item 5. The first ICHORCNA_PLOT grids were unreadable: the
ISCN string from ichorkaryo lists every arm at every copy number, including
the forty-odd normal ones, so the suptitle ran to several hundred characters
and bbox_inches="tight" widened the canvas around it until the 22 panels
occupied a third of the image.

Two changes:

  bin/ichorkaryo.py   iscn_string() lists only arms whose copy number differs
                      from their expected value (autosomes against the modal
                      baseline, chrX/chrY against sex), as an ISCN line
                      would. "seq(10p)x1 (6p, 8q, 15q)x3 (1q)x4" rather than
                      fifty labels. An unremarkable genome reads
                      "seq: no arm-level change". This is what the dashboard
                      headline shows, so it changes there too.
  bin/mmplot.py       the header is wrapped to the figure width before it is
                      drawn, in every layout; grid cells go from 5.4x4.4 in
                      to 6.4x4.8 in; grid dpi 110 -> 120.

ICHORCNA_PLOT itself is unchanged: it still composes ploidy_class +
canonical trisomies + ISCN, which is now short.

Idempotent; sentinel `cn_figure_legibility_v1`; --check; --revert.
"""

import argparse
import os
import shutil
import sys

SENTINEL = "cn_figure_legibility_v1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

EDITS = []

# ---------------------------------------------------------------------------
# ichorkaryo.py
# ---------------------------------------------------------------------------
EDITS.append(("bin/ichorkaryo.py", [
    (
        "def iscn_string(arms, sex):\n"
        "    \"\"\"Group arms by copy number into an ISCN-style seq line.\"\"\"\n"
        "    by_cn = defaultdict(list)\n"
        "    for label, value in arms.items():\n"
        "        by_cn[value[\"cn\"]].append(label)\n",
        "def iscn_string(arms, sex, baseline=2):\n"
        "    \"\"\"Group abnormal arms by copy number into an ISCN-style seq line.\n"
        "\n"
        "    cn_figure_legibility_v1: only arms whose copy number differs from\n"
        "    the expected value are listed, as an ISCN line would. Autosomes are\n"
        "    judged against the modal baseline, chrX/chrY against sex. Listing\n"
        "    every normal arm made the string several hundred characters long\n"
        "    and unusable as a figure title.\n"
        "    \"\"\"\n"
        "    by_cn = defaultdict(list)\n"
        "    for label, value in arms.items():\n"
        "        cn = value.get(\"cn\")\n"
        "        if cn is None:\n"
        "            continue\n"
        "        chrom = \"chr\" + label.rstrip(\"pq\")\n"
        "        expected = expected_copies_for(chrom, sex, baseline)\n"
        "        if expected is not None and cn == expected:\n"
        "            continue\n"
        "        by_cn[cn].append(label)\n"
        "    if not by_cn:\n"
        "        return \"seq: no arm-level change\"\n",
    ),
    (
        "        \"ISCN_karyotype\": iscn_string(arms, sex),\n",
        "        \"ISCN_karyotype\": iscn_string(arms, sex, baseline),\n",
    ),
]))

# ---------------------------------------------------------------------------
# mmplot.py
# ---------------------------------------------------------------------------
EDITS.append(("bin/mmplot.py", [
    (
        "import argparse\n"
        "import math\n"
        "import os\n"
        "import sys\n",
        "import argparse\n"
        "import math\n"
        "import os\n"
        "import sys\n"
        "import textwrap  # cn_figure_legibility_v1\n",
    ),
    (
        "    if args.subtitle:\n"
        "        header += \"   |   \" + args.subtitle\n",
        "    if args.subtitle:\n"
        "        header += \"   |   \" + args.subtitle\n"
        "\n"
        "    def wrap_header(text, figure_width_in, font_pt):\n"
        "        # cn_figure_legibility_v1. bbox_inches='tight' grows the canvas to\n"
        "        # fit the title, so an unwrapped long title shrinks every panel.\n"
        "        # Roughly 0.55 em per character at the given size.\n"
        "        chars = max(40, int(figure_width_in * 72.0 / (font_pt * 0.55)))\n"
        "        return \"\\n\".join(textwrap.wrap(text, chars)) or text\n",
    ),
    (
        "        built[0].set_title(\"%s  --  %s\" % (chrom, header), fontsize=13, pad=22)\n",
        "        built[0].set_title(wrap_header(\"%s  --  %s\" % (chrom, header), 17, 13),\n"
        "                           fontsize=13, pad=22)\n",
    ),
    (
        "        ax_depth.set_title(header, fontsize=13)\n",
        "        ax_depth.set_title(wrap_header(header, 21, 13), fontsize=13)\n",
    ),
    (
        "    figure = pyplot.figure(figsize=(5.4 * columns, 4.4 * rows))\n",
        "    figure = pyplot.figure(figsize=(6.4 * columns, 4.8 * rows))\n",
    ),
    (
        "    outer = gridspec.GridSpec(rows, columns, hspace=0.48, wspace=0.22)\n",
        "    # cn_figure_legibility_v1: explicit margins so the panels use the canvas\n"
        "    # rather than leaving a fifth of it blank above and below the grid.\n"
        "    outer = gridspec.GridSpec(rows, columns, hspace=0.48, wspace=0.22,\n"
        "                              left=0.05, right=0.99, top=0.955, bottom=0.04)\n",
    ),
    (
        "    figure.suptitle(header, fontsize=15, y=0.998)\n",
        "    figure.suptitle(wrap_header(header, 6.4 * columns, 15), fontsize=15, y=0.992)\n",
    ),
    (
        "    figure.savefig(args.output, bbox_inches=\"tight\", dpi=110)\n",
        "    figure.savefig(args.output, bbox_inches=\"tight\", dpi=120)\n",
    ),
]))


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def check():
    problems = 0
    for rel, pairs in EDITS:
        path = os.path.join(REPO, rel)
        if not os.path.isfile(path):
            print(f"MISSING  {rel}"); problems += 1; continue
        text = read(path)
        if SENTINEL in text:
            print(f"APPLIED  {rel}"); continue
        bad = [old for old, _ in pairs if text.count(old) != 1]
        if bad:
            print(f"NO-MATCH {rel}: {len(bad)} anchor(s) not found exactly once")
            for old in bad:
                print(f"         anchor starts: {old.strip().splitlines()[0]!r}")
            problems += 1
        else:
            print(f"READY    {rel}")
    return 1 if problems else 0


def apply():
    staged = []
    for rel, pairs in EDITS:
        path = os.path.join(REPO, rel)
        if not os.path.isfile(path):
            sys.exit(f"ERROR: {rel} not found; nothing changed.")
        text = read(path)
        if SENTINEL in text:
            print(f"skip     {rel} (already applied)"); continue
        for old, new in pairs:
            n = text.count(old)
            if n != 1:
                sys.exit(f"ERROR: anchor found {n} times in {rel} (expected 1); nothing changed.\n"
                         f"       anchor starts: {old.strip().splitlines()[0]!r}")
            text = text.replace(old, new, 1)
        if SENTINEL not in text:
            sys.exit(f"ERROR: internal: sentinel absent after edit of {rel}")
        staged.append((rel, path, text))
    for rel, path, text in staged:
        shutil.copy2(path, path + ".bak")
        write(path, text)
        print(f"patched  {rel} (backup: {rel}.bak)")


def revert():
    for rel, _ in EDITS:
        path = os.path.join(REPO, rel)
        if os.path.isfile(path + ".bak"):
            shutil.move(path + ".bak", path); print(f"restored {rel}")
        else:
            print(f"no .bak  {rel}")


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true")
    g.add_argument("--revert", action="store_true")
    args = ap.parse_args()
    if args.check:
        sys.exit(check())
    if args.revert:
        revert(); return
    apply()


if __name__ == "__main__":
    main()
