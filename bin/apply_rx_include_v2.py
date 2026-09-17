#!/usr/bin/env python3
"""
apply_rx_include_v2.py  (rx_include_v2)

Two fixes on top of rx_include_v1.

1. Rearrangement Include boxes ticked but nothing reached the Reporting tab.
   The Rearrangements tab's inline script runs where it sits in the body,
   before variant-browser.js (loaded at the end of body) defines the store,
   so the wiring found no store and returned; the boxes were toggling as
   bare DOM checkboxes. The wiring now waits for window load, which is also
   after $(document).ready has set the store's sample.

2. Include for the panel-window rows on the Copy number & allelic state tab.
   The genome and per-chromosome figures have an Include card control; the
   window rows (TP53+TNFSF12 with its depth + allele-fraction plot, and the
   depth-only gene windows) had none. Each window row gets an Include box;
   the Reporting tab gets an "Allelic state at panel windows" section that
   shows the row's numbers, the inline plot cloned from the tab where one
   exists, and an editable interpretation. Included in the TSV and Clear all.

Files: templates/rearrangements_tab.html.j2, templates/copy_number_tab.html.j2,
templates/sample_report.html.j2. Requires rx_include_v1 in the first and
last. Backups are <file>.rx2.bak.

Usage (from the repo root):
  python3 bin/apply_rx_include_v2.py            apply
  python3 bin/apply_rx_include_v2.py --check
  python3 bin/apply_rx_include_v2.py --revert
"""

import argparse
import shutil
import sys
from pathlib import Path

SENTINEL = "rx_include_v2"
BAK_SUFFIX = ".rx2.bak"
DASH = "bin/dashboard_builder/"

EDITS = [
    (
        DASH + "templates/rearrangements_tab.html.j2",
        [
            (
                "    // rx_include_v1: Include checkboxes, backed by the Reporting store.\n"
                "    (function () {\n"
                "      var store = window.tspipeReporting;\n"
                "      if (!store) { return; }\n",
                "    // rx_include_v1: Include checkboxes, backed by the Reporting store.\n"
                "    // rx_include_v2: this inline script runs before variant-browser.js\n"
                "    // (loaded at the end of body) defines the store, and before\n"
                "    // $(document).ready sets its sample, so the wiring waits for load.\n"
                "    window.addEventListener('load', function () {\n"
                "      var store = window.tspipeReporting;\n"
                "      if (!store) { return; }\n",
            ),
            (
                "      store.onChange(sync);\n"
                "      sync();\n"
                "    })();\n",
                "      store.onChange(sync);\n"
                "      sync();\n"
                "    });\n",
            ),
        ],
    ),
    (
        DASH + "templates/copy_number_tab.html.j2",
        [
            (
                "                  <th>Allelic state</th><th class=\"text-end\">Hets / usable</th><th>Reading</th>\n"
                "                </tr></thead>\n",
                "                  <th>Allelic state</th><th class=\"text-end\">Hets / usable</th><th>Reading</th>\n"
                "                  <th title=\"tick to include this window in the Reporting tab\">Include</th>{# rx_include_v2 #}\n"
                "                </tr></thead>\n",
            ),
            (
                "                    <td>{{ w.reading }}</td>\n"
                "                  </tr>\n"
                "                  {% if w.svg %}\n"
                "                    <tr class=\"cn-win-plot\" data-cn-win-plot=\"{{ page.chrom }}-{{ loop.index }}\" style=\"display:none\">\n"
                "                      <td colspan=\"9\" class=\"bg-white\">{{ w.svg | safe }}\n",
                "                    <td>{{ w.reading }}</td>\n"
                "                    <td class=\"text-center\">\n"
                "                      <input type=\"checkbox\" class=\"form-check-input cn-win-include\" title=\"include in the Reporting tab\"\n"
                "                             data-win-id=\"{{ w.chrom }}:{{ w.gene }}\" data-gene=\"{{ w.gene }}\" data-chrom=\"{{ w.chrom }}\" data-band=\"{{ w.band }}\"\n"
                "                             data-log2=\"{{ '%.2f'|format(w.log2_ontarget) if w.log2_ontarget is not none else '' }}\"\n"
                "                             data-cn=\"{{ w.cn if w.cn is not none else '' }}{% if w.cn_event and w.cn_event != 'none' %} ({{ w.cn_event }}){% endif %}\"\n"
                "                             data-cf=\"{{ '%.0f%%'|format(w.cell_fraction*100) if w.cell_fraction is not none else '' }}\"\n"
                "                             data-state=\"{{ w.allelic_state or 'not assessed' }}\" data-hets=\"{% if w.n_usable %}{{ w.n_het }} / {{ w.n_usable }}{% endif %}\"\n"
                "                             data-reading=\"{{ w.reading }}\" data-svg-id=\"{% if w.svg %}cn-win-svg-{{ page.chrom }}-{{ loop.index }}{% endif %}\">\n"
                "                    </td>\n"
                "                  </tr>\n"
                "                  {% if w.svg %}\n"
                "                    <tr class=\"cn-win-plot\" data-cn-win-plot=\"{{ page.chrom }}-{{ loop.index }}\" style=\"display:none\">\n"
                "                      <td colspan=\"10\" class=\"bg-white\"><div id=\"cn-win-svg-{{ page.chrom }}-{{ loop.index }}\">{{ w.svg | safe }}</div>\n",
            ),
            (
                "  <script>\n"
                "  (function () {\n"
                "    document.addEventListener(\"click\", function (event) {\n"
                "      var row = event.target.closest ? event.target.closest(\".cn-win-row\") : null;\n"
                "      if (!row) { return; }\n",
                "  <script>\n"
                "  (function () {\n"
                "    // rx_include_v2: window Include boxes, backed by the Reporting store.\n"
                "    // Waits for load: the store is defined by a script at the end of body.\n"
                "    window.addEventListener(\"load\", function () {\n"
                "      var store = window.tspipeReporting;\n"
                "      if (!store) { return; }\n"
                "      function snapshotOf(box) {\n"
                "        var s = {};\n"
                "        [\"gene\", \"chrom\", \"band\", \"log2\", \"cn\", \"cf\", \"state\", \"hets\", \"reading\"].forEach(function (k) {\n"
                "          s[k] = box.getAttribute(\"data-\" + k) || \"\";\n"
                "        });\n"
                "        s.svg_id = box.getAttribute(\"data-svg-id\") || \"\";\n"
                "        return s;\n"
                "      }\n"
                "      function sync() {\n"
                "        var ids = {};\n"
                "        store.load().forEach(function (it) { if (it.kind === \"allelic_window\") { ids[it.id] = true; } });\n"
                "        document.querySelectorAll(\".cn-win-include\").forEach(function (box) {\n"
                "          var on = !!ids[box.getAttribute(\"data-win-id\")];\n"
                "          if (box.checked !== on) { box.checked = on; }\n"
                "        });\n"
                "      }\n"
                "      document.addEventListener(\"change\", function (e) {\n"
                "        var box = e.target.closest ? e.target.closest(\".cn-win-include\") : null;\n"
                "        if (!box) { return; }\n"
                "        store.toggle(\"allelic_window\", box.getAttribute(\"data-win-id\"), snapshotOf(box));\n"
                "      });\n"
                "      store.onChange(sync);\n"
                "      sync();\n"
                "    });\n"
                "\n"
                "    document.addEventListener(\"click\", function (event) {\n"
                "      if (event.target.closest && event.target.closest(\"input\")) { return; }\n"
                "      var row = event.target.closest ? event.target.closest(\".cn-win-row\") : null;\n"
                "      if (!row) { return; }\n",
            ),
        ],
    ),
    (
        DASH + "templates/sample_report.html.j2",
        [
            (
                "              Variants ticked in the Clinical tab and CNV plots ticked in the CNV\n"
                "              tab appear here in selection order.",
                "              Rearrangement loci, panel windows, variants and copy-number plots\n"
                "              ticked elsewhere in this report appear here in selection order.",
            ),
            (
                "        {# ----- Variants section ----- #}\n"
                "        <div id=\"reporting-variants-section\" style=\"display:none\">\n",
                "        {# ----- Panel windows section (rx_include_v2) ----- #}\n"
                "        <div id=\"reporting-win-section\" class=\"mt-4\" style=\"display:none\">\n"
                "          <h5 class=\"mt-2 mb-2\">\n"
                "            Copy number and allelic state at panel windows\n"
                "            <span id=\"reporting-win-count-badge\" class=\"badge bg-secondary ms-1\">0</span>\n"
                "          </h5>\n"
                "          <p class=\"text-muted small mb-2\">\n"
                "            Windows ticked on the Copy number &amp; allelic state tab. Where the\n"
                "            window was assessed, its depth and allele-fraction plot is shown.\n"
                "            The interpretation is editable and persists in this browser.\n"
                "          </p>\n"
                "          <div id=\"reporting-win-list\"></div>\n"
                "        </div>\n"
                "\n"
                "        {# ----- Variants section ----- #}\n"
                "        <div id=\"reporting-variants-section\" style=\"display:none\">\n",
            ),
            (
                "    const rxCountBadge = document.getElementById(\"reporting-rx-count-badge\");\n",
                "    const rxCountBadge = document.getElementById(\"reporting-rx-count-badge\");\n"
                "    const winSection = document.getElementById(\"reporting-win-section\");    // rx_include_v2\n"
                "    const winList = document.getElementById(\"reporting-win-list\");\n"
                "    const winCountBadge = document.getElementById(\"reporting-win-count-badge\");\n",
            ),
            (
                "    function renderVariants(variants) {\n",
                "    // rx_include_v2. Notes share the rx-note store; window ids (chrN:GENE)\n"
                "    // and locus keys (A::B|...) cannot collide.\n"
                "    const WIN_COLUMNS = [\"Window\", \"Chromosome\", \"Band\", \"Depth log2 (on-target)\", \"Copy number\",\n"
                "                         \"Cell fraction\", \"Allelic state\", \"Hets / usable\", \"Reading\", \"Interpretation\"];\n"
                "\n"
                "    function winRowToTsv(it) {\n"
                "      const s = it.snapshot || {};\n"
                "      const note = (window.tspipeReporting.getRxNote(it.id) || \"\").replace(/\\s+/g, \" \");\n"
                "      return [s.gene || \"\", s.chrom || \"\", s.band || \"\", s.log2 || \"\", s.cn || \"\", s.cf || \"\",\n"
                "              s.state || \"\", s.hets || \"\", s.reading || \"\", note].join(\"\\t\");\n"
                "    }\n"
                "\n"
                "    function renderWins(wins) {\n"
                "      if (!winSection) return;\n"
                "      if (wins.length === 0) {\n"
                "        winSection.style.display = \"none\";\n"
                "        winList.innerHTML = \"\";\n"
                "        return;\n"
                "      }\n"
                "      winSection.style.display = \"\";\n"
                "      winCountBadge.textContent = String(wins.length);\n"
                "      winList.innerHTML = wins.map(function (it, idx) {\n"
                "        const s = it.snapshot || {};\n"
                "        const note = escapeHtml(window.tspipeReporting.getRxNote(it.id) || \"\");\n"
                "        const facts = [\n"
                "          [\"Depth log2\", s.log2], [\"Copy number\", s.cn], [\"Cell fraction\", s.cf],\n"
                "          [\"Allelic state\", s.state], [\"Hets / usable\", s.hets], [\"Reading\", s.reading]\n"
                "        ].filter(function (p) { return p[1]; }).map(function (p) {\n"
                "          return '<span class=\"me-3\"><span class=\"text-muted\">' + p[0] + '</span> ' + escapeHtml(p[1]) + '</span>';\n"
                "        }).join(\"\");\n"
                "        return '<div class=\"card mb-3\" data-win-id=\"' + escapeHtml(it.id) + '\">' +\n"
                "                 '<div class=\"card-header d-flex justify-content-between align-items-center py-2\">' +\n"
                "                   '<span class=\"fw-medium\">' + (idx + 1) + '. ' + escapeHtml(s.gene) +\n"
                "                     ' <span class=\"text-muted small\">' + escapeHtml(s.chrom) + ' ' + escapeHtml(s.band) + '</span></span>' +\n"
                "                   '<button type=\"button\" class=\"btn btn-sm btn-outline-secondary reporting-win-remove\" ' +\n"
                "                     'data-win-id=\"' + escapeHtml(it.id) + '\" title=\"Remove from report\">Remove</button>' +\n"
                "                 '</div>' +\n"
                "                 '<div class=\"card-body\">' +\n"
                "                   '<div class=\"small mb-2\">' + facts + '</div>' +\n"
                "                   (s.svg_id ? '<div data-win-svg-slot=\"' + escapeHtml(s.svg_id) + '\"></div>' : '') +\n"
                "                   '<label class=\"form-label mt-2 mb-1 small text-muted\">Interpretation (editable, saved automatically)</label>' +\n"
                "                   '<textarea class=\"form-control form-control-sm reporting-win-note-input\" data-win-id=\"' + escapeHtml(it.id) + '\" rows=\"2\" ' +\n"
                "                     'placeholder=\"e.g. 17p13.1 balanced at 173 heterozygous sites; no evidence of TP53 loss\">' + note + '</textarea>' +\n"
                "                 '</div>' +\n"
                "               '</div>';\n"
                "      }).join(\"\");\n"
                "      // The plot lives on the copy-number tab of this same page; clone it\n"
                "      // rather than store 100 KB of SVG per window in localStorage.\n"
                "      winList.querySelectorAll(\"[data-win-svg-slot]\").forEach(function (slot) {\n"
                "        const src = document.getElementById(slot.getAttribute(\"data-win-svg-slot\"));\n"
                "        if (src) slot.appendChild(src.cloneNode(true));\n"
                "      });\n"
                "    }\n"
                "\n"
                "    function renderVariants(variants) {\n",
            ),
            (
                "      const rxs = all.filter(function (it) { return it.kind === \"rearrangement\"; });   // rx_include_v1\n"
                "      const total = variants.length + excluded.length + cnvs.length + rxs.length;\n",
                "      const rxs = all.filter(function (it) { return it.kind === \"rearrangement\"; });   // rx_include_v1\n"
                "      const wins = all.filter(function (it) { return it.kind === \"allelic_window\"; }); // rx_include_v2\n"
                "      const total = variants.length + excluded.length + cnvs.length + rxs.length + wins.length;\n",
            ),
            (
                "        if (rxSection) rxSection.style.display = \"none\";\n",
                "        if (rxSection) rxSection.style.display = \"none\";\n"
                "        if (winSection) winSection.style.display = \"none\";\n",
            ),
            (
                "      renderRx(rxs);\n"
                "      renderVariants(variants);\n",
                "      renderRx(rxs);\n"
                "      renderWins(wins);\n"
                "      renderVariants(variants);\n",
            ),
            (
                "        rxs.forEach(function (it) { sections.push(rxRowToTsv(it)); });\n"
                "      }\n",
                "        rxs.forEach(function (it) { sections.push(rxRowToTsv(it)); });\n"
                "      }\n"
                "      const wins = all.filter(function (it) { return it.kind === \"allelic_window\"; });\n"
                "      if (wins.length) {\n"
                "        if (sections.length) sections.push(\"\");\n"
                "        sections.push(\"# Copy number and allelic state at panel windows\");\n"
                "        sections.push(WIN_COLUMNS.join(\"\\t\"));\n"
                "        wins.forEach(function (it) { sections.push(winRowToTsv(it)); });\n"
                "      }\n",
            ),
            (
                "          .forEach(function (it) { window.tspipeReporting.setRxNote(it.id, \"\"); });\n"
                "        window.tspipeReporting.clearAll();\n",
                "          .forEach(function (it) { window.tspipeReporting.setRxNote(it.id, \"\"); });\n"
                "        window.tspipeReporting.load()\n"
                "          .filter(function (it) { return it.kind === \"allelic_window\"; })\n"
                "          .forEach(function (it) { window.tspipeReporting.setRxNote(it.id, \"\"); });\n"
                "        window.tspipeReporting.clearAll();\n",
            ),
            (
                "    // Re-render on selection change\n"
                "    window.tspipeReporting.onChange(render);\n",
                "    // rx_include_v2: window notes persist on input; remove toggles off.\n"
                "    if (winList) {\n"
                "      winList.addEventListener(\"input\", function (ev) {\n"
                "        const ta = ev.target.closest(\".reporting-win-note-input\");\n"
                "        if (!ta) return;\n"
                "        window.tspipeReporting.setRxNote(ta.getAttribute(\"data-win-id\"), ta.value);\n"
                "      });\n"
                "      winList.addEventListener(\"click\", function (ev) {\n"
                "        const btn = ev.target.closest(\".reporting-win-remove\");\n"
                "        if (!btn) return;\n"
                "        const key = btn.getAttribute(\"data-win-id\");\n"
                "        const item = window.tspipeReporting.load().find(function (it) {\n"
                "          return it.kind === \"allelic_window\" && it.id === key;\n"
                "        });\n"
                "        if (item) {\n"
                "          window.tspipeReporting.setRxNote(key, \"\");\n"
                "          window.tspipeReporting.toggle(\"allelic_window\", key, item.snapshot);\n"
                "        }\n"
                "      });\n"
                "    }\n"
                "\n"
                "    // Re-render on selection change\n"
                "    window.tspipeReporting.onChange(render);\n",
            ),
        ],
    ),
]

PREREQ = {DASH + "templates/rearrangements_tab.html.j2": "rx_include_v1",
          DASH + "templates/sample_report.html.j2": "rx_include_v1"}


def state(root):
    out = []
    for rel, edits in EDITS:
        p = root / rel
        if not p.is_file():
            out.append((rel, "MISSING FILE"))
            continue
        text = p.read_text()
        if SENTINEL in text:
            out.append((rel, "applied"))
            continue
        if rel in PREREQ and PREREQ[rel] not in text:
            out.append((rel, "PREREQUISITE MISSING: %s" % PREREQ[rel]))
            continue
        bad = [a for a, _ in edits if text.count(a) != 1]
        if bad:
            out.append((rel, "ANCHOR MISMATCH: %d of %d anchors not unique" % (len(bad), len(edits))))
            for a in bad:
                sys.stderr.write("    count=%d  %r\n" % (text.count(a), a.splitlines()[0]))
            continue
        out.append((rel, "ready"))
    return out


def apply(root):
    st = state(root)
    if any(s not in ("ready", "applied") for _, s in st):
        for rel, s in st:
            print("  %-52s %s" % (rel, s))
        sys.exit("Refusing to apply: fix the problems above first.")
    if all(s == "applied" for _, s in st):
        print("Already applied (sentinel present in every file); nothing to do.")
        return
    for rel, edits in EDITS:
        p = root / rel
        text = p.read_text()
        if SENTINEL in text:
            print("  %-52s already applied, skipped" % rel)
            continue
        bak = p.with_name(p.name + BAK_SUFFIX)
        shutil.copy2(p, bak)
        for anchor, repl in edits:
            text = text.replace(anchor, repl, 1)
        assert SENTINEL in text, "sentinel absent after patching %s" % rel
        p.write_text(text)
        print("  %-52s patched (backup: %s)" % (rel, bak.name))


def revert(root):
    for rel, _ in EDITS:
        p = root / rel
        bak = p.with_name(p.name + BAK_SUFFIX)
        if not bak.is_file():
            print("  %-52s no %s, skipped" % (rel, BAK_SUFFIX))
            continue
        shutil.copy2(bak, p)
        bak.unlink()
        print("  %-52s restored from %s" % (rel, BAK_SUFFIX))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1].strip())
    ap.add_argument("--root", default=".")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true")
    g.add_argument("--revert", action="store_true")
    args = ap.parse_args()
    root = Path(args.root).resolve()
    if args.check:
        for rel, s in state(root):
            print("  %-52s %s" % (rel, s))
        return
    if args.revert:
        revert(root)
        return
    apply(root)


if __name__ == "__main__":
    main()
