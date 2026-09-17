#!/usr/bin/env python3
"""
apply_rx_include.py  (rx_include_v1)

Include checkbox on rearrangement loci, and a Rearrangements section on the
Reporting tab.

Why: the Reporting tab already curates variants (Incl/Excl with tier and
reason) and CNV plots (Include with an interpretation), persisted per sample
in the browser and exported as TSV. Rearrangements were the one finding class
with no path into it. This adds one, using the same store, so the tab badge,
"Clear all" and the TSV export cover all three.

What changes
------------
  parsers/rearrangements.py       each locus carries a durable `key`
                                  (anchor::partner plus both breakend
                                  positions to 100 kb) that survives rebuilds;
                                  table row ids are positional and do not
  templates/rearrangements_tab    an Include column on every locus row
                                  (reportable and other), an "Include all
                                  reportable" button, checkbox state synced
                                  with the store
  assets/js/variant-browser.js    getRxNote/setRxNote: the editable
                                  interpretation per included locus, same
                                  persistence model as CNV captions
  templates/sample_report         "Rearrangements" section on the Reporting
                                  tab (anchor, partner, named pair/tier,
                                  reads, callers, nanomonsv, interpretation,
                                  remove), included in the TSV export and in
                                  Clear all

Nothing is pre-selected: curation is an explicit act, as for variants. The
"Include all reportable" button is the one-click default when that is what
the reader wants.

Usage (from the repo root):
  python3 bin/apply_rx_include.py            apply
  python3 bin/apply_rx_include.py --check    report state, change nothing
  python3 bin/apply_rx_include.py --revert   restore the backups

Backups are written as <file>.rx.bak. The sentinel `rx_include_v1` refuses a
second application.
"""

import argparse
import shutil
import sys
from pathlib import Path

SENTINEL = "rx_include_v1"
BAK_SUFFIX = ".rx.bak"
DASH = "bin/dashboard_builder/"

EDITS = [
    (
        DASH + "parsers/rearrangements.py",
        [
            (
                "        flagged = [m for m in ms if m[\"flag\"]]\n"
                "        return {\n"
                "            \"anchor\": anchor,\n",
                "        flagged = [m for m in ms if m[\"flag\"]]\n"
                "        # rx_include_v1: durable identity for the Reporting tab's store.\n"
                "        # Table row ids are positional and shift when a filter or a\n"
                "        # rebuild reorders loci; gene pair plus both breakend positions\n"
                "        # to 100 kb is stable across rebuilds of the same calls.\n"
                "        key = \"%s::%s|%s:%d|%s:%d\" % (\n"
                "            _mode([m[\"agene\"] for m in ms]) or anchor,\n"
                "            _mode([m[\"pgene\"] for m in ms]) or \"\",\n"
                "            best[\"achrom\"], min(apos) // 100000,\n"
                "            best[\"pchrom\"], min(ppos) // 100000)\n"
                "        return {\n"
                "            \"key\": key,\n"
                "            \"anchor\": anchor,\n",
            ),
        ],
    ),
    (
        DASH + "assets/js/variant-browser.js",
        [
            (
                "    function onChange(cb) { listeners.push(cb); }\n",
                "    // rx_include_v1: interpretation note per included rearrangement locus.\n"
                "    // Same persistence model as CNV captions, separate prefix.\n"
                "    const RX_NOTE_KEY_PREFIX = \"tspipe-rx-note:\";\n"
                "\n"
                "    function getRxNote(id) {\n"
                "      if (!sampleKey || !HAS_STORAGE) return \"\";\n"
                "      try {\n"
                "        return window.localStorage.getItem(RX_NOTE_KEY_PREFIX + sampleKey + \":\" + id) || \"\";\n"
                "      } catch (e) { return \"\"; }\n"
                "    }\n"
                "\n"
                "    function setRxNote(id, text) {\n"
                "      if (!sampleKey || !HAS_STORAGE) return;\n"
                "      try {\n"
                "        if (text) {\n"
                "          window.localStorage.setItem(RX_NOTE_KEY_PREFIX + sampleKey + \":\" + id, text);\n"
                "        } else {\n"
                "          window.localStorage.removeItem(RX_NOTE_KEY_PREFIX + sampleKey + \":\" + id);\n"
                "        }\n"
                "      } catch (e) {}\n"
                "    }\n"
                "\n"
                "    function onChange(cb) { listeners.push(cb); }\n",
            ),
            (
                "      getCnvCaption: getCnvCaption,\n"
                "      setCnvCaption: setCnvCaption,\n",
                "      getCnvCaption: getCnvCaption,\n"
                "      setCnvCaption: setCnvCaption,\n"
                "      getRxNote: getRxNote,\n"
                "      setRxNote: setRxNote,\n",
            ),
        ],
    ),
    (
        DASH + "templates/rearrangements_tab.html.j2",
        [
            # filters bar: include-all button
            (
                "    <span class=\"text-muted ms-auto\" id=\"rx-count\"></span>\n"
                "  </div>\n",
                "    {# rx_include_v1 #}\n"
                "    <button type=\"button\" class=\"btn btn-sm btn-outline-primary ms-2\" id=\"rx-include-reportable\" title=\"tick every reportable locus for the Reporting tab\">Include all reportable</button>\n"
                "    <span class=\"text-muted ms-auto\" id=\"rx-count\"></span>\n"
                "  </div>\n",
            ),
            # locus row: include cell after the flag cell
            (
                "      <td>{% if l.flag %}<span class=\"badge border text-warning border-warning bg-transparent\" title=\"{{ l.flag }}\" style=\"cursor:help\">flag</span>{% endif %}</td>\n"
                "    </tr>\n"
                "    <tr class=\"rx-detail\" data-rx-detail=\"{{ lid }}\" style=\"display:none\">\n"
                "      <td></td>\n"
                "      <td colspan=\"9\" class=\"bg-white small\">\n",
                "      <td>{% if l.flag %}<span class=\"badge border text-warning border-warning bg-transparent\" title=\"{{ l.flag }}\" style=\"cursor:help\">flag</span>{% endif %}</td>\n"
                "      <td class=\"text-center\">\n"
                "        <input type=\"checkbox\" class=\"form-check-input rx-include\" title=\"include in the Reporting tab\"\n"
                "               data-rx-key=\"{{ l.key }}\" data-reportable=\"{{ 1 if l.reportable else 0 }}\"\n"
                "               data-anchor=\"{{ l.anchor_gene }} {{ l.anchor_band }}{% if l.anchor_ig %} {{ l.anchor_ig }}{% endif %}\"\n"
                "               data-partner=\"{{ l.partner_gene }} {{ l.partner_band }}{% if l.partner_ig %} {{ l.partner_ig }}{% endif %}\"\n"
                "               data-pair=\"{{ l.known_pair }}\" data-tier=\"{{ l.tier }}\"\n"
                "               data-reads=\"{{ l.reads if l.reads is not none else '' }}\" data-callers=\"{{ l.callers | join(', ') }}\"\n"
                "               data-nano=\"{{ 'confirmed' if l.nano == 'yes' else (l.nano_filter or ('not reported' if l.nano == 'no' else '')) }}\"\n"
                "               data-point-a=\"{{ l.igv.point_a }}\" data-point-b=\"{{ l.igv.point_b }}\">\n"
                "      </td>\n"
                "    </tr>\n"
                "    <tr class=\"rx-detail\" data-rx-detail=\"{{ lid }}\" style=\"display:none\">\n"
                "      <td></td>\n"
                "      <td colspan=\"10\" class=\"bg-white small\">\n",
            ),
            # member row: one more empty cell
            (
                "          <td>{% if m.flag %}<span class=\"badge border text-warning border-warning bg-transparent\" title=\"{{ m.flag }}\" style=\"cursor:help\">flag</span>{% endif %}</td>\n"
                "        </tr>\n"
                "      {% endfor %}\n",
                "          <td>{% if m.flag %}<span class=\"badge border text-warning border-warning bg-transparent\" title=\"{{ m.flag }}\" style=\"cursor:help\">flag</span>{% endif %}</td>\n"
                "          <td></td>\n"
                "        </tr>\n"
                "      {% endfor %}\n",
            ),
            # header
            (
                "            <th class=\"text-end\">Reads</th><th>Callers</th><th>Worst filter</th><th>nanomonsv</th><th>Flag</th>\n",
                "            <th class=\"text-end\">Reads</th><th>Callers</th><th>Worst filter</th><th>nanomonsv</th><th>Flag</th>\n"
                "            <th title=\"tick to include this locus in the Reporting tab\">Include</th>\n",
            ),
            # group and empty rows span one more column
            (
                "                <td colspan=\"10\" class=\"small\">\n",
                "                <td colspan=\"11\" class=\"small\">\n",
            ),
            (
                "            <tr><td colspan=\"10\" class=\"text-muted small\">none</td></tr>\n",
                "            <tr><td colspan=\"11\" class=\"text-muted small\">none</td></tr>\n",
            ),
            # row click must not toggle details when the checkbox is clicked
            (
                "      if (e.target.closest && (e.target.closest('button') || e.target.closest('a'))) { return; }\n",
                "      if (e.target.closest && (e.target.closest('button') || e.target.closest('a') || e.target.closest('input') || e.target.closest('label'))) { return; }\n",
            ),
            # checkbox wiring
            (
                "    applyFilters();\n"
                "\n"
                "    // Paired breakpoint viewer\n",
                "    applyFilters();\n"
                "\n"
                "    // rx_include_v1: Include checkboxes, backed by the Reporting store.\n"
                "    (function () {\n"
                "      var store = window.tspipeReporting;\n"
                "      if (!store) { return; }\n"
                "      function snapshotOf(box) {\n"
                "        return {\n"
                "          anchor: box.getAttribute('data-anchor') || '', partner: box.getAttribute('data-partner') || '',\n"
                "          pair: box.getAttribute('data-pair') || '', tier: box.getAttribute('data-tier') || '',\n"
                "          reads: box.getAttribute('data-reads') || '', callers: box.getAttribute('data-callers') || '',\n"
                "          nano: box.getAttribute('data-nano') || '',\n"
                "          point_a: box.getAttribute('data-point-a') || '', point_b: box.getAttribute('data-point-b') || '',\n"
                "          reportable: box.getAttribute('data-reportable') === '1'\n"
                "        };\n"
                "      }\n"
                "      function sync() {\n"
                "        var ids = {};\n"
                "        store.load().forEach(function (it) { if (it.kind === 'rearrangement') { ids[it.id] = true; } });\n"
                "        document.querySelectorAll('.rx-include').forEach(function (box) {\n"
                "          var on = !!ids[box.getAttribute('data-rx-key')];\n"
                "          if (box.checked !== on) { box.checked = on; }\n"
                "        });\n"
                "      }\n"
                "      document.addEventListener('change', function (e) {\n"
                "        var box = e.target.closest ? e.target.closest('.rx-include') : null;\n"
                "        if (!box) { return; }\n"
                "        store.toggle('rearrangement', box.getAttribute('data-rx-key'), snapshotOf(box));\n"
                "      });\n"
                "      var all = document.getElementById('rx-include-reportable');\n"
                "      if (all) {\n"
                "        all.addEventListener('click', function () {\n"
                "          document.querySelectorAll('.rx-include[data-reportable=\"1\"]').forEach(function (box) {\n"
                "            if (!box.checked) { store.toggle('rearrangement', box.getAttribute('data-rx-key'), snapshotOf(box)); }\n"
                "          });\n"
                "        });\n"
                "      }\n"
                "      store.onChange(sync);\n"
                "      sync();\n"
                "    })();\n"
                "\n"
                "    // Paired breakpoint viewer\n",
            ),
        ],
    ),
    (
        DASH + "templates/sample_report.html.j2",
        [
            # empty state text
            (
                "          Nothing triaged yet. Open the \"Variants &mdash; Clinical\" tab and mark a\n"
                "          variant <strong>Incl</strong> (include in report) or <strong>Excl</strong>\n"
                "          (reviewed, not reported), and/or open the \"CNV\" tab and tick a plot's\n"
                "          \"Include\" checkbox.\n",
                "          Nothing triaged yet. Tick a locus <strong>Include</strong> box on the\n"
                "          Rearrangements tab, mark a variant <strong>Incl</strong> or\n"
                "          <strong>Excl</strong> on the \"Variants &mdash; Clinical\" tab, or tick a\n"
                "          copy-number plot's \"Include\" checkbox.\n",
            ),
            # rearrangements section markup, before the variants section
            (
                "        {# ----- Variants section ----- #}\n"
                "        <div id=\"reporting-variants-section\" style=\"display:none\">\n",
                "        {# ----- Rearrangements section (rx_include_v1) ----- #}\n"
                "        <div id=\"reporting-rx-section\" style=\"display:none\">\n"
                "          <h5 class=\"mt-2 mb-2\">\n"
                "            Rearrangements <span id=\"reporting-rx-count-badge\" class=\"badge bg-secondary ms-1\">0</span>\n"
                "          </h5>\n"
                "          <p class=\"text-muted small mb-2\">\n"
                "            Loci ticked on the Rearrangements tab, in selection order. The\n"
                "            interpretation is editable and persists in this browser; it is\n"
                "            included in the copied TSV.\n"
                "          </p>\n"
                "          <div class=\"table-responsive\">\n"
                "            <table class=\"table table-sm table-hover align-middle\" id=\"reporting-rx-table\">\n"
                "              <thead class=\"table-light\">\n"
                "                <tr>\n"
                "                  <th scope=\"col\" style=\"width:32px\"></th>\n"
                "                  <th scope=\"col\">Anchor</th>\n"
                "                  <th scope=\"col\">Partner</th>\n"
                "                  <th scope=\"col\">Named pair / tier</th>\n"
                "                  <th scope=\"col\">Breakends</th>\n"
                "                  <th scope=\"col\">Reads</th>\n"
                "                  <th scope=\"col\">Callers</th>\n"
                "                  <th scope=\"col\">nanomonsv</th>\n"
                "                  <th scope=\"col\" style=\"min-width:220px\">Interpretation</th>\n"
                "                  <th scope=\"col\" style=\"width:32px\"></th>\n"
                "                </tr>\n"
                "              </thead>\n"
                "              <tbody id=\"reporting-rx-tbody\"></tbody>\n"
                "            </table>\n"
                "          </div>\n"
                "        </div>\n"
                "\n"
                "        {# ----- Variants section ----- #}\n"
                "        <div id=\"reporting-variants-section\" style=\"display:none\">\n",
            ),
            # element handles
            (
                "    const toast = document.getElementById(\"reporting-copied-toast\");\n"
                "    if (!tbody || !window.tspipeReporting) return;\n",
                "    const toast = document.getElementById(\"reporting-copied-toast\");\n"
                "    const rxSection = document.getElementById(\"reporting-rx-section\");      // rx_include_v1\n"
                "    const rxTbody = document.getElementById(\"reporting-rx-tbody\");\n"
                "    const rxCountBadge = document.getElementById(\"reporting-rx-count-badge\");\n"
                "    if (!tbody || !window.tspipeReporting) return;\n",
            ),
            # columns + row-to-tsv + renderRx, before renderVariants
            (
                "    function renderVariants(variants) {\n",
                "    const RX_COLUMNS = [\"Anchor\", \"Partner\", \"Named pair\", \"Tier\", \"Breakend A\", \"Breakend B\",\n"
                "                        \"Reads\", \"Callers\", \"nanomonsv\", \"Interpretation\"];\n"
                "\n"
                "    function rxRowToTsv(it) {\n"
                "      const s = it.snapshot || {};\n"
                "      const note = (window.tspipeReporting.getRxNote(it.id) || \"\").replace(/\\s+/g, \" \");\n"
                "      return [s.anchor || \"\", s.partner || \"\", s.pair || \"\", s.tier || \"\", s.point_a || \"\", s.point_b || \"\",\n"
                "              s.reads || \"\", s.callers || \"\", s.nano || \"\", note].join(\"\\t\");\n"
                "    }\n"
                "\n"
                "    function renderRx(rxs) {\n"
                "      if (!rxSection) return;\n"
                "      if (rxs.length === 0) {\n"
                "        rxSection.style.display = \"none\";\n"
                "        rxTbody.innerHTML = \"\";\n"
                "        return;\n"
                "      }\n"
                "      rxSection.style.display = \"\";\n"
                "      rxCountBadge.textContent = String(rxs.length);\n"
                "      rxTbody.innerHTML = rxs.map(function (it, idx) {\n"
                "        const s = it.snapshot || {};\n"
                "        const note = escapeHtml(window.tspipeReporting.getRxNote(it.id) || \"\");\n"
                "        const pair = escapeHtml(s.pair || \"\") + (s.tier ? ' <span class=\"badge bg-secondary\">' + escapeHtml(s.tier) + '</span>' : \"\");\n"
                "        return '<tr data-rx-key=\"' + escapeHtml(it.id) + '\">' +\n"
                "                 '<td class=\"text-muted small\">' + (idx + 1) + '</td>' +\n"
                "                 '<td>' + escapeHtml(s.anchor) + '</td>' +\n"
                "                 '<td>' + escapeHtml(s.partner) + '</td>' +\n"
                "                 '<td>' + pair + '</td>' +\n"
                "                 '<td class=\"font-monospace small\">' + escapeHtml(s.point_a) + '<br>' + escapeHtml(s.point_b) + '</td>' +\n"
                "                 '<td>' + escapeHtml(s.reads) + '</td>' +\n"
                "                 '<td class=\"small\">' + escapeHtml(s.callers) + '</td>' +\n"
                "                 '<td class=\"small\">' + escapeHtml(s.nano) + '</td>' +\n"
                "                 '<td><input type=\"text\" class=\"form-control form-control-sm reporting-rx-note-input\" ' +\n"
                "                       'data-rx-key=\"' + escapeHtml(it.id) + '\" placeholder=\"interpretation\" value=\"' + note + '\" maxlength=\"300\"></td>' +\n"
                "                 '<td><button type=\"button\" class=\"btn btn-sm btn-outline-secondary reporting-rx-remove\" ' +\n"
                "                       'data-rx-key=\"' + escapeHtml(it.id) + '\" title=\"Remove from report\">&times;</button></td>' +\n"
                "               '</tr>';\n"
                "      }).join(\"\");\n"
                "    }\n"
                "\n"
                "    function renderVariants(variants) {\n",
            ),
            # render(): count and sections
            (
                "      const cnvs = all.filter(function (it) { return it.kind === \"cnv_plot\"; });\n"
                "      const total = variants.length + excluded.length + cnvs.length;\n",
                "      const cnvs = all.filter(function (it) { return it.kind === \"cnv_plot\"; });\n"
                "      const rxs = all.filter(function (it) { return it.kind === \"rearrangement\"; });   // rx_include_v1\n"
                "      const total = variants.length + excluded.length + cnvs.length + rxs.length;\n",
            ),
            (
                "        variantsSection.style.display = \"none\";\n"
                "        excludedSection.style.display = \"none\";\n"
                "        cnvSection.style.display = \"none\";\n"
                "        copyBtn.disabled = true;\n",
                "        variantsSection.style.display = \"none\";\n"
                "        excludedSection.style.display = \"none\";\n"
                "        cnvSection.style.display = \"none\";\n"
                "        if (rxSection) rxSection.style.display = \"none\";\n"
                "        copyBtn.disabled = true;\n",
            ),
            (
                "      renderVariants(variants);\n"
                "      renderExcluded(excluded);\n"
                "      renderCnvs(cnvs);\n",
                "      renderRx(rxs);\n"
                "      renderVariants(variants);\n"
                "      renderExcluded(excluded);\n"
                "      renderCnvs(cnvs);\n",
            ),
            # buildTsv
            (
                "      const cnvs = all.filter(function (it) { return it.kind === \"cnv_plot\"; });\n"
                "      const sections = [];\n",
                "      const cnvs = all.filter(function (it) { return it.kind === \"cnv_plot\"; });\n"
                "      const rxs = all.filter(function (it) { return it.kind === \"rearrangement\"; });\n"
                "      const sections = [];\n"
                "      if (rxs.length) {\n"
                "        sections.push(\"# Rearrangements\");\n"
                "        sections.push(RX_COLUMNS.join(\"\\t\"));\n"
                "        rxs.forEach(function (it) { sections.push(rxRowToTsv(it)); });\n"
                "      }\n",
            ),
            (
                "      if (variants.length) {\n"
                "        sections.push(\"# Variants\");\n",
                "      if (variants.length) {\n"
                "        if (sections.length) sections.push(\"\");\n"
                "        sections.push(\"# Variants\");\n",
            ),
            # clear all wipes notes
            (
                "        window.tspipeReporting.clearAll();\n",
                "        window.tspipeReporting.load()\n"
                "          .filter(function (it) { return it.kind === \"rearrangement\"; })\n"
                "          .forEach(function (it) { window.tspipeReporting.setRxNote(it.id, \"\"); });\n"
                "        window.tspipeReporting.clearAll();\n",
            ),
            # note input and remove handlers
            (
                "    // Re-render on selection change\n"
                "    window.tspipeReporting.onChange(render);\n",
                "    // rx_include_v1: interpretation edits persist on input; remove toggles off.\n"
                "    if (rxTbody) {\n"
                "      rxTbody.addEventListener(\"input\", function (ev) {\n"
                "        const inp = ev.target.closest(\".reporting-rx-note-input\");\n"
                "        if (!inp) return;\n"
                "        window.tspipeReporting.setRxNote(inp.getAttribute(\"data-rx-key\"), inp.value);\n"
                "      });\n"
                "      rxTbody.addEventListener(\"click\", function (ev) {\n"
                "        const btn = ev.target.closest(\".reporting-rx-remove\");\n"
                "        if (!btn) return;\n"
                "        const key = btn.getAttribute(\"data-rx-key\");\n"
                "        const item = window.tspipeReporting.load().find(function (it) {\n"
                "          return it.kind === \"rearrangement\" && it.id === key;\n"
                "        });\n"
                "        if (item) {\n"
                "          window.tspipeReporting.setRxNote(key, \"\");\n"
                "          window.tspipeReporting.toggle(\"rearrangement\", key, item.snapshot);\n"
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
    ap.add_argument("--root", default=".", help="repository root (default: cwd)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true", help="report state only")
    g.add_argument("--revert", action="store_true", help="restore backup copies")
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
