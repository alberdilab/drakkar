"""Render ``drakkar_report.html`` from the report database.

The renderer reads ``drakkar.db`` and nothing else: every number on the page
comes from a SQL query, so the source ``.tsv``/``.tsv.xz`` tables are read
exactly once, at ingest, and never again here.

Two constraints shape the module. First, a Drakkar output directory rarely
holds every module's output, so a section is rendered only when the database
actually contains it, and whatever is absent is named on the page rather than
raised. Second, ``gene_annotation``, ``cluster_annotation`` and
``gene_expression`` can each hold tens of millions of rows, so aggregation
happens in SQL and only summary-sized results are ever pulled into memory.

The page is self-contained: the stylesheet and the small navigation script are
inlined, and the Plotly bundle is embedded once, in the first figure, so the
report opens offline. Layout-wise it is a left sidebar — report details and the
table of contents — beside a panel area that shows one section at a time; long
tables are paginated in the browser and preceded by the averages of their
numeric columns. Without JavaScript the same markup degrades to the flat page
it used to be: every section stacked, every row listed.
"""

import base64
import math
import re
import sqlite3
from datetime import datetime, timezone
from html import escape
from pathlib import Path

from drakkar.report.schema import SCHEMA_VERSION, TAXONOMIC_RANKS

# Category caps for the figures. A heatmap or a bar chart cannot be paged, so
# what it shows is bounded here rather than in SQL scattered across the section
# renderers. Tables are not capped: they carry every row and are paged in the
# browser instead.
# Rows per page in the browser, and the smallest table worth averaging.
PAGE_ROWS = 20
MIN_STAT_ROWS = 2
TOP_GENOMES = 40
TOP_MAGS = 40
# How many taxa a stacked composition chart names before pooling the rest as
# Other. Set by what its key can hold: the key is one row spread across the
# page, so each entry gets a page-width fraction, and past eight or so a long
# phylum name no longer fits in its share and is cut off mid-word. Staying
# under the palette also keeps every named taxon a colour of its own, which a
# stacked bar needs — two slices in the same colour are read as one.
TOP_TAXA = 8
TOP_CATEGORIES = 30

FIGURE_HEIGHT = 420

# Set as a figure's ``layout.meta`` to have it rendered edge to edge instead of
# inside the content column's padding.
WIDE = "wide"

# Muted palette, deliberately close to the terminal theme in drakkar.output.
PALETTE = (
    "#5f9ea0", "#d6a642", "#7fb069", "#e85d75", "#8878b0", "#4f7f9f",
    "#c98b5e", "#6fae9a", "#b0637f", "#8d9db6",
)

# The neutral a figure gives to everything it pooled, or to a category it is
# not naming. It stays out of PALETTE so no named category is ever drawn in it.
OTHER_COLOUR = "#b7c7d3"

SECTION_LABELS = {
    "preprocessing": "Preprocessing",
    "cataloging": "Cataloging",
    "dereplication": "Dereplication",
    "profiling": "Profiling",
    "taxonomy": "Taxonomy",
    "function": "Functional annotation",
    "expression": "Expression",
    "amr": "Antimicrobial resistance",
    "resources": "Runs and resources",
}

# One plain-language sentence per section, shown under its heading. The
# report is read by people who did not run the workflow, so each section
# says what its numbers are before showing them.
SECTION_INTROS = {
    "preprocessing": (
        "What happened to the raw sequencing reads of each sample before "
        "anything was assembled: adapter and quality trimming with fastp, "
        "then removal of reads matching the host genome. The metagenomic "
        "reads that survive both steps are the input to every later section."
    ),
    "cataloging": (
        "Assembly of the metagenomic reads into contigs, and the recovery of "
        "genome bins from those contigs. An assembly may be built from one "
        "sample or from several pooled together, which is why samples and "
        "assemblies are reported in separate tables."
    ),
    "dereplication": (
        "Bins recovered from different assemblies often represent the same "
        "organism. Dereplication collapses them into one representative "
        "genome per cluster, so the drop in bin count here is redundancy "
        "being removed, not genomes being lost."
    ),
    "profiling": (
        "The final genome catalogue, and how much of each sample maps to it. "
        "Quality figures come from CheckM2 estimates: completeness is how "
        "much of the expected genome was recovered, contamination how much "
        "of it likely belongs to another organism."
    ),
    "taxonomy": (
        "Which organisms the catalogue genomes correspond to, according to "
        "GTDB-Tk. A genome left unclassified at a rank has no close enough "
        "relative in the reference database, which is common for "
        "under-sampled environments and is not an error."
    ),
    "function": (
        "What the genes in the catalogue are predicted to do. Each "
        "annotation source is a different reference database, and a gene can "
        "be annotated by several of them at once, so figures from different "
        "sources are not additive."
    ),
    "expression": (
        "Metatranscriptomic reads mapped onto the catalogue genes, so these "
        "counts say how much each gene was transcribed, not whether it is "
        "present. A gene can be present in the catalogue and still have zero "
        "counts in a sample."
    ),
    "amr": (
        "Antimicrobial resistance genes called on assembled contigs by "
        "AMRFinderPlus and RGI, reconciled into one locus per gene occurrence "
        "and placed against the plasmids and viruses geNomad found. This "
        "section runs on assemblies rather than on the catalogue genomes, so "
        "its counts do not line up with the genome tables above."
    ),
    "resources": (
        "How the workflow itself ran: which commands were executed and what "
        "compute they asked for against what they actually used. Nothing "
        "here affects the biological results; it is for tuning the resource "
        "profile and diagnosing failed jobs."
    ),
}

STYLESHEET = """
:root {
  --ink: #22303c;
  --muted: #6b7a88;
  --rule: #d8dee4;
  --accent: #5f9ea0;
  --panel: #f6f8f9;
  --sidebar: 20.5rem;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
  font-size: 15px;
  line-height: 1.55;
  color: var(--ink);
  background: #ffffff;
}
.layout { display: flex; align-items: flex-start; gap: 0; }
main.content {
  flex: 1 1 auto;
  min-width: 0;
  max-width: 1180px;
  padding: 2.25rem 2rem 4rem;
}
h2 { margin: 0 0 .75rem; font-size: 1.35rem; padding-bottom: .3rem; border-bottom: 1px solid var(--rule); }
h3 { margin: 1.75rem 0 .4rem; font-size: 1.05rem; color: var(--ink); }
p { margin: .5rem 0; }
p.note { color: var(--muted); font-size: .88rem; }

/* Sidebar: report identity, navigation, and everything about the run itself. */
aside.sidebar {
  flex: 0 0 var(--sidebar);
  width: var(--sidebar);
  position: sticky;
  top: 0;
  height: 100vh;
  overflow-y: auto;
  padding: 1.75rem 1.4rem 2.5rem;
  background: var(--panel);
  border-right: 1px solid var(--rule);
}
aside.sidebar .brand {
  border-bottom: 2px solid var(--accent);
  padding-bottom: .9rem;
  text-align: center;
}
aside.sidebar .logo {
  display: block;
  width: 100%;
  max-width: 14rem;
  height: auto;
  margin: 0 auto .9rem;
  /* The mark ships on white; multiply drops that onto the sidebar's panel. */
  mix-blend-mode: multiply;
}
aside.sidebar h1 { margin: 0 0 .2rem; font-size: 1.35rem; letter-spacing: .01em; }
aside.sidebar .subtitle { margin: 0; color: var(--muted); font-size: .82rem; }
aside.sidebar .side-block { margin-top: 1.6rem; }
aside.sidebar h2 {
  margin: 0 0 .5rem;
  font-size: .72rem;
  text-transform: uppercase;
  letter-spacing: .09em;
  color: var(--muted);
  border: 0;
  padding: 0;
}
nav.toc ol { margin: 0; padding: 0; list-style: none; counter-reset: toc; }
nav.toc li { counter-increment: toc; }
nav.toc a {
  display: block;
  padding: .3rem .6rem .3rem .55rem;
  color: var(--ink);
  text-decoration: none;
  border-left: 3px solid transparent;
  font-size: .92rem;
}
nav.toc a::before { content: counter(toc) ". "; color: var(--muted); font-size: .82rem; }
nav.toc a:hover { color: var(--accent); background: #ffffff; }
nav.toc a.is-active {
  border-left-color: var(--accent);
  background: #ffffff;
  font-weight: 600;
}
dl.summary { display: grid; grid-template-columns: 1fr; gap: 0; margin: 0; }
dl.summary dt {
  color: var(--muted);
  font-size: .68rem;
  text-transform: uppercase;
  letter-spacing: .06em;
  margin-top: .6rem;
}
dl.summary dd { margin: 0; font-size: .86rem; overflow-wrap: anywhere; }
.skipped { margin-top: 1.6rem; padding-top: 1rem; border-top: 1px solid var(--rule); }
.skipped p { font-size: .82rem; color: var(--muted); margin: .4rem 0 0; }
.skipped ul { margin: .3rem 0 0; padding-left: 1.1rem; color: var(--muted); font-size: .82rem; }

/* One section at a time, once the script has taken over. */
html.js main .panel { display: none; }
html.js main .panel.is-active { display: block; }
html:not(.js) main .panel + .panel { margin-top: 3rem; }

/* Averages of the numeric columns, above the table they summarize. */
.stats { display: flex; flex-wrap: wrap; gap: .6rem; margin: 1rem 0 .4rem; }
.stat {
  flex: 1 1 10.5rem;
  max-width: 18rem;
  background: var(--panel);
  border: 1px solid var(--rule);
  border-left: 3px solid var(--accent);
  padding: .5rem .75rem;
}
.stat-label {
  display: block;
  font-size: .68rem;
  text-transform: uppercase;
  letter-spacing: .06em;
  color: var(--muted);
}
.stat-value { display: block; font-size: 1.2rem; font-variant-numeric: tabular-nums; }
.stat-sub {
  display: block;
  font-size: .78rem;
  color: var(--muted);
  font-variant-numeric: tabular-nums;
}

table { border-collapse: collapse; width: 100%; margin: .75rem 0; font-size: .9rem; }
.table-block { margin: .5rem 0 1rem; }
.scroll { overflow-x: auto; }
th, td { padding: .35rem .6rem; text-align: right; border-bottom: 1px solid var(--rule); white-space: nowrap; }
th { background: var(--panel); font-weight: 600; color: var(--muted); text-transform: uppercase; font-size: .72rem; letter-spacing: .05em; }
/* Sorting is offered only once the script has made the headers interactive. */
th.sortable { cursor: pointer; user-select: none; }
th.sortable:hover { color: var(--accent); }
th.sortable::after { content: "\2195"; margin-left: .3em; opacity: .3; font-size: .9em; }
th.sortable[aria-sort="ascending"]::after { content: "\2191"; opacity: 1; color: var(--accent); }
th.sortable[aria-sort="descending"]::after { content: "\2193"; opacity: 1; color: var(--accent); }
th.sortable:focus-visible { outline: 2px solid var(--accent); outline-offset: -2px; }
th:first-child, td:first-child { text-align: left; }
tbody tr:hover td { background: #fbfcfd; }

.pager { display: flex; align-items: center; gap: .6rem; margin: -.25rem 0 .5rem; font-size: .82rem; color: var(--muted); }
.pager button {
  font: inherit;
  color: var(--ink);
  background: #ffffff;
  border: 1px solid var(--rule);
  padding: .18rem .6rem;
  cursor: pointer;
}
.pager button:hover:not(:disabled) { border-color: var(--accent); color: var(--accent); }
.pager button:disabled { opacity: .4; cursor: default; }

/* A subsection title with a control of its own against the right edge. */
.head-row { display: flex; align-items: baseline; justify-content: space-between; gap: 1rem; }
.head-row h3 { margin-bottom: 0; }
.unit-toggle {
  flex: 0 0 auto;
  font: inherit;
  font-size: .82rem;
  color: var(--ink);
  background: #ffffff;
  border: 1px solid var(--rule);
  padding: .18rem .7rem;
  cursor: pointer;
}
.unit-toggle:hover { border-color: var(--accent); color: var(--accent); }
.unit-toggle:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
/* Nothing to switch until the script is running, and nothing to click on paper. */
html:not(.js) .unit-toggle { display: none; }

.figure { margin: 1rem 0 .25rem; }

/* The circular phylogeny. The SVG scales to its column but is capped: past
   about 46rem the rings grow without any more of the tree becoming legible,
   and the tip labels drift away from the arcs they belong to. */
.tree-figure { margin: 1rem 0 1.25rem; }
.tree-figure svg.tree {
  display: block;
  width: 100%;
  max-width: 46rem;
  height: auto;
  margin: 0 auto;
  font-family: inherit;
}
.tree-key { max-width: 46rem; margin: .25rem auto 0; }
.tree-key .key-row {
  display: flex;
  flex-wrap: wrap;
  gap: .25rem .9rem;
  margin-top: .35rem;
  font-size: .8rem;
  color: var(--muted);
}
.tree-key .key-item { display: inline-flex; align-items: center; gap: .3rem; }
.tree-key .key-swatch {
  width: .72rem;
  height: .72rem;
  border: 1px solid rgba(34, 48, 60, .12);
}
.tree-key .key-ramp {
  width: 2.4rem;
  height: .72rem;
  border: 1px solid rgba(34, 48, 60, .12);
}
.tree-key p.note { margin: .45rem 0 0; font-size: .78rem; }
/* A figure that carries a whole section reads better edge to edge, so it is
   let out of the cap the text column keeps for line length: it starts where
   the text does and runs to the far side of the page, leaving the same
   gutter there. Below 900px the sidebar is no longer beside it and the
   figure goes back to the width of its column. */
.figure.wide, .head-row.wide {
  width: calc(100vw - var(--sidebar) - 4rem);
  max-width: none;
}
footer.report { margin-top: 3.5rem; padding-top: 1rem; border-top: 1px solid var(--rule); color: var(--muted); font-size: .85rem; }
code { font-family: "SFMono-Regular", Consolas, monospace; font-size: .85em; }

@media (max-width: 900px) {
  .layout { display: block; }
  nav.toc ol { display: flex; flex-wrap: wrap; gap: .2rem .4rem; }
  nav.toc a { border-left: 0; border-bottom: 2px solid transparent; padding: .25rem .4rem; }
  nav.toc a.is-active { border-bottom-color: var(--accent); }
  aside.sidebar {
    position: static;
    width: auto;
    height: auto;
    border-right: 0;
    border-bottom: 1px solid var(--rule);
  }
  main.content { padding: 1.5rem 1.25rem 3rem; }
  .figure.wide, .head-row.wide { width: auto; }
}

/* Printing wants the whole report, not the panel that happens to be open. */
@media print {
  aside.sidebar { display: none; }
  .figure.wide, .head-row.wide { width: auto; }
  .unit-toggle { display: none; }
  html.js main .panel { display: block !important; }
  .pager { display: none; }
  main .panel + .panel { margin-top: 3rem; page-break-before: always; }
}
"""

# The page works without this script — every section stacked, every row shown,
# in the order the renderer chose. It adds the three behaviours the markup only
# hints at: the sidebar switches panels instead of scrolling to them, long
# tables are paged, and any column can be sorted on.
SCRIPT = """
(function () {
  var panels = Array.prototype.slice.call(document.querySelectorAll("main .panel"));
  var links = Array.prototype.slice.call(document.querySelectorAll("nav.toc a"));

  function resizePlots(scope) {
    if (!window.Plotly || !window.Plotly.Plots) { return; }
    var plots = scope.querySelectorAll(".js-plotly-plot");
    for (var index = 0; index < plots.length; index += 1) {
      try { window.Plotly.Plots.resize(plots[index]); } catch (error) { /* ignore */ }
    }
  }

  function show(id, scroll) {
    var target = null;
    for (var index = 0; index < panels.length; index += 1) {
      if (panels[index].id === id) { target = panels[index]; }
    }
    if (!target) { return false; }
    for (var p = 0; p < panels.length; p += 1) {
      panels[p].classList.toggle("is-active", panels[p] === target);
    }
    for (var l = 0; l < links.length; l += 1) {
      links[l].classList.toggle(
        "is-active", links[l].getAttribute("href") === "#" + id
      );
    }
    // Figures drawn inside a hidden panel have no width to measure until now.
    resizePlots(target);
    if (scroll) { window.scrollTo(0, 0); }
    return true;
  }

  for (var index = 0; index < links.length; index += 1) {
    links[index].addEventListener("click", function (event) {
      var id = this.getAttribute("href").slice(1);
      if (!show(id, true)) { return; }
      event.preventDefault();
      try { window.history.replaceState(null, "", "#" + id); } catch (error) { /* file:// */ }
    });
  }

  window.addEventListener("hashchange", function () {
    show(window.location.hash.slice(1), true);
  });

  function count(value) {
    return value.toLocaleString ? value.toLocaleString("en-US") : String(value);
  }

  // A cell sorts on its raw value when it has one, and on its text otherwise.
  // Blanks are missing values, not zeros or empty strings, so they sit at the
  // bottom whichever way the column is sorted.
  function key(row, index) {
    var cell = row.cells[index];
    if (!cell) { return null; }
    var raw = cell.getAttribute("data-sort");
    if (raw !== null && raw !== "") {
      var number = parseFloat(raw);
      return isNaN(number) ? null : number;
    }
    var text = (cell.textContent || "").trim();
    return text === "" ? null : text;
  }

  function compare(left, right) {
    if (typeof left === "number" && typeof right === "number") {
      return left - right;
    }
    return String(left).localeCompare(String(right), undefined, {
      numeric: true, sensitivity: "base"
    });
  }

  // Sorting reorders the shared row list, so paging keeps working on it.
  function sortable(table, rows, redraw) {
    var head = table.tHead;
    if (!head || !head.rows.length || rows.length < 2) { return; }
    var headers = Array.prototype.slice.call(head.rows[0].cells);

    function sort(index, header) {
      var descending = header.getAttribute("aria-sort") === "ascending";
      var decorated = rows.map(function (row, position) {
        return { row: row, position: position, value: key(row, index) };
      });
      decorated.sort(function (a, b) {
        if (a.value === null && b.value === null) { return a.position - b.position; }
        if (a.value === null) { return 1; }
        if (b.value === null) { return -1; }
        var order = compare(a.value, b.value);
        if (order === 0) { return a.position - b.position; }
        return descending ? -order : order;
      });
      var body = table.tBodies[0];
      for (var index2 = 0; index2 < decorated.length; index2 += 1) {
        rows[index2] = decorated[index2].row;
        body.appendChild(decorated[index2].row);
      }
      for (var h = 0; h < headers.length; h += 1) {
        headers[h].setAttribute("aria-sort", "none");
      }
      header.setAttribute("aria-sort", descending ? "descending" : "ascending");
      redraw();
    }

    headers.forEach(function (header, index) {
      header.classList.add("sortable");
      header.setAttribute("role", "button");
      header.setAttribute("tabindex", "0");
      header.setAttribute("aria-sort", "none");
      header.setAttribute("title", "Sort by " + (header.textContent || "").trim());
      header.addEventListener("click", function () { sort(index, header); });
      header.addEventListener("keydown", function (event) {
        if (event.key === "Enter" || event.key === " " || event.key === "Spacebar") {
          event.preventDefault();
          sort(index, header);
        }
      });
    });
  }

  // Client-side paging: the rows are all in the document, only 20 are visible.
  function paginate(table, rows) {
    var size = parseInt(table.getAttribute("data-page-size"), 10);
    if (!size || rows.length <= size) { return null; }
    var pages = Math.ceil(rows.length / size);
    var current = 1;

    var pager = document.createElement("div");
    pager.className = "pager";
    var previous = document.createElement("button");
    previous.type = "button";
    previous.textContent = "\u2039 Previous";
    var next = document.createElement("button");
    next.type = "button";
    next.textContent = "Next \u203a";
    var label = document.createElement("span");
    pager.appendChild(previous);
    pager.appendChild(next);
    pager.appendChild(label);

    function draw() {
      var first = (current - 1) * size;
      var last = Math.min(first + size, rows.length);
      for (var index = 0; index < rows.length; index += 1) {
        rows[index].style.display = (index >= first && index < last) ? "" : "none";
      }
      label.textContent = "Rows " + count(first + 1) + "\u2013" + count(last) +
        " of " + count(rows.length) + " \u00b7 page " + current + " of " + pages;
      previous.disabled = current === 1;
      next.disabled = current === pages;
    }

    previous.addEventListener("click", function () {
      if (current > 1) { current -= 1; draw(); }
    });
    next.addEventListener("click", function () {
      if (current < pages) { current += 1; draw(); }
    });

    var scroll = table.parentNode;
    scroll.parentNode.insertBefore(pager, scroll.nextSibling);
    draw();
    // Re-sorting changes which rows belong on the open page, so it returns to
    // the first one rather than leaving the reader mid-table.
    return function () { current = 1; draw(); };
  }

  function enhance(table) {
    var body = table.tBodies[0];
    if (!body) { return; }
    var rows = Array.prototype.slice.call(body.rows);
    var redraw = paginate(table, rows) || function () {};
    sortable(table, rows, redraw);
  }

  var tables = document.querySelectorAll("main table");
  for (var t = 0; t < tables.length; t += 1) { enhance(tables[t]); }

  // A figure whose traces carry the same quantity in a second unit, in
  // trace.meta.alt: clicking the button swaps the two around, so the plot is
  // redrawn in place instead of a second copy of it being shipped.
  function unitToggle(button) {
    var target = document.getElementById(button.getAttribute("data-target"));
    var plot = target ? target.querySelector(".js-plotly-plot") : null;
    if (!plot || !window.Plotly) { button.style.display = "none"; return; }

    function swappable() {
      var data = plot.data || [];
      var found = [];
      for (var index = 0; index < data.length; index += 1) {
        var meta = data[index].meta;
        if (meta && meta.alt) { found.push(index); }
      }
      return found;
    }

    if (!swappable().length) { button.style.display = "none"; return; }
    var showingAlt = false;

    button.addEventListener("click", function () {
      var data = plot.data || [];
      var indices = swappable();
      var swapped = [];
      for (var i = 0; i < indices.length; i += 1) {
        var trace = data[indices[i]];
        swapped.push(trace.meta.alt);
        trace.meta.alt = trace.y;
      }
      showingAlt = !showingAlt;
      window.Plotly.restyle(plot, { y: swapped }, indices);
      window.Plotly.relayout(plot, {
        "yaxis.title.text": button.getAttribute(
          showingAlt ? "data-alt-axis" : "data-axis"
        )
      });
      button.textContent = "Switch to " + button.getAttribute(
        showingAlt ? "data-unit" : "data-alt-unit"
      );
    });
  }

  var toggles = document.querySelectorAll("button.unit-toggle");
  for (var u = 0; u < toggles.length; u += 1) { unitToggle(toggles[u]); }

  if (!show(window.location.hash.slice(1), false) && panels.length) {
    show(panels[0].id, false);
  }
})();
"""


# ---------------------------------------------------------------------------
# Value and markup helpers
# ---------------------------------------------------------------------------

def _text(value):
    return "" if value is None else escape(str(value))


LOGO_PATH = Path(__file__).parent / "assets" / "drakkar.png"


def _logo_html():
    """The Drakkar mark, base64-inlined so the report stays one portable file."""
    try:
        encoded = base64.b64encode(LOGO_PATH.read_bytes()).decode("ascii")
    except OSError:
        return ""
    return f'<img class="logo" src="data:image/png;base64,{encoded}" alt="Drakkar">'


def _parse_timestamp(value):
    """Read an ISO timestamp as it is stored in the database, in UTC."""
    if value is None:
        return None
    try:
        moment = datetime.fromisoformat(str(value).strip().replace("Z", "+00:00"))
    except ValueError:
        return None
    if moment.tzinfo is None:
        return moment.replace(tzinfo=timezone.utc)
    return moment.astimezone(timezone.utc)


def _friendly_datetime(value):
    """``2026-08-25T13:18:19.205575+00:00`` reads as ``25 Aug 2026 at 13:18 UTC``."""
    moment = _parse_timestamp(value)
    if moment is None:
        return "" if value is None else str(value)
    return f"{moment.day} {moment:%b %Y} at {moment:%H:%M} UTC"


def _timestamp(value):
    return escape(_friendly_datetime(value))


def _integer(value):
    if value is None:
        return ""
    try:
        return f"{int(round(float(value))):,}"
    except (TypeError, ValueError):
        return escape(str(value))


def _decimal(digits):
    """Return a formatter that prints a float with a fixed number of digits."""
    def formatter(value):
        if value is None:
            return ""
        try:
            return f"{float(value):,.{digits}f}"
        except (TypeError, ValueError):
            return escape(str(value))
    return formatter


_ONE = _decimal(1)
_TWO = _decimal(2)


def _gigabases(value):
    """Base counts read better in gigabases: 1,850,000,000 becomes 1.85 GB."""
    if value is None or value == "":
        return ""
    try:
        return f"{float(value) / 1e9:,.2f} GB"
    except (TypeError, ValueError):
        return escape(str(value))


def _share(value):
    """A share as a percentage, whichever way the tool wrote it.

    SingleM and Nonpareil both report shares as fractions of one, so that is
    what ingest stores; anything above one is already a percentage and is left
    alone rather than multiplied into nonsense.
    """
    if value is None or value == "":
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return escape(str(value))
    return f"{(number * 100 if number <= 1 else number):,.1f}%"


def _cards_percent(value):
    """A highlight card's percentage, or None when there is nothing to average."""
    return None if value is None else f"{_ONE(value)}%"


def _quantity(number, singular, plural=None):
    """Format a count with a noun that agrees with it."""
    noun = singular if number == 1 else (plural or singular + "s")
    return f"{number:,} {noun}"


def _ratio(numerator, denominator):
    """Percent of numerator over denominator, or None when it is undefined."""
    if numerator is None or not denominator:
        return None
    return 100.0 * float(numerator) / float(denominator)


def _percentile(values, fraction):
    """Linear-interpolated percentile of the numeric values, blanks ignored.

    Written out rather than taken from ``statistics`` because a rule can have
    launched a single job, and the quantile helpers there refuse fewer than
    two data points.
    """
    numbers = []
    for value in values:
        if value is None or value == "":
            continue
        try:
            numbers.append(float(value))
        except (TypeError, ValueError):
            continue
    if not numbers:
        return None
    numbers.sort()
    if len(numbers) == 1:
        return numbers[0]
    position = fraction * (len(numbers) - 1)
    lower = int(position)
    upper = min(lower + 1, len(numbers) - 1)
    weight = position - lower
    return numbers[lower] * (1 - weight) + numbers[upper] * weight


def _mean(rows, index):
    """Average of a column, ignoring blanks and anything that is not a number."""
    values = []
    for row in rows:
        value = row[index]
        if value is None or value == "":
            continue
        try:
            values.append(float(value))
        except (TypeError, ValueError):
            continue
    if not values:
        return None
    return sum(values) / len(values)


def _highlights(items):
    """Render highlight cards from ``(label, value)`` triples.

    An optional third element is a quieter line under the value — the same
    quantity in other units, say — and is already formatted when it arrives.
    """
    cards = []
    for item in items:
        label, value = item[0], item[1]
        if value in (None, ""):
            continue
        sub = item[2] if len(item) > 2 else None
        under = f'<span class="stat-sub">{sub}</span>' if sub else ""
        cards.append(
            f'<div class="stat"><span class="stat-label">{escape(label)}</span>'
            f'<span class="stat-value">{value}</span>{under}</div>'
        )
    if not cards:
        return ""
    return '<div class="stats">' + "".join(cards) + "</div>"


def _stats_html(columns, rows, stats):
    """Averages of the named columns, computed over every row, not just shown ones.

    ``stats`` is a sequence of ``(label, column index)`` pairs, optionally with
    a third element overriding the column's own formatter.
    """
    if not stats or len(rows) < MIN_STAT_ROWS:
        return ""
    items = []
    for spec in stats:
        label, index = spec[0], spec[1]
        formatter = spec[2] if len(spec) > 2 else columns[index][1]
        value = _mean(rows, index)
        if value is not None:
            items.append((label, formatter(value)))
    return _highlights(items)


def _table(columns, rows, stats=()):
    """Render rows as an HTML table, dropping columns that are entirely empty.

    ``columns`` is a sequence of ``(header, formatter)`` pairs and ``rows`` a
    sequence of raw value tuples in the same order. ``stats`` names the columns
    whose average belongs above the table. Every row is rendered: tables longer
    than ``PAGE_ROWS`` carry the page size the script pages them by, and
    without JavaScript the same markup is the full listing.
    """
    if not rows:
        return ""
    keep = [
        index
        for index in range(len(columns))
        if any(row[index] is not None and row[index] != "" for row in rows)
    ]
    if not keep:
        return ""
    head = "".join(f"<th>{escape(columns[index][0])}</th>" for index in keep)
    body = []
    for row in rows:
        cells = "".join(_cell(row[index], columns[index][1]) for index in keep)
        body.append(f"<tr>{cells}</tr>")
    paging = f' data-page-size="{PAGE_ROWS}"' if len(rows) > PAGE_ROWS else ""
    parts = [
        '<div class="table-block">',
        _stats_html(columns, rows, stats),
        f'<div class="scroll"><table{paging}>',
        f"<thead><tr>{head}</tr></thead>",
        "<tbody>" + "".join(body) + "</tbody>",
        "</table></div>",
        "</div>",
    ]
    return "".join(parts)


def _cell(value, formatter):
    """One table cell, carrying its raw value when the browser can sort on it.

    The rendered text is grouped and rounded for reading — "1,234", "87.0" —
    which sorts wrongly as text, so numeric cells keep the unformatted number
    in ``data-sort`` for the sorting script to use.
    """
    rendered = formatter(value)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return f"<td>{rendered}</td>"
    return f'<td data-sort="{value}">{rendered}</td>'


def _note(message):
    return f'<p class="note">{escape(message)}</p>'


def _paragraph(message):
    return f"<p>{escape(message)}</p>"


def _section_intro(name):
    """The plain-language note shown under a section's own heading."""
    explanation = SECTION_INTROS.get(name)
    return _note(explanation) if explanation else ""


def _heading(title, explanation=None, action=None):
    """A subsection heading, with a plain-language note on how to read it.

    Neighbouring tables often hold statistics that look interchangeable but
    are not — a read-weighted rate beside a mean of per-sample rates, say —
    so each heading carries a sentence or two aimed at a reader who does not
    already know how the number was computed.

    ``action`` is markup for a control belonging to the subsection — a unit
    toggle, say — laid on the same line as the title, against the right edge
    of the figure below it.
    """
    if action:
        markup = (
            f'<div class="head-row wide"><h3>{escape(title)}</h3>{action}</div>'
        )
    else:
        markup = f"<h3>{escape(title)}</h3>"
    if explanation:
        markup += _note(explanation)
    return markup


def _unit_toggle(target, unit, alt_unit, axis_title, alt_axis_title):
    """A button that swaps a figure between the two units it carries.

    The figure is drawn in ``unit`` and each of its traces holds the same
    quantity in ``alt_unit``; the script swaps the two. Without JavaScript
    there is nothing to swap, so the button is hidden rather than dead.
    """
    return (
        f'<button type="button" class="unit-toggle" data-target="{escape(target)}"'
        f' data-unit="{escape(unit)}" data-alt-unit="{escape(alt_unit)}"'
        f' data-axis="{escape(axis_title)}" data-alt-axis="{escape(alt_axis_title)}">'
        f'Switch to {escape(alt_unit)}</button>'
    )


def _duration(started, finished):
    """Format the wall-clock span between two ISO timestamps."""
    if not started or not finished:
        return None
    try:
        start = datetime.fromisoformat(str(started))
        end = datetime.fromisoformat(str(finished))
    except ValueError:
        return None
    seconds = (end - start).total_seconds()
    if seconds < 0:
        return None
    hours, remainder = divmod(int(seconds), 3600)
    minutes = remainder // 60
    if hours:
        return f"{hours} h {minutes} min"
    if minutes:
        return f"{minutes} min"
    return f"{int(seconds)} s"


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _figure_html(figure, state):
    """Serialize one Plotly figure, embedding the bundle only in the first.

    The house style is applied here rather than in each renderer, but a figure
    that set a margin, a height or the ``WIDE`` marker keeps what it asked
    for: a legend laid above the plot needs more room at the top than the
    default leaves.
    """
    import plotly.io as pio

    margin = dict(l=70, r=30, t=40, b=70)
    for side in margin:
        asked = getattr(figure.layout.margin, side)
        if asked is not None:
            margin[side] = asked
    wide = figure.layout.meta == WIDE

    figure.update_layout(
        template="simple_white",
        font=dict(family="Helvetica Neue, Helvetica, Arial, sans-serif",
                  size=13, color="#22303c"),
        margin=margin,
        paper_bgcolor="white",
        plot_bgcolor="white",
        height=figure.layout.height or FIGURE_HEIGHT,
        colorway=list(PALETTE),
    )
    include = "inline" if not state["plotly_embedded"] else False
    state["plotly_embedded"] = True
    return f'<div class="figure{" wide" if wide else ""}">' + pio.to_html(
        figure,
        include_plotlyjs=include,
        full_html=False,
        config={"displaylogo": False, "responsive": True},
    ) + "</div>"


def _serialize(blocks, state):
    """Turn a section's blocks — HTML strings and figures — into markup."""
    rendered = []
    for block in blocks:
        if block is None or block == "":
            continue
        rendered.append(block if isinstance(block, str) else _figure_html(block, state))
    return "\n".join(rendered)


# ---------------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------------

def _query(connection, statement, parameters=()):
    """Run a query, returning an empty list when the table is not there.

    A database written by an older ingest may lack a table this renderer knows
    about; that is a missing section, not a failure.
    """
    try:
        return connection.execute(statement, parameters).fetchall()
    except sqlite3.Error:
        return []


def _scalar(connection, statement, parameters=(), default=None):
    rows = _query(connection, statement, parameters)
    if not rows or rows[0][0] is None:
        return default
    return rows[0][0]


def _placeholders(values):
    return ", ".join(["?"] * len(values))


def present_sections(connection):
    """Return the section names the database actually holds data for."""
    rows = _query(
        connection,
        "SELECT section, SUM(COALESCE(rows_ingested, 0)) FROM ingest_log GROUP BY section",
    )
    return {row[0] for row in rows if row[1]}


# ---------------------------------------------------------------------------
# Section renderers
#
# Each returns a list of blocks, or an empty list when the section holds
# nothing worth showing.
# ---------------------------------------------------------------------------

def _render_preprocessing(connection):
    rows = _query(connection, """
        SELECT sample_id, reads_pre_fastp, bases_pre_fastp, reads_post_fastp,
               bases_post_fastp, host_reads, host_bases, metagenomic_reads,
               metagenomic_bases,
               singlem_fraction, nonpareil_C, nonpareil_LR, nonpareil_LRstar,
               nonpareil_diversity
        FROM sample
        ORDER BY sample_id
    """)
    if not rows:
        return []

    blocks = [
        _paragraph(
            f"{_quantity(len(rows), 'sample')} passed through quality filtering."
        ),
        _note(
            "Reads in are the raw reads; After fastp is what survived quality "
            "and adapter trimming; Host reads were then removed by mapping "
            "against the host genome, leaving the metagenomic reads. "
            "Metagenomic % is that final count as a share of the raw reads. "
            "Each average is given as reads, with the same quantity in "
            "gigabases underneath."
        ),
        _preprocessing_highlights(rows),
    ]

    table_rows = [
        (
            row["sample_id"],
            row["reads_pre_fastp"],
            row["reads_post_fastp"],
            row["host_reads"],
            row["metagenomic_reads"],
            _ratio(row["metagenomic_reads"], row["reads_pre_fastp"]),
            row["metagenomic_bases"],
        )
        for row in rows
    ]
    blocks.append(_table(
        [
            ("Sample", _text),
            ("Reads in", _integer),
            ("After fastp", _integer),
            ("Host reads", _integer),
            ("Metagenomic reads", _integer),
            ("Metagenomic %", _ONE),
            ("Metagenomic bases", _gigabases),
        ],
        table_rows,
    ))

    figure = _read_fate_figure(rows)
    if figure is not None:
        blocks.append(_heading(
            "Read fates",
            "Where each sample's reads ended up: discarded by quality "
            "filtering, matched to the host genome, or retained as "
            "metagenomic. A tall host fraction means the extraction was "
            "dominated by host DNA, which limits how much sequencing depth "
            "was available for the microbes but is not a sequencing failure. "
            "Switching to bases plots the same stack in gigabases, weighing "
            "each fate by sequence length rather than read count: the two "
            "differ once trimming has shortened the reads that survived it.",
            action=_unit_toggle(
                "read-fates", "reads", "bases", "Reads", "Gigabases"
            ),
        ))
        blocks.append('<div id="read-fates">')
        blocks.append(figure)
        blocks.append("</div>")
    blocks.extend(_microbial_blocks(rows))
    return blocks


def _preprocessing_highlights(rows):
    """Mean reads at each stage, each with the same quantity in gigabases.

    Host reads sit immediately left of the metagenomic ones, so the two halves
    of what was sequenced are read side by side.
    """
    if len(rows) < MIN_STAT_ROWS:
        return ""

    def card(label, reads, bases):
        return label, _integer(_mean(rows, reads)), _gigabases(_mean(rows, bases))

    # _mean reads one column out of each row, so the ratios are handed to it
    # as rows of their own.
    ratios = [
        (_ratio(row["metagenomic_reads"], row["reads_pre_fastp"]),) for row in rows
    ]
    return _highlights([
        card("Mean reads in", "reads_pre_fastp", "bases_pre_fastp"),
        card("Mean host reads", "host_reads", "host_bases"),
        card("Mean metagenomic reads", "metagenomic_reads", "metagenomic_bases"),
        ("Mean metagenomic %", _ONE(_mean(ratios, 0))),
    ])


def _read_fate_figure(rows):
    """Stacked read counts per sample: low quality, host, and metagenomic.

    Every trace also carries the same fate in gigabases, which the button
    beside the heading swaps in: the two orderings are not the same once
    trimming has shortened the surviving reads.
    """
    import plotly.graph_objects as go

    samples = [row["sample_id"] for row in rows]

    def difference(before, after):
        return [
            row[before] - row[after]
            if row[before] is not None and row[after] is not None else None
            for row in rows
        ]

    def column(name):
        return [row[name] for row in rows]

    def in_gigabases(values):
        return [None if v is None else v / 1e9 for v in values]

    traces = [
        ("Removed by fastp",
         difference("reads_pre_fastp", "reads_post_fastp"),
         difference("bases_pre_fastp", "bases_post_fastp")),
        ("Host", column("host_reads"), column("host_bases")),
        ("Metagenomic", column("metagenomic_reads"), column("metagenomic_bases")),
    ]
    traces = [
        (name, reads, bases) for name, reads, bases in traces
        if any(v is not None for v in reads)
    ]
    if not traces:
        return None

    # The button swaps the whole stack at once, so it is only offered when
    # every fate knows its bases: half a stack in gigabases means nothing.
    with_bases = all(
        any(v is not None for v in bases) for _, _, bases in traces
    )

    figure = go.Figure()
    for name, reads, bases in traces:
        figure.add_bar(
            name=name, x=samples, y=reads,
            meta=dict(alt=in_gigabases(bases)) if with_bases else None,
        )
    figure.update_layout(
        barmode="stack",
        xaxis_title="Sample",
        yaxis_title="Reads",
        # One bar per sample, so the plot wants every pixel of the page: the
        # legend goes above it in a single row rather than down its right side.
        legend=dict(orientation="h", traceorder="normal", title_text="",
                    yanchor="bottom", y=1.02, xanchor="left", x=0),
        margin=dict(t=60),
        meta=WIDE,
    )
    return figure


def _microbial_blocks(rows):
    """Microbial fraction and Nonpareil coverage, when they were estimated.

    Both come from optional preprocessing steps — ``--fraction`` and
    ``--nonpareil`` — so the whole subsection is left out when neither ran.
    """
    columns = ("singlem_fraction", "nonpareil_C", "nonpareil_LR",
               "nonpareil_LRstar", "nonpareil_diversity")
    if not any(row[name] is not None for row in rows for name in columns):
        return []

    blocks = [_heading(
        "Microbial fraction and coverage",
        "Two estimates made on the metagenomic reads themselves. The microbial "
        "fraction, from SingleM, is the share of the sequenced material that "
        "came from prokaryotes: a low value means most of what was sequenced "
        "was neither bacterial nor archaeal, even after host removal. "
        "Nonpareil completeness is how much of the community's sequence space "
        "the sample already covers at this depth, and the effort for near "
        "coverage is the sequencing Nonpareil projects would be needed to "
        "reach it — compare it against the effort actually spent. Nonpareil "
        "diversity is an index of community diversity estimated from read "
        "redundancy, useful for ranking samples against each other rather "
        "than as an absolute value."
    )]

    if len(rows) >= MIN_STAT_ROWS:
        blocks.append(_highlights([
            ("Mean microbial fraction", _share(_mean(rows, "singlem_fraction"))),
            ("Mean Nonpareil completeness", _share(_mean(rows, "nonpareil_C"))),
            ("Mean Nonpareil diversity", _TWO(_mean(rows, "nonpareil_diversity"))),
        ]))

    blocks.append(_table(
        [
            ("Sample", _text),
            ("Microbial fraction", _share),
            ("Nonpareil completeness", _share),
            ("Nonpareil diversity", _TWO),
            ("Effort spent", _gigabases),
            ("Effort for near coverage", _gigabases),
        ],
        [
            (
                row["sample_id"],
                row["singlem_fraction"],
                row["nonpareil_C"],
                row["nonpareil_diversity"],
                row["nonpareil_LR"],
                row["nonpareil_LRstar"],
            )
            for row in rows
        ],
    ))

    figure = _microbial_figure(rows)
    if figure is not None:
        blocks.append(figure)
    return blocks


def _microbial_figure(rows):
    """Microbial fraction beside Nonpareil completeness, per sample.

    Both are shares of one, so they belong on a single percentage axis; a
    sample that is mostly prokaryotic but poorly covered reads at a glance.
    """
    import plotly.graph_objects as go

    series = [
        ("Microbial fraction", "singlem_fraction"),
        ("Nonpareil completeness", "nonpareil_C"),
    ]
    series = [
        (name, column) for name, column in series
        if any(row[column] is not None for row in rows)
    ]
    if not series:
        return None

    def share(value):
        """Percentages of whatever scale the value arrived on, as _share does."""
        if value is None:
            return None
        return value * 100 if value <= 1 else value

    samples = [row["sample_id"] for row in rows]
    figure = go.Figure()
    for name, column in series:
        figure.add_bar(name=name, x=samples, y=[share(row[column]) for row in rows])
    figure.update_layout(
        barmode="group",
        xaxis_title="Sample",
        yaxis_title="%",
        yaxis_range=[0, 100],
        legend=dict(orientation="h", traceorder="normal", title_text="",
                    yanchor="bottom", y=1.02, xanchor="left", x=0),
        margin=dict(t=60),
        meta=WIDE,
    )
    return figure


def _hex_to_rgba(colour, alpha):
    """A palette colour as an ``rgba()`` string, for translucent Sankey links."""
    value = colour.lstrip("#")
    red, green, blue = (int(value[index:index + 2], 16) for index in (0, 2, 4))
    return f"rgba({red}, {green}, {blue}, {alpha})"


def _render_cataloging(connection):
    """Assemblies and the bins recovered from them, as two separate tables.

    They describe different objects — one row per assembly either way, but
    contiguity and mapping belong to the assembly and quality belongs to the
    bins — and reading them side by side in one wide table meant scrolling.
    """
    rows = _query(connection, """
        SELECT assembly_id, assembly_contigs, assembly_total_length,
               assembly_largest_contig, assembly_N50, assembly_gc_percent,
               mapping_rate_percent, final_bins, high_quality_bins,
               medium_quality_bins, bin_mean_completeness,
               bin_mean_contamination
        FROM assembly
        ORDER BY assembly_id
    """)
    if not rows:
        return []

    total_bins = sum(row["final_bins"] or 0 for row in rows)
    blocks = [_paragraph(
        f"{_quantity(len(rows), 'assembly', 'assemblies')} yielded "
        f"{_quantity(total_bins, 'bin')} in total."
    )]
    blocks.extend(_assembly_blocks(rows))
    blocks.extend(_mapping_rate_blocks(connection))
    blocks.extend(_bin_blocks(rows, connection))
    blocks.extend(_binner_blocks(connection))
    blocks.extend(_bin_fate_blocks(connection))
    return blocks


def _assembly_blocks(rows):
    """Contiguity and read recruitment, per assembly."""
    return [
        _heading(
            "Assemblies",
            "One row per assembly. Contigs, total length, largest contig and "
            "N50 describe contiguity: a higher N50 means the sequence sits in "
            "fewer, longer pieces. Mapping rate % is read-weighted and pooled "
            "over every sample mapped against this assembly — all mapped "
            "reads divided by all reads — so deeply sequenced samples count "
            "for more, and a low value means the assembly represents little "
            "of the sequenced material."
        ),
        _table(
            [
                ("Assembly", _text),
                ("Contigs", _integer),
                ("Total length", _integer),
                ("Largest contig", _integer),
                ("N50", _integer),
                ("GC %", _ONE),
                ("Mapping rate %", _ONE),
            ],
            [
                (
                    row["assembly_id"], row["assembly_contigs"],
                    row["assembly_total_length"], row["assembly_largest_contig"],
                    row["assembly_N50"], row["assembly_gc_percent"],
                    row["mapping_rate_percent"],
                )
                for row in rows
            ],
            stats=[
                ("Mean contigs", 1),
                ("Mean N50", 4),
                ("Mean assembly mapping rate %", 6),
            ],
        ),
    ]


def _bin_blocks(rows, connection):
    """Bin yield and quality, per assembly.

    Low-quality bins are discarded downstream, so they are counted in
    ``final_bins`` but not given a column of their own.
    """
    final_bins = sum(row["final_bins"] or 0 for row in rows)
    produced = _scalar(connection, "SELECT SUM(bin_count) FROM assembly_binner")
    share = _ratio(final_bins, produced)
    return [
        _heading(
            "Bins",
            "One row per assembly, counting the genome bins recovered from "
            "it. High and medium quality follow the usual completeness and "
            "contamination thresholds; low-quality bins are discarded "
            "downstream, so they are included in Final bins but have no "
            "column of their own. Completeness and contamination are "
            "averaged over the final bins of each assembly, so an assembly "
            "with few bins can show a high mean without having yielded much."
        ),
        _highlights([
            (
                "Bins produced by the binners",
                _integer(produced),
                "before Binette reconciled them",
            ),
            (
                "Final bins after Binette",
                _integer(final_bins),
                f"{share:.1f}% of the bins produced" if share is not None else None,
            ),
        ]),
        _table(
            [
                ("Assembly", _text),
                ("Final bins", _integer),
                ("High quality", _integer),
                ("Medium quality", _integer),
                ("Mean completeness", _ONE),
                ("Mean contamination", _ONE),
            ],
            [
                (
                    row["assembly_id"], row["final_bins"],
                    row["high_quality_bins"], row["medium_quality_bins"],
                    row["bin_mean_completeness"], row["bin_mean_contamination"],
                )
                for row in rows
            ],
            stats=[
                ("Mean bins per assembly", 1),
                ("Mean completeness", 4),
                ("Mean contamination", 5),
            ],
        ),
    ]


def _binner_blocks(connection):
    """Bins contributed by each binner, aggregated across assemblies."""
    import plotly.graph_objects as go

    rows = _query(connection, """
        SELECT binner, SUM(bin_count) AS bins, COUNT(*) AS assemblies
        FROM assembly_binner
        GROUP BY binner
        ORDER BY bins DESC
    """)
    if not rows:
        return []
    figure = go.Figure()
    figure.add_bar(
        x=[row["binner"] for row in rows],
        y=[row["bins"] for row in rows],
        marker_color=PALETTE[0],
    )
    figure.update_layout(xaxis_title="Binner", yaxis_title="Bins produced")
    return [
        _heading(
            "Bins per binner",
            "How many bins each binning tool produced, before the tools' "
            "outputs were reconciled into the final set. The tools use "
            "different strategies, so uneven counts are expected; producing "
            "more bins is not by itself a sign of a better binner, and the "
            "same genome is usually found by several of them."
        ),
        _table(
            [("Binner", _text), ("Bins", _integer), ("Assemblies", _integer)],
            [tuple(row) for row in rows],
            stats=[("Mean bins per binner", 1)],
        ),
        figure,
    ]


def _bin_fate_counts(connection):
    """Per binner, what became of the bins it produced.

    Binette records the origin of every final bin: the binner that produced it
    when it kept the bin as it was, several binners when they all produced the
    same bin, and itself when it built a bin out of contigs no single binner
    had grouped that way. Everything a binner produced that is not named by a
    final bin was replaced — by a bin from another binner, or by one Binette
    assembled in its place.
    """
    produced = {
        row["binner"]: row["bins"] or 0
        for row in _query(connection, """
            SELECT binner, SUM(bin_count) AS bins
            FROM assembly_binner
            GROUP BY binner
        """)
    }
    bins = _query(connection, """
        SELECT bin.assembly_id, bin.bin_name,
               GROUP_CONCAT(origin.binner) AS binners
        FROM assembly_bin AS bin
        LEFT JOIN assembly_bin_origin AS origin
               ON origin.assembly_id = bin.assembly_id
              AND origin.bin_name = bin.bin_name
        GROUP BY bin.assembly_id, bin.bin_name
    """)
    if not produced or not bins:
        return None

    kept_alone = {binner: 0 for binner in produced}
    kept_shared = {binner: 0 for binner in produced}
    shared_bins = 0
    built = 0
    for row in bins:
        names = [
            name for name in (row["binners"] or "").split(",")
            if name and name in produced
        ]
        if not names:
            # Either a bin Binette assembled itself, or one whose origin names
            # a binner this run did not record a count for; both are Binette's.
            built += 1
            continue
        if len(names) == 1:
            kept_alone[names[0]] += 1
            continue
        shared_bins += 1
        for name in names:
            kept_shared[name] += 1

    replaced = {
        binner: max(count - kept_alone[binner] - kept_shared[binner], 0)
        for binner, count in produced.items()
    }
    if not any(kept_alone.values()) and not any(kept_shared.values()) and not built:
        return None
    return {
        "produced": produced,
        "kept_alone": kept_alone,
        "kept_shared": kept_shared,
        "shared_bins": shared_bins,
        "built": built,
        "replaced": replaced,
    }


def _bin_fate_figure(counts):
    """The per-binner bin fates as a Sankey diagram."""
    import plotly.graph_objects as go

    binners = sorted(counts["produced"], key=lambda name: -counts["produced"][name])
    labels = list(binners) + [
        "Built by Binette",
        "Kept as the final bin",
        "Produced by several binners",
        "Replaced by another bin",
        "Final bins",
    ]
    built_index = len(binners)
    kept_index = built_index + 1
    shared_index = built_index + 2
    replaced_index = built_index + 3
    final_index = built_index + 4

    node_colours = [PALETTE[index % len(PALETTE)] for index in range(len(binners))]
    node_colours += [PALETTE[4], PALETTE[2], PALETTE[0], PALETTE[9], PALETTE[1]]

    sources, targets, values, link_colours = [], [], [], []

    def link(source, target, value, colour):
        if value:
            sources.append(source)
            targets.append(target)
            values.append(value)
            link_colours.append(_hex_to_rgba(colour, 0.35))

    for index, binner in enumerate(binners):
        link(index, kept_index, counts["kept_alone"][binner], PALETTE[2])
        link(index, shared_index, counts["kept_shared"][binner], PALETTE[0])
        link(index, replaced_index, counts["replaced"][binner], PALETTE[9])
    link(kept_index, final_index, sum(counts["kept_alone"].values()), PALETTE[2])
    link(shared_index, final_index, counts["shared_bins"], PALETTE[0])
    link(built_index, final_index, counts["built"], PALETTE[4])

    figure = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            label=labels,
            color=node_colours,
            pad=18,
            thickness=14,
            line=dict(width=0),
        ),
        link=dict(source=sources, target=targets, value=values, color=link_colours),
    ))
    figure.update_layout(height=460, margin=dict(l=10, r=10, t=20, b=20))
    return figure


def _bin_fate_blocks(connection):
    """What happened to every bin the binners produced, per binner."""
    counts = _bin_fate_counts(connection)
    if counts is None:
        return []
    binners = sorted(counts["produced"], key=lambda name: -counts["produced"][name])
    table_rows = [
        (
            binner,
            counts["produced"][binner],
            counts["kept_alone"][binner] + counts["kept_shared"][binner],
            counts["kept_shared"][binner],
            counts["replaced"][binner],
        )
        for binner in binners
    ]
    blocks = [
        _heading(
            "What became of each binner's bins",
            "Binette does not pick one binner's output over another's: it "
            "compares every bin the binners produced, keeps the best "
            "non-overlapping ones, and builds new bins where a mix of contigs "
            "scores better than anything a single binner proposed. A bin is "
            "replaced when a better-scoring bin covers the same contigs, so a "
            "binner contributing few final bins is not necessarily a poor one "
            "— it may have found the same genomes as its neighbours and lost "
            "the tie. Bins several binners produced identically are counted "
            "once for each of them on the left, which is why that stream "
            "narrows before it reaches the final bins."
        ),
        _bin_fate_figure(counts),
        _table(
            [
                ("Binner", _text),
                ("Bins produced", _integer),
                ("Kept as a final bin", _integer),
                ("Of those, also found by another binner", _integer),
                ("Replaced", _integer),
            ],
            table_rows,
        ),
    ]
    if counts["built"]:
        blocks.append(_note(
            f"Binette built {_quantity(counts['built'], 'final bin')} of its own, "
            "out of contigs that no single binner had grouped that way."
        ))
    return blocks


def _mapping_rate_blocks(connection):
    """Per-sample mapping rates against the assemblies they contributed to.

    These are reads mapped back to the assembly, not to the bins recovered
    from it, so they belong beside the assembly statistics.
    """
    rows = _query(connection, """
        SELECT sample_id,
               COUNT(*) AS assemblies,
               AVG(mapping_rate_percent) AS mean_rate,
               MIN(mapping_rate_percent) AS min_rate,
               MAX(mapping_rate_percent) AS max_rate
        FROM assembly_sample
        WHERE mapping_rate_percent IS NOT NULL
        GROUP BY sample_id
        ORDER BY sample_id
    """)
    if not rows:
        return []
    return [
        _heading(
            "Per-sample mapping rates",
            "One row per sample: how well that sample's own reads map to the "
            "assemblies it contributed to. Mean rate % averages the sample's "
            "rate across those assemblies, counting each assembly once "
            "regardless of depth, so it is not the same statistic as Mapping "
            "rate % in the assembly table above, which pools reads across "
            "samples within one assembly. Min and max show the spread; when a "
            "sample belongs to a single assembly the three columns coincide. "
            "A low value flags a sample poorly represented by the "
            "assemblies, for instance because of residual host DNA or an "
            "unusual community."
        ),
        _table(
            [
                ("Sample", _text),
                ("Assemblies", _integer),
                ("Mean rate %", _ONE),
                ("Min rate %", _ONE),
                ("Max rate %", _ONE),
            ],
            [tuple(row) for row in rows],
            stats=[("Mean sample mapping rate %", 2)],
        ),
    ]


def _ani_percent(value):
    """Read a stored ANI as a percentage, whichever way dRep wrote it."""
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number * 100 if number <= 1 else number


def _render_dereplication(connection):
    rows = _query(connection, "SELECT * FROM dereplication")
    if not rows:
        return []
    row = rows[0]

    blocks = []
    threshold = _ani_percent(row["dereplication_ani"])
    if threshold is not None:
        blocks.append(_paragraph(
            f"Bins were collapsed into one representative genome wherever they "
            f"shared at least {threshold:g}% average nucleotide identity."
        ))

    before, after = row["input_bin_number"], row["output_bin_number"]
    blocks.extend(_dereplication_yield_blocks(before, after))
    blocks.append(_table(
        [
            ("Stage", _text),
            ("Bins", _integer),
            ("Mean completeness", _ONE),
            ("Mean contamination", _ONE),
        ],
        [
            ("Before dereplication", row["input_bin_number"],
             row["input_bin_completeness"], row["input_bin_contamination"]),
            ("After dereplication", row["output_bin_number"],
             row["output_bin_completeness"], row["output_bin_contamination"]),
        ],
    ))
    blocks.extend(_dereplication_clustering_blocks(connection, threshold))
    return blocks


def _dereplication_yield_blocks(before, after):
    """How much of the input survived, as highlights and as one stacked bar.

    One bar rather than two: the question is what share of the input the
    catalogue kept, and a part-to-whole split answers it directly, where two
    bars side by side leave the reader to do the division.
    """
    import plotly.graph_objects as go

    if before is None or after is None:
        return []
    collapsed = max(before - after, 0)
    share = _ratio(after, before)
    blocks = [_highlights([
        ("Bins dereplicated", _integer(before)),
        ("Representative genomes", _integer(after)),
        (
            "Bins retained",
            f"{share:.1f}%" if share is not None else None,
            f"{collapsed:,} collapsed into a representative",
        ),
    ])]
    if not before:
        return blocks

    figure = go.Figure()
    for label, value, colour in (
        ("Retained as a representative", after, PALETTE[0]),
        ("Collapsed into a representative", collapsed, PALETTE[3]),
    ):
        figure.add_bar(
            y=["Bins"],
            x=[value],
            name=label,
            orientation="h",
            marker_color=colour,
            text=[f"{value:,}"],
            textposition="inside",
            insidetextanchor="middle",
        )
    figure.update_layout(
        barmode="stack",
        height=210,
        margin=dict(l=30, r=30, t=60, b=40),
        xaxis_title="Bins",
        yaxis=dict(showticklabels=False),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    blocks.append(figure)
    return blocks


# dRep clusters bins at 90% MASH identity before it compares anything with ANI,
# so a bin whose nearest MASH neighbour is further away than this was never a
# candidate for dereplication: nothing was close enough to collapse it into.
PRIMARY_MASH_DISTANCE = 0.1

# Half-percent bands from 100% down to this floor, then one bucket for the rest.
ANI_BAND_FLOOR = 95.0
ANI_BAND_WIDTH = 0.5

# The histogram is drawn finer than the bands are counted, so that a pile-up
# against the threshold is visible rather than averaged into a half-point bar.
ANI_BIN_WIDTH = 0.25


def _ani_bands():
    """The ``(lower, upper)`` band edges, widest identity first."""
    bands = []
    upper = 100.0
    while upper > ANI_BAND_FLOOR + 1e-9:
        bands.append((upper - ANI_BAND_WIDTH, upper))
        upper -= ANI_BAND_WIDTH
    return bands


def _dereplication_clustering_blocks(connection, threshold):
    """How the ANI threshold acted on the pairs it was actually applied to.

    dRep only compares bins that already look alike: it clusters them by MASH
    first and runs ANI within those clusters. So the bins whose fate the
    threshold decided are the ones with a near MASH neighbour, and the pairs
    worth looking at are the ones inside a cluster. A run whose comparisons all
    sit far from the threshold was insensitive to its exact value; a dense band
    of pairs just either side of it means the opposite.
    """
    total = _scalar(connection, "SELECT COUNT(*) FROM genome_cluster", default=0)
    if not total:
        return []

    with_neighbour = _scalar(
        connection,
        "SELECT COUNT(*) FROM genome_cluster "
        "WHERE nearest_mash_distance IS NOT NULL AND nearest_mash_distance <= ?",
        (PRIMARY_MASH_DISTANCE,),
        default=0,
    )
    pairs = _query(connection, """
        SELECT comparison.ani AS ani,
               CASE
                   WHEN first.secondary_cluster IS NOT NULL
                    AND first.secondary_cluster = second.secondary_cluster
                   THEN 1 ELSE 0
               END AS same_cluster
        FROM genome_comparison AS comparison
        LEFT JOIN genome_cluster AS first ON first.genome = comparison.genome_a
        LEFT JOIN genome_cluster AS second ON second.genome = comparison.genome_b
        WHERE comparison.ani IS NOT NULL
    """)

    neighbour_share = _ratio(with_neighbour, total)
    blocks = [
        _heading(
            "How the identity threshold acted",
            "dRep never compares every bin with every other one: it groups "
            "them by MASH sketch first and computes identity only inside those "
            "groups. A bin with no MASH neighbour within "
            f"{PRIMARY_MASH_DISTANCE:.0%} distance was therefore never a "
            "candidate for collapsing, whatever the threshold had been set to. "
            "The comparisons below are the ones the threshold was applied to, "
            "each pair counted once, at the average of the two directions dRep "
            "clusters on."
        ),
        _highlights([
            ("Bins with a MASH neighbour", _integer(with_neighbour),
             f"within {PRIMARY_MASH_DISTANCE:.0%} distance — "
             f"{neighbour_share:.1f}% of {total:,} bins"
             if neighbour_share is not None else None),
            ("Pairwise comparisons", _integer(len(pairs)) if pairs else None),
        ]),
    ]
    if not pairs:
        return blocks

    identities = [
        (_ani_percent(pair["ani"]), bool(pair["same_cluster"]))
        for pair in pairs
    ]
    identities = [item for item in identities if item[0] is not None]
    if not identities:
        return blocks
    blocks.append(_ani_histogram(identities, threshold))
    blocks.append(_ani_band_table(identities))
    return blocks


def _ani_histogram(identities, threshold):
    """Pairwise identities, split by whether the pair was collapsed."""
    import plotly.graph_objects as go

    # Snap the first bin edge onto a multiple of the bin width, so the bars
    # fall on round identities instead of wherever the lowest pair happens to be.
    lowest = min(value for value, _ in identities)
    start = max(int(lowest / ANI_BIN_WIDTH) * ANI_BIN_WIDTH, 0.0)
    figure = go.Figure()
    for label, keep, colour in (
        ("Collapsed into one genome", True, PALETTE[0]),
        ("Kept as separate genomes", False, PALETTE[3]),
    ):
        values = [value for value, same in identities if same is keep]
        if not values:
            continue
        figure.add_histogram(
            x=values,
            name=label,
            marker_color=colour,
            xbins=dict(start=start, end=100.0 + ANI_BIN_WIDTH, size=ANI_BIN_WIDTH),
        )
    if threshold is not None:
        figure.add_vline(
            x=threshold,
            line_dash="dash",
            line_color="#22303c",
            annotation_text=f"{threshold:g}% threshold",
            annotation_position="top left",
        )
    figure.update_layout(
        barmode="stack",
        xaxis_title="Average nucleotide identity between the pair (%)",
        yaxis_title="Pairs of bins",
        margin=dict(l=70, r=30, t=70, b=70),
        legend=dict(orientation="h", yanchor="bottom", y=1.06, x=0),
    )
    return figure


def _ani_band_table(identities):
    """The same distribution as counts, in half-percent bands."""
    bands = _ani_bands()
    counts = [[0, 0] for _ in bands]
    below = [0, 0]
    for value, same in identities:
        index = 1 if same else 0
        for position, (lower, upper) in enumerate(bands):
            # The top band is closed at 100 so that an identical pair lands in it.
            if value > lower and (value <= upper or position == 0):
                counts[position][index] += 1
                break
        else:
            below[index] += 1

    rows = [
        (f"{upper:g} – {lower:g}%", separate + collapsed, collapsed, separate)
        for (lower, upper), (separate, collapsed) in zip(bands, counts)
    ]
    rows.append((
        f"below {ANI_BAND_FLOOR:g}%", below[0] + below[1], below[1], below[0],
    ))
    # Empty bands between two populated ones are the point of the table — an
    # empty band next to the threshold says the decision was clear-cut — so
    # only the leading and trailing runs of zeroes are trimmed.
    populated = [index for index, row in enumerate(rows) if row[1]]
    if not populated:
        return ""
    rows = rows[populated[0]:populated[-1] + 1]
    return _table(
        [
            ("Identity band", _text),
            ("Pairs", _integer),
            ("Collapsed", _integer),
            ("Kept separate", _integer),
        ],
        rows,
    )


def _render_profiling(connection):
    blocks = []
    genomes = _scalar(connection, "SELECT COUNT(*) FROM genome", default=0)
    quality = _query(connection, """
        SELECT COUNT(*) AS genomes, AVG(completeness) AS completeness,
               AVG(contamination) AS contamination, AVG(size) AS size,
               AVG(n50) AS n50, SUM(cluster_members) AS members
        FROM genome
    """)
    if genomes:
        row = quality[0]
        blocks.append(_paragraph(
            f"The catalogue holds {_quantity(genomes, 'dereplicated genome')}."
        ))
        blocks.append(_highlights([
            ("Genomes", _integer(row["genomes"])),
            ("Mean completeness %", _ONE(row["completeness"])),
            ("Mean contamination %", _ONE(row["contamination"])),
            ("Mean size", _integer(row["size"])),
            ("Mean N50", _integer(row["n50"])),
            ("Clustered members", _integer(row["members"])),
        ]))
        blocks.extend(_genome_quality_blocks(connection))

    blocks.extend(_abundance_blocks(connection))
    return blocks


def _genome_quality_blocks(connection):
    """Completeness against contamination, the usual MAG quality view."""
    import plotly.graph_objects as go

    rows = _query(connection, """
        SELECT genome_id, completeness, contamination, size
        FROM genome
        WHERE completeness IS NOT NULL AND contamination IS NOT NULL
    """)
    if not rows:
        return []
    figure = go.Figure()
    figure.add_scatter(
        x=[row["completeness"] for row in rows],
        y=[row["contamination"] for row in rows],
        mode="markers",
        marker=dict(color=PALETTE[0], size=7, opacity=0.7,
                    line=dict(width=0.5, color="#ffffff")),
        text=[row["genome_id"] for row in rows],
        hovertemplate="%{text}<br>completeness %{x:.1f}%<br>contamination %{y:.1f}%<extra></extra>",
        showlegend=False,
    )
    figure.update_layout(
        xaxis_title="Completeness (%)",
        yaxis_title="Contamination (%)",
    )
    return [
        _heading(
            "Genome quality",
            "Each point is one dereplicated genome, placed by estimated "
            "completeness (horizontal) against estimated contamination "
            "(vertical). The desirable corner is the bottom right: nearly "
            "complete and nearly clean. Points high on the vertical axis "
            "likely mix sequence from more than one organism, and are "
            "usually filtered out or treated with caution downstream."
        ),
        figure,
    ]


def _abundance_blocks(connection):
    """Relative abundance of the most abundant genomes across samples."""
    import plotly.graph_objects as go

    # The metagenomic read count comes from preprocessing, so it is only there
    # when that section was ingested too; the join is left outer for that reason
    # and the column drops out of the table on its own when it is empty.
    totals = _query(connection, """
        SELECT c.sample_id AS sample_id,
               SUM(c.read_count) AS reads,
               MAX(s.metagenomic_reads) AS metagenomic_reads
        FROM genome_count AS c
        LEFT JOIN sample AS s ON s.sample_id = c.sample_id
        WHERE c.read_count IS NOT NULL
        GROUP BY c.sample_id
        ORDER BY c.sample_id
    """)
    if not totals:
        return []
    sample_totals = {row["sample_id"]: row["reads"] for row in totals}

    top = _query(connection, """
        SELECT genome_id, SUM(read_count) AS reads
        FROM genome_count
        WHERE read_count IS NOT NULL
        GROUP BY genome_id
        ORDER BY reads DESC
        LIMIT ?
    """, (TOP_GENOMES,))
    if not top:
        return []
    genome_ids = [row["genome_id"] for row in top]

    cells = _query(connection, f"""
        SELECT genome_id, sample_id, read_count
        FROM genome_count
        WHERE genome_id IN ({_placeholders(genome_ids)}) AND read_count IS NOT NULL
    """, tuple(genome_ids))
    lookup = {(row["genome_id"], row["sample_id"]): row["read_count"] for row in cells}

    samples = [row["sample_id"] for row in totals]
    matrix = [
        [
            _ratio(lookup.get((genome_id, sample)), sample_totals.get(sample))
            for sample in samples
        ]
        for genome_id in genome_ids
    ]

    figure = go.Figure()
    figure.add_heatmap(
        z=matrix,
        x=samples,
        y=genome_ids,
        colorscale="Teal",
        colorbar=dict(title="% of sample"),
        hovertemplate="%{y}<br>%{x}<br>%{z:.2f}% of sample<extra></extra>",
    )
    figure.update_layout(
        xaxis_title="Sample",
        yaxis_title="Genome",
        height=max(FIGURE_HEIGHT, 18 * len(genome_ids) + 140),
    )

    total_genomes = _scalar(
        connection, "SELECT COUNT(DISTINCT genome_id) FROM genome_count", default=0
    )
    blocks = [_heading(
        "Genome abundance",
        "How much of each sample each genome accounts for, as a percentage "
        "of that sample's mapped reads rather than a raw count, so samples "
        "sequenced to different depths can be compared. Darker cells are "
        "more abundant. These are relative abundances: they say nothing "
        "about the absolute amount of microbial material in a sample."
    )]
    if total_genomes > len(genome_ids):
        blocks.append(_note(
            f"Showing the {len(genome_ids)} most abundant of {total_genomes:,} genomes, "
            f"as a percentage of each sample's mapped reads."
        ))
    blocks.append(figure)
    # Mapped reads on their own say nothing about how well the catalogue
    # represents a sample: the same count is a good result for a shallow
    # sample and a poor one for a deep sample, so the metagenomic reads the
    # mapping started from and the share of them that landed sit beside it.
    blocks.append(_note(
        "Mapped reads are the reads assigned to a catalogue genome; "
        "metagenomic reads are what preprocessing handed the mapper, after "
        "low-quality and host reads were removed. Mapped % is the first as a "
        "share of the second, and is how much of the sample the catalogue "
        "accounts for — a low value means most of the sample matched no "
        "genome in the catalogue."
    ))
    blocks.append(_table(
        [
            ("Sample", _text),
            ("Mapped reads", _integer),
            ("Metagenomic reads", _integer),
            ("Mapped %", _ONE),
        ],
        [
            (
                row["sample_id"],
                row["reads"],
                row["metagenomic_reads"],
                _ratio(row["reads"], row["metagenomic_reads"]),
            )
            for row in totals
        ],
        stats=[
            ("Mean mapped reads", 1),
            ("Mean metagenomic reads", 2),
            ("Mean mapped %", 3),
        ],
    ))
    return blocks


def _render_taxonomy(connection):
    import plotly.graph_objects as go

    classified = _scalar(connection, "SELECT COUNT(*) FROM genome_taxonomy", default=0)
    if not classified:
        return []

    rank_rows = []
    for rank in TAXONOMIC_RANKS:
        row = _query(connection, f"""
            SELECT COUNT(DISTINCT "{rank}") AS taxa,
                   SUM(CASE WHEN "{rank}" IS NULL OR "{rank}" = '' THEN 1 ELSE 0 END) AS unclassified
            FROM genome_taxonomy
        """)
        if not row:
            continue
        unclassified = row[0]["unclassified"] or 0
        rank_rows.append((
            rank.capitalize(), row[0]["taxa"], classified - unclassified, unclassified
        ))

    blocks = [
        _paragraph(
            f"GTDB classifications are available for {_quantity(classified, 'genome')}."
        ),
        _note(
            "One row per taxonomic rank, from broad to fine. Distinct taxa is "
            "how many different groups the genomes fall into at that rank, so "
            "it grows towards the species end; genomes unclassified grows the "
            "same way, because assigning a genome to a species needs a closer "
            "reference match than assigning it to a phylum."
        ),
        _table(
            [
                ("Rank", _text),
                ("Distinct taxa", _integer),
                ("Genomes classified", _integer),
                ("Genomes unclassified", _integer),
            ],
            rank_rows,
        ),
    ]

    colours = _phylum_colours(connection)
    blocks.extend(_phylogeny_blocks(connection, colours))
    blocks.extend(_genome_lineage_blocks(connection))

    phyla = _query(connection, """
        SELECT COALESCE(NULLIF(phylum, ''), 'Unclassified') AS name, COUNT(*) AS genomes
        FROM genome_taxonomy
        GROUP BY name
        ORDER BY genomes DESC
        LIMIT ?
    """, (TOP_CATEGORIES,))
    if phyla:
        figure = go.Figure()
        figure.add_bar(
            x=[row["genomes"] for row in phyla],
            y=[row["name"] for row in phyla],
            orientation="h",
            marker_color=[
                colours.get(row["name"], OTHER_COLOUR) for row in phyla
            ],
        )
        figure.update_layout(
            xaxis_title="Genomes",
            yaxis_title="",
            yaxis=dict(autorange="reversed"),
            height=max(FIGURE_HEIGHT, 22 * len(phyla) + 140),
        )
        blocks.append(_heading(
            "Genomes per phylum",
            "How many catalogue genomes fall in each phylum. This counts "
            "distinct genomes recovered, not how abundant they are: a phylum "
            "with many genomes can still be rare in the samples, and an "
            "abundant one can be represented by a single genome. Each phylum "
            "keeps the colour it has in the composition chart below; the ones "
            "that chart pools as Other are drawn in its grey."
        ))
        blocks.append(figure)

    blocks.extend(_composition_blocks(connection, colours))
    return blocks


# ---------------------------------------------------------------------------
# Circular phylogeny
# ---------------------------------------------------------------------------
#
# The tree is drawn as inline SVG rather than handed to Plotly. It has no zoom,
# pan or hover behaviour to gain from a plotting library, it has to stay legible
# when the panel it lives in is hidden and later shown — which is exactly what
# Plotly cannot do without being told to resize — and drawing it directly keeps
# a picture of several hundred genomes down to a few tens of kilobytes.
#
# Every length below is in the SVG's own units, on a square canvas the figure
# scales to whatever width its column gives it.

TREE_CANVAS = 1000
TREE_CENTRE = TREE_CANVAS / 2
# Where the root sits and where the tips do. The root is not pushed all the way
# in: a circular phylogram with no hole crowds every deep split into a few
# pixels.
TREE_INNER_RADIUS = 60
TREE_TIP_RADIUS = 300
# A wedge left undrawn at the top, so each ring can be named where it ends
# instead of in a key the reader has to match up by colour.
TREE_GAP_DEGREES = 34
TREE_RING_THICKNESS = 15
TREE_RING_GAP = 5
TREE_FIRST_RING_RADIUS = 312
TREE_BRANCH_COLOUR = "#a9b6c1"
TREE_TIP_COLOUR = "#8fa2b0"
# The dotted run from where a branch ends to where the rings begin.
TREE_LEADER_COLOUR = "#dde3e8"
# The empty part of a proportional ring — drawn so a low value reads as a short
# bar rather than as a missing one.
TREE_TRACK_COLOUR = "#eef1f4"
# A leaf the ring has no value for at all: left blank, and named as such in the
# key, so absent is visibly different from zero.
TREE_MISSING_COLOUR = "#ffffff"

# The proportional rings, in the order they are drawn outwards from the tree.
# Each names the ``genome`` column it reads, the two ends of its colour ramp,
# and how its value is turned into the fraction of the ring that is filled.
# Adding a ring is adding an entry here — the layout, the key and the labels in
# the gap all follow from it.
TREE_RINGS = (
    {
        "key": "size",
        "label": "Genome size",
        "column": "size",
        "ramp": ("#d9e6e7", "#3d7d80"),
    },
    {
        "key": "completeness",
        "label": "Completeness",
        "column": "completeness",
        "ramp": ("#e1ebdb", "#4e7a3f"),
        # A CheckM2 estimate is already a percentage of a whole genome, so the
        # ring is filled against 100 rather than against the best genome here:
        # a catalogue whose best genome is 82% complete should look like one.
        "maximum": 100.0,
    },
    {
        "key": "contamination",
        "label": "Contamination",
        "column": "contamination",
        "ramp": ("#f2ddd9", "#a8394a"),
        # Scaled against the worst genome, but never against less than the 5%
        # that separates high- from medium-quality MAGs: a clean catalogue
        # should not have its cleanest genomes redrawn as contaminated ones.
        "floor": 5.0,
    },
)

# Past this many tips the names no longer fit round the circle, so the tree is
# drawn without them and the lineage table below carries the naming.
TREE_MAX_LABELS = 110
TREE_LABEL_CHARS = 26


def _phylogeny_blocks(connection, colours):
    """The GTDB-Tk placement tree, drawn as a circular phylogram with rings.

    One figure per domain, because GTDB-Tk places bacteria and archaea on two
    unrelated reference trees and there is no branch joining them. The rings a
    figure carries are whatever the database can fill: the phylum ring needs
    only the taxonomy, the other three come from the CheckM2 columns the
    profiling section ingests, and each is left out when its column is empty.
    """
    trees = _query(
        connection,
        "SELECT domain, newick, tip_count FROM genome_tree ORDER BY tip_count DESC",
    )
    if not trees:
        return []

    phyla = {
        row["genome_id"]: row["phylum"]
        for row in _query(connection, """
            SELECT genome_id, COALESCE(NULLIF(phylum, ''), 'Unclassified') AS phylum
            FROM genome_taxonomy
        """)
    }
    quality = {
        row["genome_id"]: row
        for row in _query(
            connection,
            "SELECT genome_id, size, completeness, contamination FROM genome",
        )
    }

    blocks = []
    for tree in trees:
        figure = _phylogeny_figure(
            tree["newick"], tree["domain"], phyla, quality, colours
        )
        if figure:
            blocks.append(figure)
    if not blocks:
        return []

    return [
        _heading(
            "Phylogeny",
            "The GTDB-Tk placement tree, pruned back to the catalogue genomes "
            "and drawn from the centre outwards: branch lengths are substitutions "
            "per site, so two tips that sit close together are close relatives. "
            "The rings around it describe each genome — its phylum in the colours "
            "the rest of this section uses, then how large, how complete and how "
            "contaminated it is — so a whole clade of small or half-recovered "
            "genomes shows up as a run of similar arcs rather than as scattered "
            "rows in a table. Point at a genome to read its figures."
        ),
    ] + blocks


def _phylogeny_figure(newick, domain, phyla, quality, colours):
    """One domain's tree, as a ``.figure`` block holding the SVG and its key."""
    from drakkar.report.newick import NewickError, parse

    try:
        root = parse(newick)
    except NewickError:
        return ""
    leaves = root.leaves()
    if len(leaves) < 2:
        return ""

    lookup = _genome_lookup(phyla, quality)
    # Which of the proportional rings this database can actually fill.
    rings = []
    for spec in TREE_RINGS:
        values = [
            _tree_value(leaf.name, lookup, quality, spec["column"]) for leaf in leaves
        ]
        present = [value for value in values if value is not None]
        if not present:
            continue
        maximum = spec.get("maximum") or max(
            max(present), spec.get("floor", 0.0)
        )
        if not maximum:
            continue
        rings.append(dict(spec, values=values, maximum=maximum))

    angles, radii = _tree_layout(root)
    parts = [_tree_branches(root, angles, radii)]

    span = math.radians(360 - TREE_GAP_DEGREES)
    slot = span / len(leaves)
    # A hair of white between neighbouring arcs, so a run of same-coloured
    # genomes still reads as several genomes and not as one block.
    pad = min(0.004, slot * 0.12)
    labelled = len(leaves) <= TREE_MAX_LABELS
    # A tree whose tips no table can be joined to gets no phylum ring at all,
    # rather than a ring of grey with nothing behind it.
    named = any(_tree_phylum(leaf.name, lookup, phyla) for leaf in leaves)
    label_radius = _tree_ring_radius(len(rings))[1] + 10
    font = 11 if len(leaves) <= 40 else 9 if len(leaves) <= 70 else 8

    # How far the longest tip name reaches past the outer ring, so the canvas
    # is cropped to the picture instead of to the longest name it could have
    # held. The per-character width is an over-estimate for Helvetica: a name
    # cropped at the edge of the canvas is worse than a little air around it,
    # and a reader without the page's fonts gets a wider substitute.
    widest = max(
        (len(_tree_text(_tree_genome(leaf.name, lookup) or leaf.name))
         for leaf in leaves),
        default=0,
    )
    reach = label_radius + (widest * font * 0.65 + 12 if labelled else 8)

    for index, leaf in enumerate(leaves):
        angle = angles[id(leaf)]
        start, end = angle - slot / 2 + pad, angle + slot / 2 - pad
        segments = []
        # A tip is named by the genome id the tables use, so a name read off
        # the circle can be found in the lineage table underneath it.
        name = _tree_genome(leaf.name, lookup) or leaf.name

        phylum = phyla.get(name)
        if phylum is not None:
            inner, outer = _tree_ring_radius(0)
            segments.append(_tree_arc(
                inner, outer, start, end, colours.get(phylum, OTHER_COLOUR)
            ))
        for position, ring in enumerate(rings, start=1):
            inner, outer = _tree_ring_radius(position)
            value = ring["values"][index]
            if value is None:
                segments.append(_tree_arc(
                    inner, outer, start, end, TREE_MISSING_COLOUR
                ))
                continue
            fraction = max(0.0, min(1.0, value / ring["maximum"]))
            segments.append(_tree_arc(inner, outer, start, end, TREE_TRACK_COLOUR))
            if fraction > 0:
                segments.append(_tree_arc(
                    inner,
                    inner + (outer - inner) * fraction,
                    start, end,
                    _mix_colours(ring["ramp"][0], ring["ramp"][1], fraction),
                ))
        if labelled:
            segments.append(_tree_label(name, angle, label_radius, font))
        parts.append(
            f'<g><title>{escape(_tree_tooltip(name, phylum, rings, index))}</title>'
            + "".join(segments) + "</g>"
        )

    parts.append(_tree_ring_labels(rings, named))
    parts.append(_tree_centre_label(len(leaves), domain))

    box = max(TREE_TIP_RADIUS, min(reach, TREE_CENTRE))
    view = f"{TREE_CENTRE - box:.0f} {TREE_CENTRE - box:.0f} {2 * box:.0f} {2 * box:.0f}"
    svg = (
        f'<svg class="tree" viewBox="{view}" role="img" '
        f'aria-label="Circular phylogeny of {escape(str(domain))}" '
        'xmlns="http://www.w3.org/2000/svg">'
        f'<g transform="translate({TREE_CENTRE:.0f},{TREE_CENTRE:.0f})">'
        + "".join(parts) +
        "</g></svg>"
    )
    return (
        '<div class="figure tree-figure">' + svg
        + _tree_key(leaves, lookup, phyla, colours, rings) + "</div>"
    )


def _genome_lookup(phyla, quality):
    """Map every name a tip might carry to the genome id the tables use.

    GTDB-Tk is given the MAG names, so the tips normally match ``genome_id``
    exactly. A catalogue whose file names reached GTDB-Tk instead carries the
    FASTA suffix on the tip, so the stripped name is registered as well.
    """
    lookup = {}
    for genome_id in set(phyla) | set(quality):
        lookup[genome_id] = genome_id
        stripped = re.sub(r"\.(fa|fna|fasta|fa\.gz|fna\.gz|fasta\.gz)$", "", genome_id)
        lookup.setdefault(stripped, genome_id)
    return lookup


def _tree_genome(name, lookup):
    """The genome id a tip stands for, or None when nothing matches it."""
    if name in lookup:
        return lookup[name]
    stripped = re.sub(r"\.(fa|fna|fasta|fa\.gz|fna\.gz|fasta\.gz)$", "", name)
    return lookup.get(stripped)


def _tree_phylum(name, lookup, phyla):
    genome_id = _tree_genome(name, lookup)
    return phyla.get(genome_id) if genome_id else None


def _tree_value(name, lookup, quality, column):
    genome_id = _tree_genome(name, lookup)
    row = quality.get(genome_id) if genome_id else None
    return row[column] if row is not None else None


def _tree_layout(root):
    """Angles and radii for every node of the circular phylogram.

    Tips are spread evenly round the circle in the order the Newick file lists
    them, which is the order that keeps sister clades adjacent, and an internal
    node sits at the mean angle of its children. Radius is cumulative branch
    length from the root, so the picture is a phylogram; a tree written without
    branch lengths — GTDB-Tk emits one when placement was by topology alone —
    falls back to node depth, which draws the same topology as a cladogram.
    """
    angles = {}
    radii = {}
    leaves = root.leaves()
    span = math.radians(360 - TREE_GAP_DEGREES)
    for index, leaf in enumerate(leaves):
        angles[id(leaf)] = (index + 0.5) / len(leaves) * span

    # Parents before children, so a node's distance from the root is known by
    # the time its children are reached; reversed, the same order puts every
    # child before its parent, which is what averaging angles upwards needs.
    order = list(root.walk())
    distance = {id(root): root.length}
    depth = {id(root): 0}
    for node in order:
        for child in node.children:
            distance[id(child)] = distance[id(node)] + child.length
            depth[id(child)] = depth[id(node)] + 1
    for node in reversed(order):
        if node.children:
            angles[id(node)] = (
                sum(angles[id(child)] for child in node.children) / len(node.children)
            )

    reach = max(distance[id(leaf)] for leaf in leaves)
    if reach > 0:
        scale = {key: value / reach for key, value in distance.items()}
    else:
        deepest = max(depth[id(leaf)] for leaf in leaves) or 1
        scale = {key: value / deepest for key, value in depth.items()}
    for key, fraction in scale.items():
        radii[key] = TREE_INNER_RADIUS + fraction * (TREE_TIP_RADIUS - TREE_INNER_RADIUS)
    return angles, radii


def _tree_branches(root, angles, radii):
    """Every branch of the tree: a radial line per node, an arc per split.

    A tip whose branch ends short of the rings — most of them, in a phylogram
    — is carried the rest of the way by a leader, so the reader can follow a
    genome from its place in the tree out to the arcs describing it without
    the branch lengths being falsified to reach.
    """
    lines = []
    arcs = []
    tips = []
    leaders = []
    for node in root.walk():
        angle = angles[id(node)]
        radius = radii[id(node)]
        if node.children:
            children = sorted(angles[id(child)] for child in node.children)
            arcs.append(_tree_arc_path(radius, children[0], children[-1]))
            for child in node.children:
                lines.append(_tree_radial(angle_of=angles[id(child)],
                                          inner=radius, outer=radii[id(child)]))
        else:
            tips.append(
                f'<circle cx="{_x(radius, angle):.1f}" '
                f'cy="{_y(radius, angle):.1f}" r="1.6"/>'
            )
            if TREE_TIP_RADIUS - radius > 1:
                leaders.append(_tree_radial(
                    angle_of=angle, inner=radius, outer=TREE_TIP_RADIUS
                ))
    return (
        f'<g fill="none" stroke="{TREE_LEADER_COLOUR}" stroke-width="0.8" '
        'stroke-dasharray="1.5 3">' + "".join(leaders) + "</g>"
        + f'<g fill="none" stroke="{TREE_BRANCH_COLOUR}" stroke-width="1.3">'
        + "".join(arcs) + "".join(lines) + "</g>"
        + f'<g fill="{TREE_TIP_COLOUR}">' + "".join(tips) + "</g>"
    )


def _x(radius, angle):
    """Angles run clockwise from the top, the way a circular tree is read."""
    return radius * math.sin(angle)


def _y(radius, angle):
    return -radius * math.cos(angle)


def _tree_radial(angle_of, inner, outer):
    return (
        f'<path d="M{_x(inner, angle_of):.1f},{_y(inner, angle_of):.1f}'
        f'L{_x(outer, angle_of):.1f},{_y(outer, angle_of):.1f}"/>'
    )


def _tree_arc_path(radius, start, end):
    """A stroked arc at one radius, joining the children of a split."""
    if end - start < 1e-9:
        return ""
    large = 1 if end - start > math.pi else 0
    return (
        f'<path d="M{_x(radius, start):.1f},{_y(radius, start):.1f}'
        f'A{radius:.1f},{radius:.1f} 0 {large} 1 '
        f'{_x(radius, end):.1f},{_y(radius, end):.1f}"/>'
    )


def _tree_arc(inner, outer, start, end, fill):
    """A filled annulus sector — one genome's slice of one ring."""
    if end - start < 1e-9 or outer - inner < 1e-9:
        return ""
    large = 1 if end - start > math.pi else 0
    return (
        f'<path d="M{_x(outer, start):.1f},{_y(outer, start):.1f}'
        f'A{outer:.1f},{outer:.1f} 0 {large} 1 {_x(outer, end):.1f},{_y(outer, end):.1f}'
        f'L{_x(inner, end):.1f},{_y(inner, end):.1f}'
        f'A{inner:.1f},{inner:.1f} 0 {large} 0 {_x(inner, start):.1f},{_y(inner, start):.1f}'
        f'Z" fill="{fill}"/>'
    )


def _tree_ring_radius(position):
    """The inner and outer radius of the ring at ``position``, counting out."""
    inner = TREE_FIRST_RING_RADIUS + position * (TREE_RING_THICKNESS + TREE_RING_GAP)
    return inner, inner + TREE_RING_THICKNESS


def _tree_text(name):
    """A tip name cut to what fits beside the circle, ellipsis and all."""
    if len(name) <= TREE_LABEL_CHARS:
        return name
    return name[: TREE_LABEL_CHARS - 1] + "\u2026"


def _tree_label(name, angle, radius, font):
    """A tip name, laid along its own radius and kept the right way up."""
    degrees = math.degrees(angle)
    flipped = degrees > 180
    text = _tree_text(name)
    turn = " rotate(180)" if flipped else ""
    return (
        f'<text transform="rotate({degrees - 90:.2f}) translate({radius:.0f},0){turn}" '
        f'text-anchor="{"end" if flipped else "start"}" '
        f'dx="{-4 if flipped else 4}" dominant-baseline="middle" '
        f'font-size="{font}" fill="#5a6a78">{escape(text)}</text>'
    )


def _tree_ring_labels(rings, has_phylum):
    """The ring names, written into the wedge left empty at the top."""
    items = []
    if has_phylum:
        items.append((0, "Phylum", "#6b7a88"))
    for position, ring in enumerate(rings, start=1):
        items.append((position, ring["label"], ring["ramp"][1]))
    labels = []
    for position, label, colour in items:
        inner, outer = _tree_ring_radius(position)
        labels.append(
            f'<text x="-6" y="{-(inner + outer) / 2:.1f}" text-anchor="end" '
            f'dominant-baseline="middle" font-size="11" fill="{colour}">'
            f'{escape(label)}</text>'
        )
    return "".join(labels)


def _tree_centre_label(tips, domain):
    name = str(domain).capitalize()
    return (
        f'<text x="0" y="-6" text-anchor="middle" dominant-baseline="middle" '
        f'font-size="20" fill="#8fa2b0">{tips}</text>'
        f'<text x="0" y="14" text-anchor="middle" dominant-baseline="middle" '
        f'font-size="11" fill="#b0bcc6">{escape(name)}</text>'
    )


def _tree_tooltip(name, phylum, rings, index):
    parts = [name]
    if phylum:
        parts.append(f"Phylum: {phylum}")
    for ring in rings:
        value = ring["values"][index]
        if value is None:
            continue
        if ring["key"] == "size":
            parts.append(f"Genome size: {value / 1e6:.2f} Mbp")
        elif ring["key"] == "completeness":
            parts.append(f"Completeness: {value:.1f}%")
        else:
            parts.append(f"{ring['label']}: {value:.2f}%")
    return "\n".join(parts)


def _tree_key(leaves, lookup, phyla, colours, rings):
    """What the colours mean, under the tree they are used in.

    The phylum swatches list only the phyla this tree holds, in the order the
    rest of the section colours them, so the key doubles as the reader's route
    from an arc back to a name.
    """
    present = []
    pooled = False
    for leaf in leaves:
        phylum = _tree_phylum(leaf.name, lookup, phyla)
        if phylum is None:
            continue
        if phylum in colours:
            if phylum not in present:
                present.append(phylum)
        else:
            pooled = True
    order = list(colours)
    present.sort(key=order.index)

    swatches = "".join(
        f'<span class="key-item"><span class="key-swatch" '
        f'style="background:{colours[name]}"></span>{escape(name)}</span>'
        for name in present
    )
    if pooled:
        swatches += (
            f'<span class="key-item"><span class="key-swatch" '
            f'style="background:{OTHER_COLOUR}"></span>Other phyla</span>'
        )

    scales = []
    for ring in rings:
        if ring["key"] == "size":
            top = f"{ring['maximum'] / 1e6:.1f} Mbp"
        else:
            top = f"{ring['maximum']:.0f}%"
        scales.append(
            f'<span class="key-item"><span class="key-ramp" style="background:'
            f'linear-gradient(to right,{ring["ramp"][0]},{ring["ramp"][1]})"></span>'
            f'{escape(ring["label"])} 0–{top}</span>'
        )

    # Named only when a slice is actually blank, so the reader is not sent
    # looking for gaps a complete catalogue does not have.
    gaps = any(value is None for ring in rings for value in ring["values"])
    note = (
        "A blank slice in a ring is a genome that has no value for it."
        if gaps else ""
    )
    return (
        '<div class="tree-key">'
        + (f'<div class="key-row">{swatches}</div>' if swatches else "")
        + (f'<div class="key-row">{"".join(scales)}</div>' if scales else "")
        + (f'<p class="note">{note}</p>' if note else "")
        + "</div>"
    )


def _mix_colours(start, end, fraction):
    """A point on a two-stop ramp, as ``#rrggbb``."""
    fraction = max(0.0, min(1.0, fraction))
    channels = []
    for offset in (1, 3, 5):
        low = int(start[offset:offset + 2], 16)
        high = int(end[offset:offset + 2], 16)
        channels.append(round(low + (high - low) * fraction))
    return "#" + "".join(f"{value:02x}" for value in channels)


def _phylum_colours(connection):
    """One colour per phylum, shared by every figure in the taxonomy section.

    The two phylum figures rank their bars differently — one by how many
    genomes a phylum contributed, the other by how many reads it accounts for
    — so a colour assigned per figure would give the same phylum a different
    colour in each. The assignment is made once here, in abundance order so
    the stacked composition chart still runs from its largest phylum down, and
    both figures look their colours up by name. Phyla past ``TOP_TAXA`` are
    the ones the composition chart pools as Other, and are left out so that
    both figures draw them in ``OTHER_COLOUR``.
    """
    rows = _query(connection, """
        SELECT COALESCE(NULLIF(t.phylum, ''), 'Unclassified') AS name,
               SUM(c.read_count) AS weight
        FROM genome_count AS c
        JOIN genome_taxonomy AS t ON t.genome_id = c.genome_id
        WHERE c.read_count IS NOT NULL
        GROUP BY name
        ORDER BY weight DESC
    """)
    if not rows:
        # No abundances ingested: fall back to how many genomes each phylum
        # holds, which is the only ordering the section can still offer.
        rows = _query(connection, """
            SELECT COALESCE(NULLIF(phylum, ''), 'Unclassified') AS name,
                   COUNT(*) AS weight
            FROM genome_taxonomy
            GROUP BY name
            ORDER BY weight DESC
        """)
    return {
        row["name"]: PALETTE[index % len(PALETTE)]
        for index, row in enumerate(rows[:TOP_TAXA])
    }


def _genome_lineage_blocks(connection):
    """One row per genome, with its GTDB lineage spread across the ranks.

    The rank summary above counts classifications; this is the classification
    itself, which is what a reader who wants to know what a particular bin is
    comes here for. Ranks GTDB-Tk left empty are named rather than blanked, so
    a gap in a lineage is legible as an unclassified rank instead of as a
    missing value.
    """
    # Every rank is a quoted identifier because one of them is "order".
    quoted = [f'"{rank}"' for rank in TAXONOMIC_RANKS]
    ranks = ", ".join(
        f"COALESCE(NULLIF({name}, ''), 'Unclassified') AS {name}"
        for name in quoted
    )
    order = ", ".join(quoted)
    rows = _query(connection, f"""
        SELECT genome_id, {ranks}
        FROM genome_taxonomy
        ORDER BY {order}, genome_id
    """)
    if not rows:
        return []
    return [
        _heading(
            "Lineage per genome",
            "One row per genome in the catalogue — a dereplicated bin — with "
            "its GTDB lineage split across the seven ranks. Rows are ordered by "
            "lineage, so related genomes sit together; click any column heading "
            "to sort by that rank or by genome name instead. Unclassified means "
            "GTDB-Tk found no reference close enough to name the genome at that "
            "rank, and everything finer than it is unclassified too."
        ),
        _table(
            [("Genome", _text)] + [(rank.capitalize(), _text) for rank in TAXONOMIC_RANKS],
            [tuple(row) for row in rows],
        ),
    ]


def _composition_blocks(connection, colours=None):
    """Phylum-level composition per sample, when abundances are also present."""
    import plotly.graph_objects as go

    if colours is None:
        colours = _phylum_colours(connection)

    rows = _query(connection, """
        SELECT COALESCE(NULLIF(t.phylum, ''), 'Unclassified') AS name,
               c.sample_id AS sample_id,
               SUM(c.read_count) AS reads
        FROM genome_count AS c
        JOIN genome_taxonomy AS t ON t.genome_id = c.genome_id
        WHERE c.read_count IS NOT NULL
        GROUP BY name, c.sample_id
    """)
    if not rows:
        return []

    totals = {}
    per_sample = {}
    for row in rows:
        totals[row["name"]] = totals.get(row["name"], 0) + (row["reads"] or 0)
        per_sample[(row["name"], row["sample_id"])] = row["reads"] or 0

    samples = sorted({row["sample_id"] for row in rows})
    sample_totals = {
        sample: sum(per_sample.get((name, sample), 0) for name in totals)
        for sample in samples
    }
    ordered = sorted(totals, key=lambda name: totals[name], reverse=True)
    shown, remainder = ordered[:TOP_TAXA], ordered[TOP_TAXA:]

    figure = go.Figure()
    for name in shown:
        figure.add_bar(
            name=name,
            x=samples,
            y=[_ratio(per_sample.get((name, s), 0), sample_totals.get(s)) for s in samples],
            marker_color=colours.get(name, OTHER_COLOUR),
        )
    if remainder:
        figure.add_bar(
            name="Other",
            x=samples,
            y=[
                _ratio(sum(per_sample.get((name, s), 0) for name in remainder),
                       sample_totals.get(s))
                for s in samples
            ],
            marker_color=OTHER_COLOUR,
        )
    entries = len(shown) + (1 if remainder else 0)
    figure.update_layout(
        barmode="stack",
        xaxis_title="Sample",
        yaxis_title="Reads (% of sample)",
        # One bar per sample and a dozen phyla in the key: the plot takes the
        # full page and the key goes underneath it, laid out as a single row.
        # Giving every entry the same fraction of the width spreads them
        # across it instead of letting Plotly pack them into stacked rows.
        legend=dict(
            orientation="h",
            traceorder="normal",
            title_text="",
            xref="container",
            yref="container",
            yanchor="bottom",
            y=0.01,
            xanchor="left",
            x=0.005,
            entrywidthmode="fraction",
            entrywidth=1 / entries,
            font=dict(size=11),
        ),
        margin=dict(l=45, r=10, t=20, b=110),
        meta=WIDE,
    )
    blocks = [_heading(
        "Composition per sample",
        "The makeup of each sample, as the percentage of its mapped reads "
        "assigned to each phylum. Every bar sums to 100%, so this shows "
        "relative composition only — a phylum can appear to rise in a sample "
        "simply because another one fell."
    )]
    if remainder:
        blocks.append(_note(
            f"The {len(shown)} most abundant phyla are shown individually; the "
            f"remaining {len(remainder)} are pooled as Other."
        ))
    blocks.append(figure)
    return blocks


def _render_function(connection):
    blocks = []
    blocks.extend(_annotation_source_blocks(connection))
    blocks.extend(_annotation_coverage_blocks(connection))
    blocks.extend(_cluster_blocks(connection))
    blocks.extend(_annotation_qc_blocks(connection))
    return blocks


def _annotation_source_blocks(connection):
    """Hits, annotated genes and distinct terms per annotation source."""
    import plotly.graph_objects as go

    total_genes = _scalar(connection, "SELECT COUNT(*) FROM gene", default=0)
    rows = _query(connection, """
        SELECT source,
               COUNT(*) AS hits,
               COUNT(DISTINCT mag || char(31) || gene) AS genes,
               COUNT(DISTINCT annotation_id) AS terms,
               SUM(CASE WHEN is_primary = 1 THEN 1 ELSE 0 END) AS primary_hits
        FROM gene_annotation
        GROUP BY source
        ORDER BY hits DESC
    """)
    if not rows:
        return []

    blocks = []
    if total_genes:
        annotated = _scalar(connection, """
            SELECT COUNT(*) FROM (
                SELECT DISTINCT mag, gene FROM gene_annotation
            )
        """, default=0)
        blocks.append(_paragraph(
            f"{annotated:,} of {total_genes:,} predicted genes "
            f"({_ratio(annotated, total_genes):.1f}%) carry at least one functional hit."
        ))

    blocks.append(_note(
        "Hits counts every annotation record, so one gene matched several times "
        "by the same source is counted several times; Annotated genes counts the "
        "genes behind those hits, and Primary hits the single best match kept per "
        "gene. Genes annotated % is out of all predicted genes in the catalogue."
    ))
    blocks.append(_table(
        [
            ("Source", _text),
            ("Hits", _integer),
            ("Annotated genes", _integer),
            ("Distinct terms", _integer),
            ("Primary hits", _integer),
            ("Genes annotated %", _ONE),
        ],
        [
            tuple(row) + (_ratio(row["genes"], total_genes),)
            for row in rows
        ],
        stats=[
            ("Mean hits per source", 1),
            ("Mean annotated genes", 2),
            ("Mean genes annotated %", 5),
        ],
    ))

    figure = go.Figure()
    figure.add_bar(
        x=[row["source"] for row in rows],
        y=[row["genes"] for row in rows],
        marker_color=PALETTE[0],
    )
    figure.update_layout(xaxis_title="Source", yaxis_title="Annotated genes")
    blocks.append(figure)
    return blocks


def _annotation_coverage_blocks(connection):
    """Per-MAG annotation coverage, as a MAG-by-source heatmap."""
    import plotly.graph_objects as go

    mags = _query(connection, """
        SELECT mag, COUNT(*) AS genes
        FROM gene
        GROUP BY mag
        ORDER BY genes DESC
        LIMIT ?
    """, (TOP_MAGS,))
    if not mags:
        return []
    mag_ids = [row["mag"] for row in mags]
    gene_totals = {row["mag"]: row["genes"] for row in mags}

    cells = _query(connection, f"""
        SELECT mag, source, COUNT(DISTINCT gene) AS annotated
        FROM gene_annotation
        WHERE mag IN ({_placeholders(mag_ids)})
        GROUP BY mag, source
    """, tuple(mag_ids))
    if not cells:
        return []

    sources = sorted({row["source"] for row in cells})
    lookup = {(row["mag"], row["source"]): row["annotated"] for row in cells}
    matrix = [
        [_ratio(lookup.get((mag, source)), gene_totals.get(mag)) for source in sources]
        for mag in mag_ids
    ]

    figure = go.Figure()
    figure.add_heatmap(
        z=matrix,
        x=sources,
        y=mag_ids,
        colorscale="Teal",
        zmin=0,
        colorbar=dict(title="% of genes"),
        hovertemplate="%{y}<br>%{x}<br>%{z:.1f}% of genes<extra></extra>",
    )
    figure.update_layout(
        xaxis_title="Source",
        yaxis_title="MAG",
        height=max(FIGURE_HEIGHT, 18 * len(mag_ids) + 140),
    )

    total_mags = _scalar(connection, "SELECT COUNT(DISTINCT mag) FROM gene", default=0)
    blocks = [_heading(
        "Annotation coverage per MAG",
        "The percentage of each genome's predicted genes that got at least "
        "one hit from each annotation source. Sources differ enormously in "
        "scope by design — a general database annotates most genes, a "
        "specialised one only a handful — so compare genomes down a column "
        "rather than sources across a row."
    )]
    if total_mags > len(mag_ids):
        blocks.append(_note(
            f"Showing the {len(mag_ids)} MAGs with the most predicted genes, "
            f"of {total_mags:,}."
        ))
    blocks.append(figure)
    return blocks


def _cluster_blocks(connection):
    """Gene clusters and regions retained per source and type."""
    rows = _query(connection, """
        SELECT source,
               COALESCE(NULLIF(type, ''), 'unspecified') AS type,
               COUNT(*) AS clusters,
               COUNT(DISTINCT mag) AS mags,
               AVG(gene_count) AS mean_genes
        FROM cluster_annotation
        GROUP BY source, type
        ORDER BY clusters DESC
        LIMIT ?
    """, (TOP_CATEGORIES,))
    if not rows:
        return []
    return [
        _heading(
            "Gene clusters and regions",
            "Features made of several neighbouring genes, such as "
            "biosynthetic gene clusters, grouped by the tool that reported "
            "them and the type it assigned. Clusters counts the features "
            "found and MAGs how many genomes carry at least one of them."
        ),
        _table(
            [
                ("Source", _text),
                ("Type", _text),
                ("Clusters", _integer),
                ("MAGs", _integer),
                ("Mean genes per cluster", _ONE),
            ],
            [tuple(row) for row in rows],
            stats=[("Mean clusters", 2), ("Mean genes per cluster", 4)],
        ),
    ]


def _annotation_qc_blocks(connection):
    """Records reported, retained and rejected by the annotation filters."""
    rows = _query(connection, """
        SELECT source, level,
               SUM(reported_records) AS reported,
               SUM(retained_records) AS retained,
               SUM(rejected_records) AS rejected,
               SUM(unmapped_records) AS unmapped,
               COUNT(DISTINCT mag) AS mags
        FROM annotation_qc
        GROUP BY source, level
        ORDER BY source, level
    """)
    if not rows:
        return []
    return [
        _heading(
            "Annotation filtering",
            "What each annotation tool reported and how much of it Drakkar "
            "kept. Rejected records failed a confidence threshold; unmapped "
            "ones could not be matched back to a predicted gene. A low "
            "Retained % usually means the tool and the filter settings are "
            "mismatched rather than that the underlying data is bad."
        ),
        _table(
            [
                ("Source", _text),
                ("Level", _text),
                ("Reported", _integer),
                ("Retained", _integer),
                ("Rejected", _integer),
                ("Unmapped", _integer),
                ("MAGs", _integer),
                ("Retained %", _ONE),
            ],
            [tuple(row) + (_ratio(row["retained"], row["reported"]),) for row in rows],
            stats=[
                ("Mean reported", 2),
                ("Mean retained", 3),
                ("Mean retained %", 7),
            ],
        ),
    ]


def _render_expression(connection):
    import plotly.graph_objects as go

    rows = _query(connection, """
        SELECT sample_id,
               SUM(count) AS total_counts,
               SUM(CASE WHEN count > 0 THEN 1 ELSE 0 END) AS detected_genes,
               COUNT(*) AS observed_genes
        FROM gene_expression
        GROUP BY sample_id
        ORDER BY sample_id
    """)
    if not rows:
        return []

    genes = _scalar(connection, "SELECT COUNT(*) FROM expressed_gene", default=0)
    blocks = [_paragraph(
        f"Quantification covers {_quantity(genes, 'gene')} across "
        f"{_quantity(len(rows), 'sample')}."
    )]
    blocks.append(_table(
        [
            ("Sample", _text),
            ("Total counts", _integer),
            ("Genes detected", _integer),
            ("Genes quantified", _integer),
            ("Genes detected %", _ONE),
        ],
        [
            tuple(row) + (_ratio(row["detected_genes"], row["observed_genes"]),)
            for row in rows
        ],
        stats=[
            ("Mean total counts", 1),
            ("Mean genes detected", 2),
            ("Mean genes detected %", 4),
        ],
    ))

    samples = [row["sample_id"] for row in rows]
    figure = go.Figure()
    figure.add_bar(
        name="Total counts",
        x=samples,
        y=[row["total_counts"] for row in rows],
        marker_color=PALETTE[0],
    )
    figure.add_scatter(
        name="Genes detected",
        x=samples,
        y=[row["detected_genes"] for row in rows],
        mode="markers",
        marker=dict(color=PALETTE[1], size=9),
        yaxis="y2",
    )
    figure.update_layout(
        xaxis_title="Sample",
        yaxis_title="Assigned counts",
        yaxis2=dict(title="Genes detected", overlaying="y", side="right",
                    showgrid=False, rangemode="tozero"),
        legend_title_text="",
    )
    blocks.append(_heading(
        "Expression per sample",
        "Assigned counts per sample (bars, left axis) against the number of "
        "genes with at least one read (points, right axis). The two usually "
        "rise together; a sample with far fewer detected genes than the rest "
        "is more often shallowly sequenced than biologically different."
    ))
    blocks.append(figure)

    lengths = _query(connection, """
        SELECT COUNT(*) AS genes, AVG(length) AS mean_length,
               MIN(length) AS min_length, MAX(length) AS max_length
        FROM expressed_gene
        WHERE length IS NOT NULL
    """)
    if lengths and lengths[0]["genes"]:
        blocks.append(_heading(
            "Quantified genes",
            "The size of the gene set that counts were assigned to. Gene "
            "length matters when comparing raw counts, because a longer gene "
            "collects more reads than a shorter one expressed just as "
            "strongly."
        ))
        blocks.append(_highlights([
            ("Genes with a length", _integer(lengths[0]["genes"])),
            ("Mean length", _integer(lengths[0]["mean_length"])),
            ("Shortest", _integer(lengths[0]["min_length"])),
            ("Longest", _integer(lengths[0]["max_length"])),
        ]))
    return blocks


# ---------------------------------------------------------------------------
# Antimicrobial resistance
# ---------------------------------------------------------------------------

# How the reconciler labels a locus by the callers that backed it, and the
# order the labels are stacked in: the strongest support first.
AMR_SUPPORT_LABELS = (
    ("amrfinder_and_rgi", "AMRFinderPlus and RGI"),
    ("amrfinder_only", "AMRFinderPlus only"),
    ("rgi_only", "RGI only"),
    ("other", "Other"),
)

# How two callers that both hit a locus agreed, weakest agreement last.
AMR_CONCORDANCE_LABELS = {
    "gene_match": "Same gene",
    "drug_class_match": "Same drug class only",
    "unresolved": "No shared gene or drug class",
    "single_source": "Only one caller",
}

# geNomad's mobile-element calls, plus the label for a locus none of them hit.
AMR_CONTEXT_LABELS = {
    "plasmid": "Plasmid",
    "virus": "Virus",
    "provirus": "Provirus",
}
AMR_CHROMOSOMAL = "No mobile-element call"


def _render_amr(connection):
    assemblies = _query(connection, """
        SELECT assembly_id, assembly_type, organism, contig_count, total_length,
               amrfinder_hits, rgi_hits, mobility_regions, amr_loci,
               multi_tool_loci, mobile_loci,
               amrfinder_hits_without_coordinates,
               rgi_hits_without_coordinates
        FROM amr_assembly
        ORDER BY assembly_id
    """)
    if not assemblies:
        return []

    blocks = _amr_overview_blocks(connection, assemblies)
    blocks.extend(_amr_support_blocks(connection))
    blocks.extend(_amr_drug_class_blocks(connection))
    blocks.extend(_amr_gene_blocks(connection))
    blocks.extend(_amr_context_blocks(connection))
    blocks.extend(_amr_qc_blocks(assemblies))
    return blocks


def _amr_overview_blocks(connection, assemblies):
    """What was screened, and how much resistance came out of it."""
    loci = _scalar(connection, "SELECT COUNT(*) FROM amr_locus", default=0) or sum(
        row["amr_loci"] or 0 for row in assemblies
    )
    both_tools = _scalar(
        connection,
        "SELECT COUNT(*) FROM amr_locus WHERE support_status = 'amrfinder_and_rgi'",
        default=0,
    ) or sum(row["multi_tool_loci"] or 0 for row in assemblies)
    mobile = _scalar(connection, """
        SELECT COUNT(*) FROM (
            SELECT DISTINCT assembly_id, locus_id FROM amr_mobility
        )
    """, default=0) or sum(row["mobile_loci"] or 0 for row in assemblies)
    genes = _scalar(connection, """
        SELECT COUNT(DISTINCT primary_gene) FROM amr_locus
        WHERE primary_gene IS NOT NULL AND primary_gene <> ''
    """, default=0)
    drug_classes = _scalar(
        connection, "SELECT COUNT(DISTINCT drug_class) FROM amr_drug_class", default=0
    )

    blocks = [_paragraph(
        f"AMRFinderPlus and RGI screened {_quantity(len(assemblies), 'assembly', 'assemblies')}, "
        f"and their calls were reconciled into {_quantity(loci, 'resistance locus', 'resistance loci')}."
    )]
    blocks.append(_highlights([
        ("Assemblies screened", _integer(len(assemblies))),
        ("Resistance loci", _integer(loci)),
        ("Called by both tools", _integer(both_tools),
         _cards_percent(_ratio(both_tools, loci))),
        ("On a mobile element", _integer(mobile),
         _cards_percent(_ratio(mobile, loci))),
        ("Distinct genes", _integer(genes) if genes else None),
        ("Drug classes", _integer(drug_classes) if drug_classes else None),
    ]))

    blocks.append(_heading(
        "Resistance per assembly",
        "AMRFinderPlus hits and RGI hits are what each caller reported on its "
        "own, so an occurrence both of them found is counted twice between "
        "them; Loci is the reconciled count, one per gene occurrence. Both "
        "tools % is the share of loci that survived as agreement between the "
        "two callers rather than as a single-caller call, and Mobile % the "
        "share sitting inside a plasmid or a (pro)virus geNomad called."
    ))
    blocks.append(_table(
        [
            ("Assembly", _text),
            ("Type", _text),
            ("Organism", _text),
            ("Contigs", _integer),
            ("Length", _integer),
            ("AMRFinderPlus hits", _integer),
            ("RGI hits", _integer),
            ("Mobile regions", _integer),
            ("Loci", _integer),
            ("Both tools", _integer),
            ("Mobile loci", _integer),
            ("Both tools %", _ONE),
            ("Mobile %", _ONE),
        ],
        [
            (
                row["assembly_id"], row["assembly_type"], row["organism"],
                row["contig_count"], row["total_length"], row["amrfinder_hits"],
                row["rgi_hits"], row["mobility_regions"], row["amr_loci"],
                row["multi_tool_loci"], row["mobile_loci"],
                _ratio(row["multi_tool_loci"], row["amr_loci"]),
                _ratio(row["mobile_loci"], row["amr_loci"]),
            )
            for row in assemblies
        ],
        stats=[
            ("Mean loci per assembly", 8),
            ("Mean both tools %", 11),
            ("Mean mobile %", 12),
        ],
    ))
    return blocks


def _amr_support_blocks(connection):
    """Which callers backed each locus, and how well they agreed."""
    import plotly.graph_objects as go

    rows = _query(connection, """
        SELECT assembly_id, support_status, COUNT(*) AS loci
        FROM amr_locus
        GROUP BY assembly_id, support_status
    """)
    if not rows:
        return []

    totals = {}
    for row in rows:
        totals.setdefault(row["assembly_id"], {})[row["support_status"]] = row["loci"]
    assemblies = sorted(totals, key=lambda name: -sum(totals[name].values()))[:TOP_GENOMES]

    blocks = [_heading(
        "Caller support",
        "AMRFinderPlus and RGI use different reference databases and different "
        "acceptance thresholds, so neither is a subset of the other. A locus "
        "both of them called is the most defensible; a single-caller locus is "
        "not wrong, but rests on one database alone."
    )]
    figure = go.Figure()
    for index, (status, label) in enumerate(AMR_SUPPORT_LABELS):
        values = [totals[name].get(status, 0) for name in assemblies]
        if not any(values):
            continue
        figure.add_bar(
            name=label, x=assemblies, y=values, marker_color=PALETTE[index]
        )
    figure.update_layout(
        barmode="stack",
        xaxis_title="Assembly",
        yaxis_title="Resistance loci",
        legend_title_text="",
    )
    blocks.append(figure)
    if len(totals) > len(assemblies):
        blocks.append(_note(
            f"The {len(assemblies)} assemblies with the most loci are shown; "
            f"{len(totals) - len(assemblies)} more are in the table above."
        ))

    concordance = _query(connection, """
        SELECT concordance, COUNT(*) AS loci
        FROM amr_locus
        WHERE support_status = 'amrfinder_and_rgi'
        GROUP BY concordance
        ORDER BY loci DESC
    """)
    if concordance:
        total = sum(row["loci"] for row in concordance)
        blocks.append(_heading(
            "Agreement on the loci both callers found",
            "Overlapping calls are merged on coordinates, so two callers can "
            "land on the same locus and still name different genes. Same gene "
            "is full agreement; same drug class only means they disagree on "
            "the gene but agree on what it confers resistance to; no shared "
            "gene or drug class marks a locus worth inspecting by hand."
        ))
        blocks.append(_table(
            [("Agreement", _text), ("Loci", _integer), ("Share %", _ONE)],
            [
                (
                    AMR_CONCORDANCE_LABELS.get(row["concordance"], row["concordance"]),
                    row["loci"],
                    _ratio(row["loci"], total),
                )
                for row in concordance
            ],
        ))
    return blocks


def _amr_drug_class_blocks(connection):
    """Which drug classes the catalogue of loci confers resistance to."""
    import plotly.graph_objects as go

    rows = _query(connection, """
        SELECT drug_class,
               COUNT(DISTINCT assembly_id || char(31) || COALESCE(locus_id, hit_id)) AS loci,
               COUNT(DISTINCT assembly_id) AS assemblies
        FROM amr_drug_class
        WHERE drug_class IS NOT NULL AND drug_class <> ''
        GROUP BY drug_class
        ORDER BY loci DESC
    """)
    if not rows:
        return []

    blocks = [_heading(
        "Drug classes",
        "A locus that confers resistance to several drug classes is counted "
        "once under each of them, so these counts add up to more than the "
        "number of loci. Read them as how often each class is covered, not as "
        "a partition of the loci."
    )]
    shown = rows[:TOP_CATEGORIES]
    figure = go.Figure()
    figure.add_bar(
        x=[row["loci"] for row in reversed(shown)],
        y=[row["drug_class"] for row in reversed(shown)],
        orientation="h",
        marker_color=PALETTE[0],
    )
    figure.update_layout(
        xaxis_title="Resistance loci",
        yaxis_title="",
        height=max(FIGURE_HEIGHT, 24 * len(shown) + 120),
        margin=dict(l=260),
    )
    blocks.append(figure)
    if len(rows) > len(shown):
        blocks.append(_note(
            f"The {len(shown)} most frequent of {len(rows)} drug classes are "
            "drawn; the table below carries all of them."
        ))
    blocks.append(_table(
        [("Drug class", _text), ("Loci", _integer), ("Assemblies", _integer)],
        [(row["drug_class"], row["loci"], row["assemblies"]) for row in rows],
    ))
    return blocks


def _amr_gene_blocks(connection):
    """The genes behind the loci, ranked by how often they were recovered."""
    rows = _query(connection, """
        SELECT l.primary_gene AS gene,
               COUNT(*) AS loci,
               COUNT(DISTINCT l.assembly_id) AS assemblies,
               SUM(CASE WHEN l.support_status = 'amrfinder_and_rgi' THEN 1 ELSE 0 END)
                   AS both_tools,
               SUM(CASE WHEN m.locus_id IS NOT NULL THEN 1 ELSE 0 END) AS mobile_loci,
               MAX(l.drug_classes) AS drug_classes
        FROM amr_locus AS l
        LEFT JOIN (
            SELECT DISTINCT assembly_id, locus_id FROM amr_mobility
        ) AS m ON m.assembly_id = l.assembly_id AND m.locus_id = l.locus_id
        WHERE l.primary_gene IS NOT NULL AND l.primary_gene <> ''
        GROUP BY l.primary_gene
        ORDER BY loci DESC, gene
    """)
    if not rows:
        return []

    return [
        _heading(
            "Resistance genes",
            "One row per gene, counting the loci it was recovered from. A gene "
            "in many assemblies is either widespread in the sample set or "
            "sitting on something that moves between them, which is what the "
            "mobile count next to it is for. Drug classes lists every class "
            "the gene was annotated with across its loci."
        ),
        _table(
            [
                ("Gene", _text),
                ("Loci", _integer),
                ("Assemblies", _integer),
                ("Both tools", _integer),
                ("Mobile loci", _integer),
                ("Drug classes", _text),
            ],
            [
                (row["gene"], row["loci"], row["assemblies"], row["both_tools"],
                 row["mobile_loci"], row["drug_classes"])
                for row in rows
            ],
        ),
    ]


def _amr_context_blocks(connection):
    """Where the loci sit: on a plasmid, in a (pro)virus, or on neither."""
    import plotly.graph_objects as go

    regions = _query(connection, """
        SELECT context_type,
               COUNT(*) AS regions,
               AVG(length) AS mean_length,
               AVG(score) AS mean_score
        FROM amr_mobility_region
        GROUP BY context_type
        ORDER BY regions DESC
    """)
    contexts = _query(connection, """
        SELECT context_type,
               COUNT(DISTINCT assembly_id || char(31) || locus_id) AS loci
        FROM amr_mobility
        GROUP BY context_type
        ORDER BY loci DESC
    """)
    if not regions and not contexts:
        return []

    total_loci = _scalar(connection, "SELECT COUNT(*) FROM amr_locus", default=0)
    mobile = _scalar(connection, """
        SELECT COUNT(*) FROM (
            SELECT DISTINCT assembly_id, locus_id FROM amr_mobility
        )
    """, default=0)

    blocks = [_heading(
        "Mobile-element context",
        "geNomad calls plasmids and viruses on the same contigs, and a locus "
        "is called mobile when it overlaps one of those calls. A locus can "
        "overlap both a plasmid and a provirus, so the per-context counts can "
        "exceed the number of mobile loci. This is positional evidence, not "
        "proof of transfer: it says the gene sits inside a region that looks "
        "mobile, not that it has moved."
    )]

    if total_loci:
        labels = [
            AMR_CONTEXT_LABELS.get(row["context_type"], row["context_type"])
            for row in contexts
        ]
        values = [row["loci"] for row in contexts]
        colours = [PALETTE[index % len(PALETTE)] for index in range(len(contexts))]
        labels.append(AMR_CHROMOSOMAL)
        values.append(max(total_loci - mobile, 0))
        colours.append(OTHER_COLOUR)
        figure = go.Figure()
        figure.add_bar(x=labels, y=values, marker_color=colours)
        figure.update_layout(xaxis_title="", yaxis_title="Resistance loci")
        blocks.append(figure)

    if contexts:
        blocks.append(_table(
            [("Context", _text), ("Resistance loci", _integer), ("Share of loci %", _ONE)],
            [
                (
                    AMR_CONTEXT_LABELS.get(row["context_type"], row["context_type"]),
                    row["loci"],
                    _ratio(row["loci"], total_loci),
                )
                for row in contexts
            ],
        ))

    if regions:
        blocks.append(_heading(
            "Mobile elements found",
            "Every plasmid and (pro)virus geNomad called, whether or not a "
            "resistance gene landed in it. The score is geNomad's own "
            "confidence in the call; a low-scoring region carrying a locus is "
            "the one to check before calling that locus mobile."
        ))
        blocks.append(_table(
            [
                ("Context", _text),
                ("Regions", _integer),
                ("Mean length", _integer),
                ("Mean score", _TWO),
            ],
            [
                (
                    AMR_CONTEXT_LABELS.get(row["context_type"], row["context_type"]),
                    row["regions"], row["mean_length"], row["mean_score"],
                )
                for row in regions
            ],
        ))
    return blocks


def _amr_qc_blocks(assemblies):
    """Hits the callers reported without usable coordinates.

    A hit with no contig and no start/end cannot be reconciled against the
    other caller or placed inside a mobile element, so it is counted here
    rather than silently dropped from the locus figures above.
    """
    rows = [
        (
            row["assembly_id"],
            row["amrfinder_hits_without_coordinates"],
            row["rgi_hits_without_coordinates"],
        )
        for row in assemblies
        if (row["amrfinder_hits_without_coordinates"] or 0)
        or (row["rgi_hits_without_coordinates"] or 0)
    ]
    if not rows:
        return []
    return [
        _heading(
            "Hits without coordinates",
            "These hits were reported by a caller but carry no position on a "
            "contig, so they could not be merged into a locus or tested "
            "against a mobile element. They are excluded from every count "
            "above; a large number here means the locus figures understate "
            "what the callers found."
        ),
        _table(
            [
                ("Assembly", _text),
                ("AMRFinderPlus", _integer),
                ("RGI", _integer),
            ],
            rows,
        ),
    ]


def _render_resources(connection):
    blocks = _run_blocks(connection)
    blocks.extend(_failure_blocks(connection))
    blocks.extend(_benchmark_blocks(connection))
    return blocks


def _run_blocks(connection):
    rows = _query(connection, """
        SELECT run_id, drakkar_version, command, modules, started_at,
               finished_at, status
        FROM run
        ORDER BY run_id
    """)
    if not rows:
        return []

    blocks = [_paragraph(
        f"This output directory records {_quantity(len(rows), 'Drakkar run')}."
    )]
    blocks.append(_table(
        [
            ("Run", _text),
            ("Version", _text),
            ("Command", _text),
            ("Modules", _text),
            ("Started", _timestamp),
            ("Finished", _timestamp),
            ("Duration", _text),
            ("Status", _text),
        ],
        [
            (
                row["run_id"], row["drakkar_version"], row["command"], row["modules"],
                row["started_at"], row["finished_at"],
                _duration(row["started_at"], row["finished_at"]), row["status"],
            )
            for row in rows
        ],
    ))
    return blocks


# What the two statuses drakkar.failures records read as on the page. A job
# Snakemake resubmitted until it succeeded left the run's results complete, so
# it is named for what happened rather than for the failure it started as.
FAILURE_OUTCOMES = {
    "recovered": "Recovered on retry",
    "failed": "Never succeeded",
}


def _failure_outcome(status):
    return FAILURE_OUTCOMES.get(str(status or "").lower(), status or "")


def _failure_blocks(connection):
    """Every rule that lost a job, grouped by the rule and the failing wildcards.

    The rows come from the failure table ``drakkar.failures`` writes after a
    run, which already groups the log's individual error events by rule and
    wildcards, so a job that failed twice and then succeeded is one row here
    carrying both facts: how many attempts it cost, and that a retry recovered
    it. A run in which nothing failed writes no such table, so an absent
    subsection means a clean run rather than a missing file.
    """
    rows = _query(connection, """
        SELECT run_id, rule, target, attempts, status, category, reason,
               slurm_state, detail, action, job_log
        FROM run_failure
        ORDER BY CASE WHEN LOWER(COALESCE(status, '')) = 'recovered'
                      THEN 1 ELSE 0 END,
                 COALESCE(attempts, 0) DESC, rule, target, run_id
    """)
    if not rows:
        return []

    unresolved = [row for row in rows if str(row["status"] or "").lower() != "recovered"]
    recovered = [row for row in rows if str(row["status"] or "").lower() == "recovered"]

    blocks = [
        _heading(
            "Errors",
            "One row per rule and per set of wildcards that failed, not per "
            "failed attempt: Snakemake resubmits a failed job, so a rule that "
            "failed on one sample and succeeded on the retry appears once, "
            "with the number of attempts it cost. The target is the wildcards "
            "the failing job carried — usually the sample or assembly it was "
            "running on. Attempts counts the failed ones: a recovered row "
            "succeeded on the attempt after those and left the run's results "
            "complete, while a row that never succeeded is one whose outputs, "
            "and everything downstream of them, are missing."
        ),
        _highlights([
            ("Errors recorded", _integer(len(rows))),
            ("Recovered on retry", _integer(len(recovered))),
            ("Never succeeded", _integer(len(unresolved))),
            ("Rules affected", _integer(len({row["rule"] for row in rows}))),
            ("Failed attempts", _integer(sum(row["attempts"] or 0 for row in rows))),
        ]),
        _table(
            [
                ("Run", _text),
                ("Rule", _text),
                ("Target", _text),
                ("Attempts", _integer),
                ("Outcome", _text),
                ("Problem", _text),
                ("Scheduler state", _text),
                ("What went wrong", _text),
                ("Job log", _text),
            ],
            [
                (
                    row["run_id"], row["rule"], row["target"], row["attempts"],
                    _failure_outcome(row["status"]), row["category"],
                    row["slurm_state"], row["detail"], row["job_log"],
                )
                for row in rows
            ],
        ),
    ]
    blocks.extend(_failure_action_blocks(unresolved, len(recovered)))
    return blocks


def _failure_action_blocks(unresolved, recovered):
    """What to do about the failures that a retry did not settle.

    The advice is one line per kind of problem rather than one per row: every
    out-of-memory job in a run is answered by the same larger multiplier, and
    repeating that sentence down a column of a hundred rows hides it.
    """
    if not unresolved:
        if not recovered:
            return []
        return [_note(
            f"Every failure above was recovered by a retry: "
            f"{_quantity(recovered, 'job')} cost more attempts than it should "
            "have, but the run's results are complete and nothing needs "
            "rerunning."
        )]

    grouped = {}
    for row in unresolved:
        key = row["category"] or "unknown"
        entry = grouped.setdefault(
            key, {"jobs": 0, "reason": row["reason"], "action": row["action"]}
        )
        entry["jobs"] += 1
        # An older run may have stored one of the two texts and not the other.
        entry["reason"] = entry["reason"] or row["reason"]
        entry["action"] = entry["action"] or row["action"]

    return [
        _heading(
            "What to do next",
            "One line per kind of problem still outstanding, with the "
            "recommendation Drakkar made when the run finished. Rerunning the "
            "same command resumes where the workflow stopped: what already "
            "completed is not recomputed."
        ),
        _table(
            [
                ("Problem", _text),
                ("Jobs", _integer),
                ("What it means", _text),
                ("Suggested action", _text),
            ],
            [
                (name, entry["jobs"], entry["reason"], entry["action"])
                for name, entry in sorted(
                    grouped.items(), key=lambda item: (-item[1]["jobs"], item[0])
                )
            ],
        ),
    ]


# Why a run carries no usage figures. Keyed by the status drakkar.benchmark
# stamps into drakkar_<run_id>.resources.yaml.
BENCHMARK_STATUS_NOTES = {
    "skipped": "Resource benchmarking was skipped for this run (--skip-benchmark).",
    "unsupported_profile": (
        "Resource benchmarking currently covers only runs launched with the "
        "slurm profile, so no usage figures were collected."
    ),
    "log_missing": (
        "The Snakemake log was not found, so no resource benchmark could be built."
    ),
    "no_submitted_jobs": (
        "No submitted cluster jobs were detected in the Snakemake log."
    ),
    "accounting_unavailable": (
        "SLURM accounting could not be queried, so requested resources are "
        "recorded without the matching actual usage."
    ),
}


def _percent(count, total):
    """A count as a percentage of a total, or None when the total is zero."""
    return _ratio(count, total) if total else None


def _minutes(seconds):
    return None if seconds is None else float(seconds) / 60.0


def _efficiency(fraction):
    """Benchmark efficiencies are stored as fractions; the page shows percent."""
    return None if fraction is None else 100.0 * float(fraction)


def _benchmark_blocks(connection):
    """Requested versus used resources, and how many jobs actually succeeded."""
    blocks = []
    blocks.extend(_job_outcome_blocks(connection))
    blocks.extend(_rule_resource_blocks(connection))
    blocks.extend(_job_resource_blocks(connection))
    if not blocks:
        blocks.extend(_benchmark_status_blocks(connection))
    return blocks


def _benchmark_status_blocks(connection):
    """Explain the absence of usage figures, when a run recorded a reason."""
    rows = _query(connection, """
        SELECT run_id, status FROM run_benchmark ORDER BY run_id
    """)
    notes = []
    for row in rows:
        message = BENCHMARK_STATUS_NOTES.get(row["status"])
        if message:
            notes.append(_note(f"{row['run_id']}: {message}"))
    if not notes:
        return []
    return [
        _heading(
            "Resource benchmark",
            "No compute usage figures were recorded for this run. The reason "
            "reported by each run follows."
        ),
        *notes,
    ]


def _job_outcome_blocks(connection):
    """How many submitted jobs completed, failed, or never reported back."""
    totals = _query(connection, """
        SELECT COUNT(*) AS launches,
               COUNT(DISTINCT logical_job_key) AS logical_jobs,
               SUM(CASE WHEN UPPER(COALESCE(state, '')) = 'COMPLETED'
                        THEN 1 ELSE 0 END) AS completed,
               SUM(CASE WHEN state IS NULL OR TRIM(state) = ''
                        THEN 1 ELSE 0 END) AS unknown,
               SUM(CASE WHEN COALESCE(oom, 0) = 1 THEN 1 ELSE 0 END) AS oom,
               SUM(CASE WHEN COALESCE(timeout, 0) = 1 THEN 1 ELSE 0 END) AS timeout,
               SUM(CASE WHEN COALESCE(attempt, 1) > 1 THEN 1 ELSE 0 END) AS retries,
               SUM(COALESCE(alloc_cpus, 0) * COALESCE(elapsed_sec, 0)) AS allocated,
               SUM(COALESCE(cpu_time_sec, 0)) AS used,
               MAX(max_rss_mb) AS peak_rss,
               SUM(COALESCE(elapsed_sec, 0)) AS elapsed
        FROM benchmark_job
    """)
    if not totals or not totals[0]["launches"]:
        return []
    total = totals[0]
    launches = total["launches"]
    # A launch with no accounting row is neither a success nor a failure: it is
    # a job sacct could not speak for, and is counted on its own.
    failed = launches - (total["completed"] or 0) - (total["unknown"] or 0)

    blocks = [
        _heading(
            "Job outcomes",
            "Whether the individual cluster jobs behind this run succeeded. "
            "Failures here are workflow-level events — a job that ran out of "
            "memory or hit its time limit — and Snakemake retries them, so a "
            "run can finish successfully with failed jobs on record."
        ),
        _paragraph(
            f"{_quantity(launches, 'job')} were submitted to the scheduler for "
            f"{_quantity(total['logical_jobs'] or 0, 'distinct workflow job')}; "
            "a job submitted more than once appears once per attempt."
        ),
        _highlights([
            ("Jobs submitted", _integer(launches)),
            ("Successful", f"{_integer(total['completed'])} "
                           f"({_ONE(_percent(total['completed'], launches))}%)"),
            ("Failed", f"{_integer(failed)} "
                       f"({_ONE(_percent(failed, launches))}%)"),
            ("Relaunches", _integer(total["retries"])),
            ("CPU efficiency", f"{_ONE(_percent(total['used'], total['allocated']))}%"
                               if total["allocated"] else None),
            ("Peak memory", f"{_integer(total['peak_rss'])} MB"
                            if total["peak_rss"] else None),
            ("Job wall time", f"{_TWO((total['elapsed'] or 0) / 3600.0)} h"
                              if total["elapsed"] else None),
        ]),
    ]

    states = _query(connection, """
        SELECT COALESCE(NULLIF(TRIM(state), ''), 'No accounting record') AS state,
               COUNT(*) AS launches,
               COUNT(DISTINCT rule) AS rules
        FROM benchmark_job
        GROUP BY 1
        ORDER BY launches DESC
    """)
    figure = _job_state_figure(states, launches)
    if figure is not None:
        blocks.append(figure)
    blocks.append(_table(
        [
            ("Final state", _text),
            ("Jobs", _integer),
            ("% of jobs", _ONE),
            ("Rules affected", _integer),
        ],
        [
            (row["state"], row["launches"], _percent(row["launches"], launches),
             row["rules"])
            for row in states
        ],
    ))

    failure_notes = []
    if total["oom"]:
        failure_notes.append(
            f"{_quantity(total['oom'], 'job')} ran out of memory"
        )
    if total["timeout"]:
        failure_notes.append(
            f"{_quantity(total['timeout'], 'job')} hit the requested time limit"
        )
    if total["unknown"]:
        failure_notes.append(
            f"{_quantity(total['unknown'], 'job')} had no accounting record"
        )
    if failure_notes:
        blocks.append(_note("Of those: " + ", ".join(failure_notes) + "."))
    return blocks


# The colour a final job state is drawn in, matched on the leading word so
# that sacct's qualified states ("CANCELLED by 501") land with their family.
# Success is the one green in the palette; the failure modes a run is tuned
# against — running out of memory, running out of time — are kept apart from
# each other so a bar says at a glance which one to fix.
JOB_STATE_COLOURS = (
    ("COMPLETED", PALETTE[2]),
    ("OUT_OF_MEMORY", PALETTE[3]),
    ("TIMEOUT", PALETTE[1]),
    ("FAILED", PALETTE[8]),
    ("CANCELLED", PALETTE[6]),
    ("NODE_FAIL", PALETTE[4]),
    ("PREEMPTED", PALETTE[5]),
)


def _job_state_colour(state):
    upper = str(state or "").upper()
    for name, colour in JOB_STATE_COLOURS:
        if upper.startswith(name):
            return colour
    return OTHER_COLOUR


def _job_state_figure(states, launches):
    """The mix of final job states as one stacked bar across the page.

    The table below carries the same numbers, but a run with a handful of
    failures among hundreds of jobs reads as a table of near-identical
    percentages; as a single bar, the failing slice is visible without being
    looked for.
    """
    import plotly.graph_objects as go

    if not states or not launches:
        return None

    figure = go.Figure()
    for row in states:
        share = _percent(row["launches"], launches)
        figure.add_bar(
            name=row["state"],
            y=["Jobs"],
            x=[row["launches"]],
            orientation="h",
            marker_color=_job_state_colour(row["state"]),
            # Only a slice with room for it is labelled; the rest are in the
            # key and the hover.
            text=[f"{row['launches']:,}"],
            textposition="inside",
            insidetextanchor="middle",
            cliponaxis=False,
            hovertemplate=(
                f"{row['state']}<br>%{{x:,}} jobs"
                + (f" ({share:.1f}% of all jobs)" if share is not None else "")
                + "<extra></extra>"
            ),
        )
    figure.update_layout(
        barmode="stack",
        xaxis_title="Jobs",
        # One bar and no categories to name on it, so the whole vertical axis
        # goes; the key sits above the bar, clear of the axis title below it.
        yaxis=dict(visible=False),
        legend=dict(
            orientation="h",
            traceorder="normal",
            title_text="",
            yanchor="bottom",
            y=1.06,
            xanchor="left",
            x=0,
            entrywidthmode="fraction",
            entrywidth=1 / len(states),
        ),
        margin=dict(l=30, t=60, b=60),
        height=240,
        uniformtext=dict(mode="hide", minsize=11),
        meta=WIDE,
    )
    return figure


# The workflow sizes most runtime requests from the input size but never lets
# the result fall below a floor: ``cap_runtime(max(15, int(input.size_mb / 10)
# ...))`` asks for a quarter of an hour however small the input is. Jobs that
# sit on that floor are unaffected by the size coefficient, so a rule that
# looks over-provisioned is only tunable through the floor, and the table
# names the floor beside the request rather than leaving it to be looked up in
# the rule files.
WORKFLOW_DIR = Path(__file__).resolve().parents[1] / "workflow"
_RULE_HEADER = re.compile(r"^(?:rule|checkpoint)\s+(\w+)\s*:")
_RUNTIME_FLOOR = re.compile(r"max\(\s*(\d+)\s*,")
_RUNTIME_CONSTANT = re.compile(r"cap_runtime\(\s*(\d+)\s*[)*]")


def _runtime_floors():
    """The smallest runtime, in minutes, each rule of the workflow can request.

    Read from the rule files of the *installed* Drakkar rather than from the
    database, which stores what the scheduler was asked for and not the
    expression that produced it. A rule name defined in more than one module
    file with more than one floor is left out: only one of those modules ran,
    and the rule files alone do not say which.
    """
    floors = {}
    paths = [WORKFLOW_DIR / "Snakefile", *sorted((WORKFLOW_DIR / "rules").glob("*.smk"))]
    for path in paths:
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        rule = None
        for line in text.splitlines():
            header = _RULE_HEADER.match(line)
            if header:
                rule = header.group(1)
                continue
            if rule is None or "runtime=" not in line:
                continue
            expression = line.split("runtime=", 1)[1]
            match = (
                _RUNTIME_FLOOR.search(expression)
                or _RUNTIME_CONSTANT.search(expression)
            )
            if match:
                floors.setdefault(rule, set()).add(int(match.group(1)))
            rule = None
    return {
        rule: values.pop()
        for rule, values in floors.items()
        if len(values) == 1
    }


def _rule_job_statistics(connection):
    """Per-rule figures the rule roll-up does not carry: totals and ceilings.

    ``benchmark_rule`` stores medians, which is what a rule's typical job
    needs; a resource request has to cover its heaviest jobs instead. The 95th
    percentile of what the rule's jobs actually used is that ceiling, and it,
    like the rule's total runtime, can only be computed from the per-job rows.
    """
    rows = _query(connection, """
        SELECT rule, elapsed_sec, memory_efficiency, runtime_efficiency,
               requested_runtime_min
        FROM benchmark_job
    """)
    grouped = {}
    for row in rows:
        entry = grouped.setdefault(
            row["rule"],
            {"elapsed": [], "memory": [], "runtime": [], "requested": []},
        )
        entry["elapsed"].append(row["elapsed_sec"])
        entry["memory"].append(row["memory_efficiency"])
        entry["runtime"].append(row["runtime_efficiency"])
        if row["requested_runtime_min"] is not None:
            entry["requested"].append(float(row["requested_runtime_min"]))
    return {
        rule: {
            "jobs": len(entry["elapsed"]),
            "total_elapsed_sec": sum(value or 0 for value in entry["elapsed"]),
            "memory_p95": _percentile(entry["memory"], 0.95),
            "runtime_p95": _percentile(entry["runtime"], 0.95),
            **_smallest_request(entry["requested"]),
        }
        for rule, entry in grouped.items()
    }


def _smallest_request(requested):
    """The rule's smallest runtime request and the share of jobs asking for it.

    The share is what says whether the floor matters: a rule whose jobs all
    ask for the same smallest amount is pinned there, while one where a single
    job happens to be the smallest is sized by its inputs.
    """
    if not requested:
        return {"min_requested_runtime_min": None, "jobs_at_min_request": None}
    smallest = min(requested)
    at_smallest = sum(1 for value in requested if value == smallest)
    return {
        "min_requested_runtime_min": smallest,
        "jobs_at_min_request": 100.0 * at_smallest / len(requested),
    }


def _rule_resource_blocks(connection):
    """Median requested versus actually used resources, per Snakemake rule."""
    rows = _query(connection, """
        SELECT rule,
               SUM(COALESCE(launches, 0)) AS launches,
               SUM(COALESCE(failed_launches, 0)) AS failed,
               AVG(median_requested_cpus) AS requested_cpus,
               AVG(median_alloc_cpus) AS alloc_cpus,
               AVG(median_requested_mem_mb) AS requested_mem_mb,
               AVG(median_max_rss_mb) AS max_rss_mb,
               AVG(median_memory_efficiency) AS memory_efficiency,
               AVG(median_requested_runtime_min) AS requested_runtime_min,
               AVG(median_elapsed_sec) AS elapsed_sec,
               AVG(median_runtime_efficiency) AS runtime_efficiency,
               SUM(COALESCE(allocated_cpu_sec, 0)) AS allocated_cpu_sec,
               SUM(COALESCE(used_cpu_sec, 0)) AS used_cpu_sec
        FROM benchmark_rule
        GROUP BY rule
        ORDER BY allocated_cpu_sec DESC, rule
    """)
    if not rows:
        return []

    per_job = _rule_job_statistics(connection)
    floors = _runtime_floors()

    table_rows = [
        (
            row["rule"],
            row["launches"],
            row["failed"],
            row["requested_cpus"],
            row["alloc_cpus"],
            row["requested_mem_mb"],
            row["max_rss_mb"],
            _efficiency(row["memory_efficiency"]),
            _efficiency(per_job.get(row["rule"], {}).get("memory_p95")),
            row["requested_runtime_min"],
            floors.get(row["rule"]),
            per_job.get(row["rule"], {}).get("min_requested_runtime_min"),
            per_job.get(row["rule"], {}).get("jobs_at_min_request"),
            _minutes(row["elapsed_sec"]),
            _efficiency(row["runtime_efficiency"]),
            _efficiency(per_job.get(row["rule"], {}).get("runtime_p95")),
            per_job.get(row["rule"], {}).get("total_elapsed_sec", 0) / 3600.0,
            (row["allocated_cpu_sec"] or 0) / 3600.0,
            (row["used_cpu_sec"] or 0) / 3600.0,
            _percent(row["used_cpu_sec"], row["allocated_cpu_sec"]),
        )
        for row in rows
    ]

    # Two different runtime questions get asked of this table and read as one
    # number if they share a card: how long a single job of the rule takes,
    # which is what the runtime request has to cover, and how much wall time
    # the rule consumed altogether, which is what says whether tuning it is
    # worth the effort at all. They are given separate, explicitly named cards.
    total_jobs = sum(entry["jobs"] for entry in per_job.values())
    total_elapsed_sec = sum(entry["total_elapsed_sec"] for entry in per_job.values())
    mean_job_minutes = (
        total_elapsed_sec / total_jobs / 60.0 if total_jobs else None
    )
    # The three column means only say something once there is more than one
    # rule to average, but the two runtime totals hold for a single rule too.
    averages = [
        ("Mean memory used %", _cards_percent(_mean(table_rows, 7))),
        ("Mean runtime used %", _cards_percent(_mean(table_rows, 14))),
        ("Mean CPU efficiency %", _cards_percent(_mean(table_rows, 19))),
    ] if len(table_rows) >= MIN_STAT_ROWS else []
    highlights = _highlights([
        *averages,
        (
            "Mean runtime per job",
            f"{_ONE(mean_job_minutes)} min" if mean_job_minutes is not None else None,
            f"averaged over {_quantity(total_jobs, 'job')}",
        ),
        (
            "Total runtime, all jobs",
            f"{_TWO(total_elapsed_sec / 3600.0)} h" if total_elapsed_sec else None,
            "wall time summed across every job",
        ),
    ])

    blocks = [
        _heading(
            "Requested versus used resources, per rule",
            "One row per workflow step, as medians across its launches. "
            "Requested figures are what the workflow asked the scheduler for; "
            "used figures are what the scheduler's accounting reports. "
            "Percentages well below 100 mean the step reserved far more than it "
            "needed, which wastes queue time; percentages near 100 mean it ran "
            "close to its limit and may fail on a larger dataset. The 95th "
            "percentile columns are the ceiling the heaviest jobs of the rule "
            "reached, and are what a request has to cover — a rule whose median "
            "is 20% and whose 95th percentile is 95% is not over-provisioned."
        ),
        _note(
            "The runtime floor is the smallest runtime the rule can ask for, "
            "read from the rule definitions of the installed Drakkar: requests "
            "are scaled from the input size but never drop below it. Where the "
            "smallest request equals the floor — times the run's time "
            "multiplier, if one was set — and the share of jobs sitting at it "
            "is high, the rule is priced by its floor and not by its size "
            "coefficient, so lowering the coefficient will not shorten a "
            "single request. The floor is blank for rules whose runtime is not "
            "written as a floored expression, and for the few rule names that "
            "several workflow modules define with different floors."
        ),
        highlights,
        _table(
            [
                ("Rule", _text),
                ("Jobs", _integer),
                ("Failed", _integer),
                ("CPUs requested", _ONE),
                ("CPUs allocated", _ONE),
                ("Memory requested (MB)", _integer),
                ("Peak memory (MB)", _integer),
                ("Memory used %, median", _ONE),
                ("Memory used %, 95th pct", _ONE),
                ("Runtime requested (min)", _ONE),
                ("Runtime floor (min)", _integer),
                ("Smallest runtime request (min)", _ONE),
                ("Jobs at smallest request %", _ONE),
                ("Runtime per job (min), median", _ONE),
                ("Runtime used %, median", _ONE),
                ("Runtime used %, 95th pct", _ONE),
                ("Total runtime (h)", _TWO),
                ("CPU-hours allocated", _TWO),
                ("CPU-hours used", _TWO),
                ("CPU efficiency %", _ONE),
            ],
            table_rows,
        ),
    ]
    figure = _rule_efficiency_figure(rows, per_job)
    if figure is not None:
        blocks.append(figure)
        blocks.append(_note(
            "The bar is the rule's median job; the whisker runs out to its "
            "95th percentile. A short bar with a long whisker is a rule whose "
            "request is sized for its typical job and is being carried by its "
            "heaviest ones — raise the request towards the whisker's end. A "
            "short bar with no whisker is simply over-provisioned throughout, "
            "and the request can come down to a little above the bar."
        ))
    return blocks


def _rule_efficiency_figure(rows, per_job):
    """Memory and runtime headroom for the rules that burn the most CPU time.

    Each bar is the rule's median job and each whisker reaches its 95th
    percentile, because those are the two figures a request has to be set
    between: the median says how much of the reservation a normal job leaves
    unused, and the 95th percentile says how far the request can actually be
    cut before the rule's heaviest jobs start failing.
    """
    import plotly.graph_objects as go

    ranked = [
        row for row in rows[:TOP_CATEGORIES]
        if row["memory_efficiency"] is not None or row["runtime_efficiency"] is not None
    ]
    if not ranked:
        return None
    ranked = list(reversed(ranked))
    rules = [row["rule"] for row in ranked]

    def whisker(medians, ceilings):
        """The span from each median out to its 95th percentile, never negative."""
        return [
            max(0.0, ceiling - median)
            if median is not None and ceiling is not None else 0.0
            for median, ceiling in zip(medians, ceilings)
        ]

    figure = go.Figure()
    for name, column, key, colour in (
        ("Memory used %", "memory_efficiency", "memory_p95", PALETTE[0]),
        ("Runtime used %", "runtime_efficiency", "runtime_p95", PALETTE[1]),
    ):
        medians = [_efficiency(row[column]) for row in ranked]
        ceilings = [
            _efficiency(per_job.get(row["rule"], {}).get(key)) for row in ranked
        ]
        figure.add_bar(
            name=name,
            y=rules,
            x=medians,
            orientation="h",
            marker_color=colour,
            error_x=dict(
                type="data",
                symmetric=False,
                array=whisker(medians, ceilings),
                arrayminus=[0.0] * len(rules),
                color="#22303c",
                thickness=1.2,
                width=4,
            ),
            customdata=[
                [0.0 if value is None else value] for value in ceilings
            ],
            hovertemplate=(
                f"%{{y}}<br>{name} — median %{{x:.1f}}%"
                "<br>95th percentile %{customdata[0]:.1f}%<extra></extra>"
            ),
        )
    figure.update_layout(
        barmode="group",
        xaxis_title="Percent of the requested resource actually used "
                    "(bar: median job, whisker: 95th percentile)",
        yaxis_title="",
        legend_title_text="",
        height=max(FIGURE_HEIGHT, 22 * len(rules) + 160),
    )
    return figure


def _job_resource_blocks(connection):
    """The individual launches, heaviest first.

    Every launch is listed. ``benchmark_job`` holds one row per submitted job,
    so it is bounded by how much the run was split up, not by the size of the
    sequencing data.
    """
    rows = _query(connection, """
        SELECT run_id, rule, wildcards, attempt, state,
               requested_cpus, alloc_cpus,
               requested_mem_mb, max_rss_mb, memory_efficiency,
               requested_runtime_min, elapsed_sec, runtime_efficiency,
               cpu_efficiency
        FROM benchmark_job
        ORDER BY COALESCE(elapsed_sec, -1) DESC, run_id, launch_index
    """)
    if not rows:
        return []

    blocks = [
        _heading(
            "Requested versus used resources, per job",
            "The same comparison for individual jobs, which is where a "
            "single oversized sample shows up as one job far heavier than "
            "the rest of its rule."
        ),
        _paragraph(
            f"{_quantity(len(rows), 'submitted job')}, longest-running first."
        ),
        _table(
            [
                ("Run", _text),
                ("Rule", _text),
                ("Wildcards", _text),
                ("Attempt", _integer),
                ("State", _text),
                ("CPUs requested", _integer),
                ("CPUs allocated", _integer),
                ("Memory requested (MB)", _integer),
                ("Peak memory (MB)", _integer),
                ("Memory used %", _ONE),
                ("Runtime requested (min)", _ONE),
                ("Runtime used (min)", _ONE),
                ("Runtime used %", _ONE),
                ("CPU efficiency %", _ONE),
            ],
            [
                (
                    row["run_id"], row["rule"], row["wildcards"], row["attempt"],
                    row["state"], row["requested_cpus"], row["alloc_cpus"],
                    row["requested_mem_mb"], row["max_rss_mb"],
                    _efficiency(row["memory_efficiency"]),
                    row["requested_runtime_min"], _minutes(row["elapsed_sec"]),
                    _efficiency(row["runtime_efficiency"]),
                    _efficiency(row["cpu_efficiency"]),
                )
                for row in rows
            ],
        ),
    ]
    return blocks


def _render_provenance(connection):
    """The ingest log, so every number on the page can be traced to a file."""
    rows = _query(connection, """
        SELECT table_name, section, rows_ingested, source_file, ingested_at
        FROM ingest_log
        ORDER BY section, table_name
    """)
    if not rows:
        return ""
    return "".join([
        '<section class="panel" id="section-provenance">',
        "<h2>Provenance</h2>",
        _paragraph(
            "Every table in the report database, the file it was read from, and "
            "when it was ingested."
        ),
        _table(
            [
                ("Table", _text),
                ("Section", _text),
                ("Rows", _integer),
                ("Source file", _text),
                ("Ingested", _timestamp),
            ],
            [tuple(row) for row in rows],
        ),
        "</section>",
    ])


SECTION_RENDERERS = {
    "preprocessing": _render_preprocessing,
    "cataloging": _render_cataloging,
    "dereplication": _render_dereplication,
    "profiling": _render_profiling,
    "taxonomy": _render_taxonomy,
    "function": _render_function,
    "expression": _render_expression,
    "amr": _render_amr,
    "resources": _render_resources,
}


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------

def _summary_html(connection, db_path, rendered, skipped, not_selected):
    """The sidebar's report details: versions, runs, and section outcome."""
    stamp = _query(
        connection,
        "SELECT version, drakkar_version, created_at FROM schema_version "
        "ORDER BY rowid DESC LIMIT 1",
    )
    schema_version = stamp[0]["version"] if stamp else SCHEMA_VERSION
    drakkar_version = stamp[0]["drakkar_version"] if stamp else None

    run_ids = [row[0] for row in _query(connection, "SELECT run_id FROM run ORDER BY run_id")]

    entries = [
        ("Drakkar version", drakkar_version or "unknown"),
        ("Report schema", f"version {schema_version}"),
        ("Database", str(db_path)),
        ("Runs", ", ".join(run_ids) if run_ids else "no run metadata recorded"),
        ("Sections rendered", ", ".join(SECTION_LABELS[name] for name in rendered)),
    ]
    if skipped:
        entries.append((
            "Sections unavailable",
            ", ".join(SECTION_LABELS[name] for name in skipped),
        ))
    if not_selected:
        entries.append((
            "Sections not selected",
            ", ".join(SECTION_LABELS[name] for name in not_selected),
        ))

    items = "".join(
        f"<dt>{escape(label)}</dt><dd>{_text(value)}</dd>" for label, value in entries
    )
    return f'<dl class="summary">{items}</dl>'


def _skipped_html(skipped, not_selected):
    if not skipped and not not_selected:
        return ""
    parts = []
    if skipped:
        items = "".join(
            f"<li>{escape(SECTION_LABELS[name])} — not present in the report database</li>"
            for name in skipped
        )
        parts.append(
            "<p><strong>Sections omitted.</strong> Drakkar is modular, so these "
            "sections have no data in this output directory:</p>"
            f"<ul>{items}</ul>"
        )
    if not_selected:
        items = "".join(
            f"<li>{escape(SECTION_LABELS[name])}</li>" for name in not_selected
        )
        parts.append(
            "<p>These sections were excluded by <code>--sections</code>:</p>"
            f"<ul>{items}</ul>"
        )
    return '<div class="skipped">' + "".join(parts) + "</div>"


def _toc_html(targets):
    """The sidebar navigation; the first entry is the panel that opens first."""
    if not targets:
        return ""
    links = []
    for index, (target, label) in enumerate(targets):
        active = ' class="is-active"' if index == 0 else ""
        links.append(f'<li><a href="#{target}"{active}>{escape(label)}</a></li>')
    return '<nav class="toc side-block"><h2>Contents</h2><ol>' + "".join(links) + "</ol></nav>"


def _sidebar_html(connection, db_path, rendered, skipped, not_selected,
                  targets, generated):
    """Everything about the report itself, kept beside the sections it describes."""
    return "".join([
        '<aside class="sidebar">',
        '<div class="brand">',
        _logo_html(),
        "<h1>Analysis Report</h1>",
        f'<p class="subtitle">Rendered {escape(generated)}</p>',
        "</div>",
        _toc_html(targets),
        '<section class="side-block"><h2>Report details</h2>',
        _summary_html(connection, db_path, rendered, skipped, not_selected),
        "</section>",
        _skipped_html(skipped, not_selected),
        "</aside>",
    ])


def render_report(db_path, html_path, sections=None):
    """Render the HTML report from a report database.

    ``sections`` is the selection already resolved by the command; sections
    that the database does not hold are reported on the page instead of
    raising. Returns a dictionary with the ``rendered``, ``skipped`` and
    ``not_selected`` section names.
    """
    from drakkar import __version__

    requested = tuple(sections) if sections is not None else tuple(SECTION_RENDERERS)
    connection = sqlite3.connect(str(db_path))
    connection.row_factory = sqlite3.Row
    try:
        available = present_sections(connection)
        state = {"plotly_embedded": False}

        rendered = []
        skipped = []
        bodies = []
        for name in requested:
            if name not in available:
                skipped.append(name)
                continue
            blocks = SECTION_RENDERERS[name](connection)
            if not blocks:
                skipped.append(name)
                continue
            rendered.append(name)
            bodies.append(
                f'<section class="panel" id="section-{name}">'
                f"<h2>{escape(SECTION_LABELS[name])}</h2>"
                + _section_intro(name)
                + _serialize(blocks, state)
                + "</section>"
            )

        not_selected = [
            name for name in SECTION_RENDERERS if name not in requested
        ]
        provenance = _render_provenance(connection)
        generated = _friendly_datetime(datetime.now(timezone.utc))

        panels = list(bodies)
        targets = [(f"section-{name}", SECTION_LABELS[name]) for name in rendered]
        if provenance:
            panels.append(provenance)
            targets.append(("section-provenance", "Provenance"))
        if not panels:
            panels.append(
                '<section class="panel">'
                + _paragraph(
                    "The report database holds none of the requested sections, "
                    "so there is nothing to show."
                )
                + "</section>"
            )
        # Whichever panel comes first is the one the page opens on, and the
        # first entry of the table of contents points at it.
        panels[0] = panels[0].replace(
            '<section class="panel"', '<section class="panel is-active"', 1
        )

        document = "\n".join([
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1">',
            "<title>Drakkar report</title>",
            f"<style>{STYLESHEET}</style>",
            # Set before the body renders, so panels are never shown and then
            # hidden; without scripting the class is absent and all of them stay.
            "<script>document.documentElement.className += ' js';</script>",
            "</head>",
            "<body>",
            '<div class="layout">',
            _sidebar_html(connection, db_path, rendered, skipped, not_selected,
                          targets, generated),
            '<main class="content">',
            "\n".join(panels),
            '<footer class="report">',
            f"Generated by Drakkar {escape(str(__version__))} from the report "
            f"database at {escape(str(db_path))}; no source table was read at "
            f"render time.",
            "</footer>",
            "</main>",
            "</div>",
            f"<script>{SCRIPT}</script>",
            "</body>",
            "</html>",
        ])
    finally:
        connection.close()

    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(document, encoding="utf-8")
    return {"rendered": rendered, "skipped": skipped, "not_selected": not_selected}
