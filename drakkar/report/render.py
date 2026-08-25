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

import sqlite3
from datetime import datetime, timezone
from html import escape

from drakkar.report.schema import SCHEMA_VERSION, TAXONOMIC_RANKS

# Row and category caps. Nothing on the page is allowed to grow with the size
# of the largest tables, so every listing is bounded here rather than in SQL
# scattered across the section renderers.
TABLE_ROW_LIMIT = 100
# Rows per page in the browser, and the smallest table worth averaging.
PAGE_ROWS = 20
MIN_STAT_ROWS = 2
TOP_GENOMES = 40
TOP_MAGS = 40
TOP_TAXA = 12
TOP_CATEGORIES = 30

FIGURE_HEIGHT = 420

# Muted palette, deliberately close to the terminal theme in drakkar.output.
PALETTE = (
    "#5f9ea0", "#d6a642", "#7fb069", "#e85d75", "#8878b0", "#4f7f9f",
    "#c98b5e", "#6fae9a", "#b0637f", "#8d9db6",
)

SECTION_LABELS = {
    "preprocessing": "Preprocessing",
    "cataloging": "Cataloging",
    "dereplication": "Dereplication",
    "profiling": "Profiling",
    "taxonomy": "Taxonomy",
    "function": "Functional annotation",
    "expression": "Expression",
    "resources": "Runs and resources",
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
  max-height: 100vh;
  overflow-y: auto;
  padding: 1.75rem 1.4rem 2.5rem;
  background: var(--panel);
  border-right: 1px solid var(--rule);
}
aside.sidebar .brand { border-bottom: 2px solid var(--accent); padding-bottom: .9rem; }
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

.figure { margin: 1rem 0 .25rem; }
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
    max-height: none;
    border-right: 0;
    border-bottom: 1px solid var(--rule);
  }
  main.content { padding: 1.5rem 1.25rem 3rem; }
}

/* Printing wants the whole report, not the panel that happens to be open. */
@media print {
  aside.sidebar { display: none; }
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


def _quantity(number, singular, plural=None):
    """Format a count with a noun that agrees with it."""
    noun = singular if number == 1 else (plural or singular + "s")
    return f"{number:,} {noun}"


def _ratio(numerator, denominator):
    """Percent of numerator over denominator, or None when it is undefined."""
    if numerator is None or not denominator:
        return None
    return 100.0 * float(numerator) / float(denominator)


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
    """Render ``(label, formatted value)`` pairs as a row of highlight cards."""
    cards = [
        f'<div class="stat"><span class="stat-label">{escape(label)}</span>'
        f'<span class="stat-value">{value}</span></div>'
        for label, value in items
        if value not in (None, "")
    ]
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


def _table(columns, rows, limit=TABLE_ROW_LIMIT, stats=()):
    """Render rows as an HTML table, dropping columns that are entirely empty.

    ``columns`` is a sequence of ``(header, formatter)`` pairs and ``rows`` a
    sequence of raw value tuples in the same order. ``stats`` names the columns
    whose average belongs above the table. Tables longer than ``PAGE_ROWS``
    carry the page size the script pages them by; the rows themselves are all
    in the markup, so the table stays complete without JavaScript.
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
    shown = rows[:limit]
    head = "".join(f"<th>{escape(columns[index][0])}</th>" for index in keep)
    body = []
    for row in shown:
        cells = "".join(_cell(row[index], columns[index][1]) for index in keep)
        body.append(f"<tr>{cells}</tr>")
    paging = f' data-page-size="{PAGE_ROWS}"' if len(shown) > PAGE_ROWS else ""
    parts = [
        '<div class="table-block">',
        _stats_html(columns, rows, stats),
        f'<div class="scroll"><table{paging}>',
        f"<thead><tr>{head}</tr></thead>",
        "<tbody>" + "".join(body) + "</tbody>",
        "</table></div>",
    ]
    if len(rows) > limit:
        parts.append(
            f'<p class="note">Showing the first {limit:,} of {len(rows):,} rows; '
            f"query the database for the rest.</p>"
        )
    parts.append("</div>")
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


def _heading(title):
    return f"<h3>{escape(title)}</h3>"


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
    """Serialize one Plotly figure, embedding the bundle only in the first."""
    import plotly.io as pio

    figure.update_layout(
        template="simple_white",
        font=dict(family="Helvetica Neue, Helvetica, Arial, sans-serif",
                  size=13, color="#22303c"),
        margin=dict(l=70, r=30, t=40, b=70),
        paper_bgcolor="white",
        plot_bgcolor="white",
        height=FIGURE_HEIGHT,
        colorway=list(PALETTE),
    )
    include = "inline" if not state["plotly_embedded"] else False
    state["plotly_embedded"] = True
    return '<div class="figure">' + pio.to_html(
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
        SELECT sample_id, reads_pre_fastp, reads_post_fastp, host_reads,
               metagenomic_reads, metagenomic_bases, singlem_fraction,
               nonpareil_diversity
        FROM sample
        ORDER BY sample_id
    """)
    if not rows:
        return []

    blocks = [_paragraph(
        f"{_quantity(len(rows), 'sample')} passed through quality filtering."
    )]

    table_rows = [
        (
            row["sample_id"],
            row["reads_pre_fastp"],
            row["reads_post_fastp"],
            row["host_reads"],
            row["metagenomic_reads"],
            _ratio(row["metagenomic_reads"], row["reads_pre_fastp"]),
            row["metagenomic_bases"],
            row["singlem_fraction"],
            row["nonpareil_diversity"],
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
            ("Metagenomic bases", _integer),
            ("SingleM fraction", _TWO),
            ("Nonpareil diversity", _TWO),
        ],
        table_rows,
        stats=[
            ("Mean reads in", 1),
            ("Mean metagenomic reads", 4),
            ("Mean metagenomic %", 5),
            ("Mean SingleM fraction", 7),
            ("Mean Nonpareil diversity", 8),
        ],
    ))

    figure = _read_fate_figure(rows)
    if figure is not None:
        blocks.append(_heading("Read fates"))
        blocks.append(figure)
    return blocks


def _read_fate_figure(rows):
    """Stacked read counts per sample: low quality, host, and metagenomic."""
    import plotly.graph_objects as go

    samples = [row["sample_id"] for row in rows]
    discarded = []
    for row in rows:
        before, after = row["reads_pre_fastp"], row["reads_post_fastp"]
        discarded.append(before - after if before is not None and after is not None else None)
    host = [row["host_reads"] for row in rows]
    metagenomic = [row["metagenomic_reads"] for row in rows]

    traces = [
        ("Removed by fastp", discarded),
        ("Host", host),
        ("Metagenomic", metagenomic),
    ]
    traces = [(name, values) for name, values in traces if any(v is not None for v in values)]
    if not traces:
        return None

    figure = go.Figure()
    for name, values in traces:
        figure.add_bar(name=name, x=samples, y=values)
    figure.update_layout(
        barmode="stack",
        xaxis_title="Sample",
        yaxis_title="Reads",
        legend_title_text="",
    )
    return figure


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
    blocks.extend(_bin_blocks(rows))
    blocks.extend(_binner_blocks(connection))
    return blocks


def _assembly_blocks(rows):
    """Contiguity and read recruitment, per assembly."""
    return [
        _heading("Assemblies"),
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
                ("Mean mapping rate %", 6),
            ],
        ),
    ]


def _bin_blocks(rows):
    """Bin yield and quality, per assembly.

    Low-quality bins are discarded downstream, so they are counted in
    ``final_bins`` but not given a column of their own.
    """
    return [
        _heading("Bins"),
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
        _heading("Bins per binner"),
        _table(
            [("Binner", _text), ("Bins", _integer), ("Assemblies", _integer)],
            [tuple(row) for row in rows],
            stats=[("Mean bins per binner", 1)],
        ),
        figure,
    ]


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
        _heading("Per-sample mapping rates"),
        _table(
            [
                ("Sample", _text),
                ("Assemblies", _integer),
                ("Mean rate %", _ONE),
                ("Min rate %", _ONE),
                ("Max rate %", _ONE),
            ],
            [tuple(row) for row in rows],
            stats=[("Mean mapping rate %", 2)],
        ),
    ]


def _render_dereplication(connection):
    import plotly.graph_objects as go

    rows = _query(connection, "SELECT * FROM dereplication")
    if not rows:
        return []
    row = rows[0]

    blocks = []
    ani = row["dereplication_ani"]
    if ani is not None:
        blocks.append(_paragraph(
            f"Bins were dereplicated at {float(ani):g} ANI."
        ))
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

    before, after = row["input_bin_number"], row["output_bin_number"]
    if before is not None and after is not None:
        figure = go.Figure()
        figure.add_bar(
            x=["Before", "After"],
            y=[before, after],
            marker_color=[PALETTE[0], PALETTE[1]],
        )
        figure.update_layout(xaxis_title="", yaxis_title="Bins")
        blocks.append(figure)
        if before:
            blocks.append(_note(
                f"Dereplication retained {_ratio(after, before):.1f}% of the input bins."
            ))
    return blocks


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
    return [_heading("Genome quality"), figure]


def _abundance_blocks(connection):
    """Relative abundance of the most abundant genomes across samples."""
    import plotly.graph_objects as go

    totals = _query(connection, """
        SELECT sample_id, SUM(read_count) AS reads
        FROM genome_count
        WHERE read_count IS NOT NULL
        GROUP BY sample_id
        ORDER BY sample_id
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
    blocks = [_heading("Genome abundance")]
    if total_genomes > len(genome_ids):
        blocks.append(_note(
            f"Showing the {len(genome_ids)} most abundant of {total_genomes:,} genomes, "
            f"as a percentage of each sample's mapped reads."
        ))
    blocks.append(figure)
    blocks.append(_table(
        [("Sample", _text), ("Mapped reads", _integer)],
        [(row["sample_id"], row["reads"]) for row in totals],
        stats=[("Mean mapped reads", 1)],
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
            marker_color=PALETTE[0],
        )
        figure.update_layout(
            xaxis_title="Genomes",
            yaxis_title="",
            yaxis=dict(autorange="reversed"),
            height=max(FIGURE_HEIGHT, 22 * len(phyla) + 140),
        )
        blocks.append(_heading("Genomes per phylum"))
        blocks.append(figure)

    blocks.extend(_composition_blocks(connection))
    return blocks


def _composition_blocks(connection):
    """Phylum-level composition per sample, when abundances are also present."""
    import plotly.graph_objects as go

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
            marker_color="#b7c7d3",
        )
    figure.update_layout(
        barmode="stack",
        xaxis_title="Sample",
        yaxis_title="Reads (% of sample)",
        legend_title_text="Phylum",
    )
    blocks = [_heading("Composition per sample")]
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
    blocks = [_heading("Annotation coverage per MAG")]
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
        _heading("Gene clusters and regions"),
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
        _heading("Annotation filtering"),
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
    blocks.append(_heading("Expression per sample"))
    blocks.append(figure)

    lengths = _query(connection, """
        SELECT COUNT(*) AS genes, AVG(length) AS mean_length,
               MIN(length) AS min_length, MAX(length) AS max_length
        FROM expressed_gene
        WHERE length IS NOT NULL
    """)
    if lengths and lengths[0]["genes"]:
        blocks.append(_heading("Quantified genes"))
        blocks.append(_highlights([
            ("Genes with a length", _integer(lengths[0]["genes"])),
            ("Mean length", _integer(lengths[0]["mean_length"])),
            ("Shortest", _integer(lengths[0]["min_length"])),
            ("Longest", _integer(lengths[0]["max_length"])),
        ]))
    return blocks


def _render_resources(connection):
    blocks = _run_blocks(connection)
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
            ("Started", _text),
            ("Finished", _text),
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


# Why a run carries no usage figures. Keyed by the status drakkar.benchmark
# stamps into drakkar_<run_id>_resources.yaml.
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
    return [_heading("Resource benchmark"), *notes]


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
        _heading("Job outcomes"),
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
            row["requested_runtime_min"],
            _minutes(row["elapsed_sec"]),
            _efficiency(row["runtime_efficiency"]),
            (row["allocated_cpu_sec"] or 0) / 3600.0,
            (row["used_cpu_sec"] or 0) / 3600.0,
            _percent(row["used_cpu_sec"], row["allocated_cpu_sec"]),
        )
        for row in rows
    ]
    blocks = [
        _heading("Requested versus used resources, per rule"),
        _paragraph(
            "Medians across the launches of each rule. Requested figures are "
            "what the workflow asked the scheduler for; used figures are what "
            "the scheduler's accounting reports, so the ratio between them is "
            "what a resource profile should be tuned against."
        ),
        _table(
            [
                ("Rule", _text),
                ("Jobs", _integer),
                ("Failed", _integer),
                ("CPUs requested", _ONE),
                ("CPUs allocated", _ONE),
                ("Memory requested (MB)", _integer),
                ("Peak memory (MB)", _integer),
                ("Memory used %", _ONE),
                ("Runtime requested (min)", _ONE),
                ("Runtime used (min)", _ONE),
                ("Runtime used %", _ONE),
                ("CPU-hours allocated", _TWO),
                ("CPU-hours used", _TWO),
                ("CPU efficiency %", _ONE),
            ],
            table_rows,
            stats=[
                ("Mean memory used %", 7),
                ("Mean runtime used %", 10),
                ("Mean CPU efficiency %", 13),
            ],
        ),
    ]
    figure = _rule_efficiency_figure(rows)
    if figure is not None:
        blocks.append(figure)
    return blocks


def _rule_efficiency_figure(rows):
    """Memory and runtime headroom for the rules that burn the most CPU time."""
    import plotly.graph_objects as go

    ranked = [
        row for row in rows[:TOP_CATEGORIES]
        if row["memory_efficiency"] is not None or row["runtime_efficiency"] is not None
    ]
    if not ranked:
        return None
    ranked = list(reversed(ranked))
    rules = [row["rule"] for row in ranked]

    figure = go.Figure()
    figure.add_bar(
        name="Memory used %",
        y=rules,
        x=[_efficiency(row["memory_efficiency"]) for row in ranked],
        orientation="h",
        marker_color=PALETTE[0],
    )
    figure.add_bar(
        name="Runtime used %",
        y=rules,
        x=[_efficiency(row["runtime_efficiency"]) for row in ranked],
        orientation="h",
        marker_color=PALETTE[1],
    )
    figure.update_layout(
        barmode="group",
        xaxis_title="Percent of the requested resource actually used",
        yaxis_title="",
        legend_title_text="",
        height=max(FIGURE_HEIGHT, 22 * len(rules) + 160),
    )
    return figure


def _job_resource_blocks(connection):
    """The individual launches, heaviest first."""
    rows = _query(connection, f"""
        SELECT run_id, rule, wildcards, attempt, state,
               requested_cpus, alloc_cpus,
               requested_mem_mb, max_rss_mb, memory_efficiency,
               requested_runtime_min, elapsed_sec, runtime_efficiency,
               cpu_efficiency
        FROM benchmark_job
        ORDER BY COALESCE(elapsed_sec, -1) DESC, run_id, launch_index
        LIMIT {TABLE_ROW_LIMIT + 1}
    """)
    if not rows:
        return []

    total = _scalar(connection, "SELECT COUNT(*) FROM benchmark_job", default=0)
    blocks = [
        _heading("Requested versus used resources, per job"),
        _paragraph(
            "One row per submitted job, longest-running first."
            if total <= TABLE_ROW_LIMIT else
            f"The {TABLE_ROW_LIMIT:,} longest-running of {total:,} submitted jobs; "
            "the full listing stays in benchmark/drakkar_<run_id>.jobs.tsv."
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
                for row in rows[:TABLE_ROW_LIMIT]
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
                ("Ingested", _text),
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
    "resources": _render_resources,
}


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------

def _summary_html(connection, db_path, rendered, skipped, not_selected):
    """The sidebar's report details: versions, runs, ingest window, outcome."""
    stamp = _query(
        connection,
        "SELECT version, drakkar_version, created_at FROM schema_version "
        "ORDER BY rowid DESC LIMIT 1",
    )
    schema_version = stamp[0]["version"] if stamp else SCHEMA_VERSION
    drakkar_version = stamp[0]["drakkar_version"] if stamp else None

    run_ids = [row[0] for row in _query(connection, "SELECT run_id FROM run ORDER BY run_id")]
    window = _query(
        connection, "SELECT MIN(ingested_at), MAX(ingested_at) FROM ingest_log"
    )
    first = window[0][0] if window else None
    last = window[0][1] if window else None

    entries = [
        ("Drakkar version", drakkar_version or "unknown"),
        ("Report schema", f"version {schema_version}"),
        ("Database", str(db_path)),
        ("Runs", ", ".join(run_ids) if run_ids else "no run metadata recorded"),
        ("Ingested", first if first == last or not last else f"{first} to {last}"),
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
        "<h1>Drakkar report</h1>",
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
                + _serialize(blocks, state)
                + "</section>"
            )

        not_selected = [
            name for name in SECTION_RENDERERS if name not in requested
        ]
        provenance = _render_provenance(connection)
        generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

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
