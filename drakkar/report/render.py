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

The page is self-contained: the stylesheet is inlined, and the Plotly bundle is
embedded once, in the first figure, so the report opens offline.
"""

import sqlite3
from datetime import datetime, timezone
from html import escape

from drakkar.report.schema import SCHEMA_VERSION, TAXONOMIC_RANKS

# Row and category caps. Nothing on the page is allowed to grow with the size
# of the largest tables, so every listing is bounded here rather than in SQL
# scattered across the section renderers.
TABLE_ROW_LIMIT = 100
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
}
* { box-sizing: border-box; }
body {
  margin: 0;
  padding: 0 1.5rem 4rem;
  font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
  font-size: 15px;
  line-height: 1.55;
  color: var(--ink);
  background: #ffffff;
}
.wrap { max-width: 1080px; margin: 0 auto; }
header.report { border-bottom: 2px solid var(--accent); padding: 2.5rem 0 1.25rem; }
header.report h1 { margin: 0 0 .25rem; font-size: 1.9rem; letter-spacing: .01em; }
header.report .subtitle { margin: 0; color: var(--muted); }
h2 { margin: 2.75rem 0 .5rem; font-size: 1.35rem; padding-bottom: .3rem; border-bottom: 1px solid var(--rule); }
h3 { margin: 1.75rem 0 .4rem; font-size: 1.05rem; color: var(--ink); }
p { margin: .5rem 0; }
p.note { color: var(--muted); font-size: .88rem; }
nav.toc { margin: 1.5rem 0 0; padding: .9rem 1.1rem; background: var(--panel); border: 1px solid var(--rule); }
nav.toc h2 { margin: 0 0 .4rem; font-size: .8rem; text-transform: uppercase; letter-spacing: .08em; color: var(--muted); border: 0; padding: 0; }
nav.toc ol { margin: 0; padding-left: 1.2rem; }
nav.toc a { color: var(--ink); text-decoration: none; }
nav.toc a:hover { color: var(--accent); text-decoration: underline; }
dl.summary { display: grid; grid-template-columns: 13rem 1fr; gap: .2rem 1rem; margin: 1.2rem 0 0; }
dl.summary dt { color: var(--muted); font-size: .88rem; }
dl.summary dd { margin: 0; }
table { border-collapse: collapse; width: 100%; margin: .75rem 0; font-size: .9rem; }
.scroll { overflow-x: auto; }
th, td { padding: .35rem .6rem; text-align: right; border-bottom: 1px solid var(--rule); white-space: nowrap; }
th { background: var(--panel); font-weight: 600; color: var(--muted); text-transform: uppercase; font-size: .72rem; letter-spacing: .05em; }
th:first-child, td:first-child { text-align: left; }
tbody tr:hover td { background: #fbfcfd; }
.figure { margin: 1rem 0 .25rem; }
.skipped { margin: 1.2rem 0 0; padding: .8rem 1.1rem; background: var(--panel); border-left: 3px solid var(--rule); }
.skipped ul { margin: .3rem 0 0; padding-left: 1.2rem; color: var(--muted); font-size: .9rem; }
footer.report { margin-top: 3.5rem; padding-top: 1rem; border-top: 1px solid var(--rule); color: var(--muted); font-size: .85rem; }
code { font-family: "SFMono-Regular", Consolas, monospace; font-size: .85em; }
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


def _table(columns, rows, limit=TABLE_ROW_LIMIT):
    """Render rows as an HTML table, dropping columns that are entirely empty.

    ``columns`` is a sequence of ``(header, formatter)`` pairs and ``rows`` a
    sequence of raw value tuples in the same order.
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
        cells = "".join(f"<td>{columns[index][1](row[index])}</td>" for index in keep)
        body.append(f"<tr>{cells}</tr>")
    parts = [
        '<div class="scroll"><table>',
        f"<thead><tr>{head}</tr></thead>",
        "<tbody>" + "".join(body) + "</tbody>",
        "</table></div>",
    ]
    if len(rows) > limit:
        parts.append(
            f'<p class="note">Showing the first {limit:,} of {len(rows):,} rows; '
            f"query the database for the rest.</p>"
        )
    return "".join(parts)


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
    rows = _query(connection, """
        SELECT assembly_id, assembly_contigs, assembly_total_length,
               assembly_largest_contig, assembly_N50, assembly_gc_percent,
               mapping_rate_percent, final_bins, high_quality_bins,
               medium_quality_bins, low_quality_bins, bin_mean_completeness,
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
    blocks.append(_table(
        [
            ("Assembly", _text),
            ("Contigs", _integer),
            ("Total length", _integer),
            ("Largest contig", _integer),
            ("N50", _integer),
            ("GC %", _ONE),
            ("Mapping rate %", _ONE),
            ("Final bins", _integer),
            ("High quality", _integer),
            ("Medium quality", _integer),
            ("Low quality", _integer),
            ("Mean completeness", _ONE),
            ("Mean contamination", _ONE),
        ],
        [tuple(row) for row in rows],
    ))

    blocks.extend(_binner_blocks(connection))
    blocks.extend(_mapping_rate_blocks(connection))
    return blocks


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
        ),
        figure,
    ]


def _mapping_rate_blocks(connection):
    """Per-sample mapping rates against the assemblies they contributed to."""
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
        blocks.append(_table(
            [
                ("Metric", _text),
                ("Value", _text),
            ],
            [
                ("Genomes", _integer(row["genomes"])),
                ("Mean completeness", _ONE(row["completeness"])),
                ("Mean contamination", _ONE(row["contamination"])),
                ("Mean size", _integer(row["size"])),
                ("Mean N50", _integer(row["n50"])),
                ("Clustered members", _integer(row["members"])),
            ],
        ))
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
        blocks.append(_table(
            [("Metric", _text), ("Value", _text)],
            [
                ("Genes with a length", _integer(lengths[0]["genes"])),
                ("Mean length", _integer(lengths[0]["mean_length"])),
                ("Shortest", _integer(lengths[0]["min_length"])),
                ("Longest", _integer(lengths[0]["max_length"])),
            ],
        ))
    return blocks


def _render_resources(connection):
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
        '<h2 id="section-provenance">Provenance</h2>',
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
    """The header block: versions, runs, ingest window, and section outcome."""
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


def _toc_html(rendered, with_provenance):
    links = [
        f'<li><a href="#section-{name}">{escape(SECTION_LABELS[name])}</a></li>'
        for name in rendered
    ]
    if with_provenance:
        links.append('<li><a href="#section-provenance">Provenance</a></li>')
    if not links:
        return ""
    return (
        '<nav class="toc"><h2>Contents</h2><ol>' + "".join(links) + "</ol></nav>"
    )


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
                f'<h2 id="section-{name}">{escape(SECTION_LABELS[name])}</h2>'
                + _serialize(blocks, state)
            )

        not_selected = [
            name for name in SECTION_RENDERERS if name not in requested
        ]
        provenance = _render_provenance(connection)
        generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

        document = "\n".join([
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1">',
            "<title>Drakkar report</title>",
            f"<style>{STYLESHEET}</style>",
            "</head>",
            "<body>",
            '<div class="wrap">',
            '<header class="report">',
            "<h1>Drakkar report</h1>",
            f'<p class="subtitle">Rendered {escape(generated)} from '
            f"{escape(str(db_path))}</p>",
            "</header>",
            _summary_html(connection, db_path, rendered, skipped, not_selected),
            _skipped_html(skipped, not_selected),
            _toc_html(rendered, bool(provenance)),
            "<main>",
            "\n".join(bodies) if bodies
            else _paragraph(
                "The report database holds none of the requested sections, so "
                "there is nothing to show."
            ),
            provenance,
            "</main>",
            '<footer class="report">',
            f"Generated by Drakkar {escape(str(__version__))} from the report "
            f"database; no source table was read at render time.",
            "</footer>",
            "</div>",
            "</body>",
            "</html>",
        ])
    finally:
        connection.close()

    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(document, encoding="utf-8")
    return {"rendered": rendered, "skipped": skipped, "not_selected": not_selected}
