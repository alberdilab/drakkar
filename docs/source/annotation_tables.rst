.. _annotation-tables:

Annotation table reference
==========================

This page documents the Drakkar 2.0 functional annotation outputs, with an
emphasis on ``annotating/gene_annotations.tsv.xz`` and migration from the 1.x
wide gene table. The 2.0 table is a lossless, long-form evidence table: one
gene can have several rows because every accepted hit from every enabled
source is retained.

The essential rules are:

- Use ``(mag, gene)`` as the gene key.
- Use ``(mag, gene, source, hit_rank)`` as the row key.
- ``hit_rank`` starts at 1 independently for each ``(mag, gene, source)``.
- ``is_primary`` marks rank 1; lower-ranked accepted hits remain in the table.
- Every predicted gene has a ``source=prodigal`` row, including genes with no
  functional annotation.
- Exclude ``source=prodigal`` when counting functionally annotated genes.
- Empty TSV fields mean that the source did not supply that value. They do not
  mean zero.

Inspect the compressed table without extracting it permanently:

.. code-block:: console

   $ xz -dc annotating/gene_annotations.tsv.xz | head

Why the table changed
---------------------

The 1.x table had one row per gene and one column per database. It retained a
single selected value for fields such as ``kegg``, ``pfam``, and ``cazy``.
That layout was convenient for a quick spreadsheet view, but it discarded
additional qualifying hits and could not preserve source-specific scores,
alignment coordinates, or filter provenance without continually adding new
columns.

Drakkar 2.0 stores each accepted source hit as its own row. Shared identifiers
and coordinates stay in stable columns, while source-specific native data is
kept in the JSON ``details`` object. This makes the table larger vertically,
but prevents evidence loss and makes filtering and ranking auditable.

Before and after
----------------

The following mock rows show the structural difference. Values are
illustrative and only the most useful columns are displayed.

In 1.x, one row combined the selected annotations from several sources:

.. code-block:: text

   gene   start  end  strand  kegg    ec       pfam    cazy
   c1_1   1      900  +       K00001  1.1.1.1  PF00005 GH5

In 2.0, the same gene occupies one Prodigal row plus one row for every
accepted source hit. A second accepted KEGG hit is retained rather than
discarded:

.. code-block:: text

   mag    gene  contig  start  end  strand  source    hit_rank  is_primary  annotation_id  annotation
   MAG_A  c1_1  c1      1      900  +       prodigal  1         True        CDS
   MAG_A  c1_1  c1      1      900  +       kegg      1         True        K00001         1.1.1.1
   MAG_A  c1_1  c1      1      900  +       kegg      2         False       K00002         2.7.7.7
   MAG_A  c1_1  c1      1      900  +       pfam      1         True        PF00005        ABC_tran
   MAG_A  c1_1  c1      1      900  +       cazy      1         True        GH5            GH5

An unannotated gene still has one row:

.. code-block:: text

   mag    gene  contig  start  end   strand  source    hit_rank  is_primary  annotation_id
   MAG_A  c1_2  c1      1000   1600  -       prodigal  1         True        CDS

Ordered gene-table schema
-------------------------

The columns below appear in this exact order. TSV is a text format and does
not embed data types; the listed types describe how consumers should parse the
values. Numeric fields are blank when unavailable.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Column
     - Type or unit
     - Meaning
   * - ``mag``
     - String
     - MAG identifier written explicitly on every row.
   * - ``gene``
     - String
     - Prodigal-derived gene identifier, unique within ``mag``.
   * - ``contig``
     - String
     - Contig containing the gene.
   * - ``start``
     - Integer, 1-based
     - Inclusive gene start copied from the Prodigal GFF.
   * - ``end``
     - Integer, 1-based
     - Inclusive gene end copied from the Prodigal GFF.
   * - ``strand``
     - ``+`` or ``-``
     - Gene strand copied from the Prodigal GFF.
   * - ``source``
     - String vocabulary
     - Annotation database or predictor; see `Source values`_.
   * - ``method``
     - String vocabulary
     - Tool or procedure that generated the hit.
   * - ``evidence``
     - String vocabulary
     - Broad evidence class, such as sequence homology or gene prediction.
   * - ``hit_rank``
     - Positive integer
     - Rank within one ``(mag, gene, source)`` group.
   * - ``is_primary``
     - Boolean
     - ``True`` exactly when ``hit_rank`` is 1.
   * - ``rank_score``
     - Source-specific numeric
     - Value used to order accepted hits within a source.
   * - ``rank_score_type``
     - String
     - Defines the meaning and direction of ``rank_score``.
   * - ``annotation_id``
     - String
     - Principal source identifier, such as a KO, Pfam accession, or CAZy family.
   * - ``annotation``
     - String
     - Principal human-readable label or mapped function. It can be blank.
   * - ``annotation_type``
     - String
     - Source-specific class of annotation.
   * - ``evalue``
     - Floating point
     - Native expectation value when supplied.
   * - ``bitscore``
     - Floating point
     - Native full-sequence or alignment bit score when supplied.
   * - ``score``
     - Source-specific numeric
     - Native score used for the hit; interpret it with ``score_type``.
   * - ``score_type``
     - String
     - Unit or meaning of ``score``.
   * - ``threshold``
     - Source-specific numeric
     - Acceptance threshold in the unit identified by ``score_type``.
   * - ``identity``
     - Percent, 0--100
     - Sequence identity for supported alignment sources, currently VFDB.
   * - ``coverage``
     - Fraction, 0--1
     - Source-specific combined coverage: CAZy HMM coverage or the minimum of
       VFDB query and target coverage.
   * - ``query_coverage``
     - Fraction, 0--1
     - Fraction of the query covered by an alignment.
   * - ``target_coverage``
     - Fraction, 0--1
     - Fraction of the database target covered by an alignment.
   * - ``confidence``
     - Source-specific numeric
     - Native predictor confidence, currently used by SignalP.
   * - ``alignment_length``
     - Integer, residues
     - Native alignment length when supplied.
   * - ``query_start``
     - Integer, source-native
     - Alignment start on the query; for CAZy this is the protein target start.
   * - ``query_end``
     - Integer, source-native
     - Alignment end on the query; for CAZy this is the protein target end.
   * - ``target_start``
     - Integer, source-native
     - Alignment start on the database target.
   * - ``target_end``
     - Integer, source-native
     - Alignment end on the database target.
   * - ``model_start``
     - Integer, source-native
     - Start on an HMM or other model, currently populated for CAZy.
   * - ``model_end``
     - Integer, source-native
     - End on an HMM or other model, currently populated for CAZy.
   * - ``details``
     - JSON object
     - Remaining native fields, mappings, filter rules, and provenance. An
       empty object is written as ``{}``.

Do not compare ``score`` or ``rank_score`` across sources. For example, a KEGG
rank score can be a bit-score margin above a KOfam cutoff, whereas a CAZy rank
score is HMM coverage. ``score_type``, ``rank_score_type``, and ``threshold``
make those meanings explicit.

Source values
-------------

CLI target names and output ``source`` names are not always identical. The
following table lists every source produced by the supported 2.0 CLI.

.. list-table::
   :header-rows: 1
   :widths: 18 20 22 22 18

   * - CLI target
     - ``source``
     - ``method``
     - ``evidence``
     - Principal identifier
   * - Implicit gene calling
     - ``prodigal``
     - ``prodigal``
     - ``gene_prediction``
     - GFF feature type, normally ``CDS``
   * - ``kegg``
     - ``kegg``
     - ``hmmscan``
     - ``sequence_homology``
     - KOfam KO
   * - ``pfam``
     - ``pfam``
     - ``hmmscan``
     - ``sequence_homology``
     - Unversioned Pfam accession
   * - ``cazy``
     - ``cazy``
     - ``run_dbcan_hmm``
     - ``sequence_homology``
     - CAZy family
   * - ``virulence`` or ``vfdb``
     - ``vfdb``
     - ``mmseqs_easy_search``
     - ``sequence_homology``
     - VFDB entry
   * - ``amr``
     - ``ncbi_amrfinder``
     - ``hmmscan``
     - ``sequence_homology``
     - NCBIfam-AMRFinder accession
   * - ``signalp``
     - ``signalp``
     - ``signalp6``
     - ``protein_feature_prediction``
     - Signal-peptide class
   * - ``defense``
     - ``defensefinder``
     - ``defense_finder``
     - ``defense_system_prediction``
     - DefenseFinder gene or profile

Foldseek/ProstT5 structure annotation remains work in progress and is not an
available annotation target in Drakkar 2.0. Do not expect a supported 2.0 CLI
run to emit structure-hit rows.

Source-specific ``details``
---------------------------

``details`` is always valid compact JSON. Consumers should ignore unfamiliar
keys so that additional native provenance can be added without changing the
stable TSV columns. Important current keys include:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - ``source``
     - Important ``details`` contents
   * - ``prodigal``
     - Original GFF attributes, feature type, GFF source, and phase.
   * - ``kegg``
     - EC mappings, HMM description, accession, domain bit score, and HMMER
       overlap/region numbers.
   * - ``pfam``
     - All GOLD Pfam-to-EC associations, model name, original versioned
       accession, domain bit score, and the Pfam gathering-threshold rule.
   * - ``cazy``
     - HMM and target lengths, HMM file, and dbCAN coverage/e-value thresholds.
   * - ``vfdb``
     - All mapping records, VFC, query/target lengths, mismatch/gap counts, and
       identity, coverage, and e-value filter thresholds.
   * - ``ncbi_amrfinder``
     - All NCBIfam mapping records, HMM model metadata, domain bit score, and
       trusted-cutoff provenance.
   * - ``signalp``
     - Currently ``{}``; the stable score and confidence columns contain the
       native confidence.
   * - ``defensefinder``
     - The complete native DefenseFinder gene record, including ``activity``
       values used to distinguish defense from antidefense hits.

.. _migrating-gene-tables-2:

Migrating from the 1.x gene table
---------------------------------

Regenerate annotations after upgrading. Do not concatenate or append a 1.x
table to a 2.0 table, and do not treat this as a mechanical reshape of old
results. Drakkar 2.0 retains additional accepted hits and also changes several
scientific filters, including coverage-aware dbCAN CAZy calls and VFDB coverage
requirements. A fresh run is therefore the reproducible migration path.

The 1.x columns map as follows:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 1.x field
     - 2.0 equivalent or migration rule
   * - ``gene``
     - Remains ``gene``, but use ``(mag, gene)`` as the key. The new ``contig``
       column records the contig separately.
   * - ``start``, ``end``, ``strand``
     - Remain stable coordinate columns and repeat on every hit row.
   * - ``kegg``
     - Filter ``source == "kegg"`` and read ``annotation_id``. Use
       ``is_primary`` if one KO per gene is required.
   * - ``ec``
     - No longer a single gene-wide field. KEGG EC values are in ``annotation``
       and ``details.ec``; all Pfam associations are in
       ``details.ec_associations``. Preserve the source when extracting ECs.
   * - ``pfam``
     - Filter ``source == "pfam"`` and read ``annotation_id``.
   * - ``cazy``
     - Filter ``source == "cazy"`` and read ``annotation_id``. Multiple
       non-overlapping domains, including repeated families, can produce
       multiple rows.
   * - ``resistance_type``
     - Filter ``source == "ncbi_amrfinder"``. The mapped subtype or type is in
       ``annotation_type``.
   * - ``resistance_target``
     - Read the mapped subclass from the records in ``details.mappings`` on
       ``ncbi_amrfinder`` rows.
   * - ``vf``
     - Filter ``source == "vfdb"`` and read ``annotation``. ``annotation_id``
       contains the native VFDB entry.
   * - ``vf_type``
     - Read ``annotation_type`` on ``vfdb`` rows. Rebuild the VFDB mapping as
       described below because legacy mappings can contain an incorrect value.
   * - ``signalp``
     - Filter ``source == "signalp"`` and read ``annotation_id`` or
       ``annotation``; the predictor confidence is now retained explicitly.
   * - ``defense``
     - Filter ``source == "defensefinder"`` and
       ``annotation_type == "Defense"``. The former value is in
       ``annotation_id`` and the DefenseFinder system type is in ``annotation``.
   * - ``defense_type``
     - Read ``annotation`` on the corresponding DefenseFinder defense row.
   * - ``antidefense``
     - Filter ``source == "defensefinder"`` and
       ``annotation_type == "Antidefense"``. The former value is in
       ``annotation_id``.
   * - ``antidefense_type``
     - Read ``annotation`` on the corresponding DefenseFinder antidefense row.
   * - ``evidence``
     - The old gene-wide ``sequence``/``structure`` label becomes row-level
       ``method`` and ``evidence`` provenance. Structure annotation is not a
       supported target in 2.0.

Drakkar 1.x VFDB mappings must also be rebuilt. Install a fresh dated release
before rerunning virulence annotation:

.. code-block:: console

   $ drakkar database vfdb --directory /path/to/vfdb --set-default

Drakkar 2.0 writes ``mapping_schema=drakkar-vfdb-v2`` and refuses a legacy
mapping instead of silently emitting an incorrect virulence-factor type.

Common analysis recipes
-----------------------

Python/pandas
^^^^^^^^^^^^^

Read the compressed TSV directly, separate functional evidence from Prodigal
gene-presence rows, and select one representative hit per gene and source:

.. code-block:: python

   import pandas as pd

   annotations = pd.read_csv("annotating/gene_annotations.tsv.xz", sep="\t")
   functional = annotations.loc[annotations["source"] != "prodigal"].copy()
   primary = functional.loc[functional["is_primary"]].copy()

Count genes with at least one functional annotation:

.. code-block:: python

   annotated_gene_count = functional[["mag", "gene"]].drop_duplicates().shape[0]

Create a convenience wide view containing primary identifiers without
discarding rows from the archival long-form table:

.. code-block:: python

   wide = (
       primary.pivot(
           index=["mag", "gene", "contig", "start", "end", "strand"],
           columns="source",
           values="annotation_id",
       )
       .reset_index()
       .rename_axis(columns=None)
   )

R
^

.. code-block:: r

   annotations <- read.delim(
       xzfile("annotating/gene_annotations.tsv.xz"),
       check.names = FALSE,
       na.strings = ""
   )

   functional <- subset(annotations, source != "prodigal")
   primary <- subset(functional, is_primary)
   annotated_gene_count <- nrow(unique(functional[c("mag", "gene")]))

   wide <- reshape(
       primary[c("mag", "gene", "contig", "start", "end", "strand",
                 "source", "annotation_id")],
       idvar = c("mag", "gene", "contig", "start", "end", "strand"),
       timevar = "source",
       direction = "wide"
   )

Provenance and QC sidecars
--------------------------

Functional annotation runs also write:

- ``annotating/annotation_manifest.yaml``: Drakkar version, enabled sources,
  configured filters and tools, database paths/releases/checksums, environment
  dependencies, and checksums of the final annotation tables.
- ``annotating/annotation_qc.tsv``: per-MAG and per-source reported, retained,
  rejected, unmapped, and unique-entity counts plus the filtering stage. A
  blank rejected count means the upstream tool emitted only accepted calls,
  so the number rejected upstream is unavailable.

Keep both sidecars with the annotation tables when archiving or transferring
results. They provide the context needed to reproduce and audit the rows.
