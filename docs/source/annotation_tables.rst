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
     - ``amrfinderplus``
     - ``sequence_homology``
     - AMRFinderPlus element symbol, such as ``blaOXA-48``
   * - ``card`` (alias: ``rgi``)
     - ``card``
     - ``rgi_main``
     - ``sequence_homology``
     - CARD ARO accession
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

.. _annotation-thresholds:

Acceptance thresholds and their evidence
----------------------------------------

Drakkar applies a two-tier acceptance policy. Where an annotation resource
ships curated, model-specific cutoffs, Drakkar defers to them and does not
layer a second global filter on top. Only where a resource provides no
calibrated cutoff does Drakkar apply its own global thresholds, which are
exposed as CLI options and recorded in ``annotation_manifest.yaml``.

Which rule applies to which source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 16 44 20 20

   * - ``source``
     - Acceptance rule
     - Applied by
     - Tunable
   * - ``kegg``
     - Per-KO KOfam bit score cutoff from ``ko_list``, compared against the
       full or domain score as the KO's ``score_type`` requires. KOs with no
       published cutoff fall back to ``--annotation-evalue``.
     - Drakkar, at merge
     - Fallback only
   * - ``pfam``
     - Pfam per-family gathering threshold (``hmmscan --cut_ga``).
     - Upstream, native
     - No
   * - ``cazy``
     - dbCAN HMM E-value < ``1e-15`` and HMM coverage > ``0.35``.
     - Upstream, native
     - No, fixed in the rule
   * - ``vfdb``
     - E-value, percent identity, query coverage and target coverage, all four
       required simultaneously.
     - Drakkar, at merge
     - Yes, all four
   * - ``ncbi_amrfinder``
     - AMRFinderPlus's own acceptance: per-gene curated identity and coverage
       cutoffs on its BLASTP arm, NCBIfam trusted cutoffs on its HMM arm.
     - Upstream, native
     - No
   * - ``card``
     - RGI's per-model curated bit score cutoffs. Only Perfect and Strict
       calls are emitted; Loose is off.
     - Upstream, native
     - No
   * - ``signalp``
     - SignalP 6 native model decision; every call it emits is retained.
     - Upstream, native
     - No
   * - ``defensefinder``
     - MacSyFinder profile GA scores plus system co-localisation rules.
     - Upstream, native
     - No
   * - ``genomad`` (cluster)
     - Virus score >= ``0.95`` and marker enrichment >= ``5``.
     - Drakkar, post-processing
     - No, fixed in the rule
   * - ``antismash`` (cluster)
     - antiSMASH native rule-based cluster detection.
     - Upstream, native
     - No
   * - ``dbcan`` (cluster)
     - ``run_dbcan`` CGC-finder defaults.
     - Upstream, native
     - No

The ``filter_stage`` column of ``annotation_qc.tsv`` separates the two tiers:
``upstream_native`` means the tool emitted only accepted calls, so Drakkar
cannot report how many were rejected, and ``drakkar`` means Drakkar applied
the filter during merging and counted both sides.

Why per-model cutoffs are preferred over a global E-value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The score that separates true from false members of a profile HMM family is
family-specific. A broad, deeply sampled, highly conserved family and a short,
shallow, divergent one do not share a usable cutoff, so a single E-value
applied across a whole profile database is simultaneously too permissive for
some families and too strict for others. Every profile database Drakkar uses
publishes its own solution to this, and Drakkar uses it rather than
substituting one number.

**KOfam / KEGG.** Aramaki *et al.* (2020) derive one threshold per KO by
maximising the F-measure against positive and negative training sets,
repeating the split three times and averaging the result. Benchmarked on 20
prokaryotic genomes, KofamScan reached F = 0.875, comparable to GhostKOALA
(0.886) and BlastKOALA (0.846) and clearly above KAAS (0.786). The cutoff
table also records *which* score to compare: in a recent ``ko_list`` release,
8,722 of 26,530 thresholded KOs are ``domain``-scored rather than
``full``-scored, and Drakkar selects the matching bit score per KO.

Because acceptance is by bit score, Drakkar deliberately runs ``hmmscan`` with
``-E 10 --domE 10`` for KEGG. That permissive *reporting* boundary exists so
that no hit is discarded by HMMER before its KO-specific cutoff can be
applied. It is not an acceptance threshold, and the raw
``annotating/kegg/<mag>.tsv`` files are not filtered output.

Per-KO reliability is not uniform, and the ``ko_list`` F-measure column makes
that measurable. In the same release the median per-KO F-measure is 0.986, but
9.4% of thresholded KOs score below 0.80 and 1.4% below 0.50. A KO assignment
is therefore only as good as its model; the ``score`` and ``threshold``
columns in the gene table let you re-filter on the margin above cutoff when a
downstream claim depends on a specific KO.

**Pfam.** Pfam curators set a gathering threshold (GA) per family, and every
sequence region scoring above it enters the family's full alignment. Mistry
*et al.* (2021) state that per-family gathering thresholds yield fewer false
positive matches than a single E-value threshold applied across all Pfam HMMs.
Drakkar uses ``--cut_ga`` and adds nothing.

**AMRFinderPlus.** The ``amr`` target runs AMRFinderPlus itself in combined
nucleotide/protein/GFF mode over the MAG's Prodigal calls, so both of its
detection arms contribute. The BLASTP arm searches the Pathogen Detection
Reference Gene Catalog under per-gene curated identity and coverage cutoffs,
and the HMM arm applies the NCBIfam trusted cutoffs (``--cut_tc``) that NCBI
maintains per model (Feldgarden *et al.*, 2021). Drakkar adds no threshold of
its own and filters only on element type, keeping this source to ``AMR`` and
leaving the stress and virulence "plus" genes out; the dropped rows are
counted in ``annotation_qc.tsv``.

The ``method`` value in ``details`` records which arm produced a call, and
``hit_rank`` follows AMRFinderPlus's own evidence ranking rather than raw
alignment score: ``ALLELE > EXACT > BLAST > INTERNAL_STOP > PARTIAL_CONTIG_END
> PARTIAL > HMM``. An HMM-only call is the weakest tier, so a gene whose only
evidence is ``method=HMM`` is a family-level assignment, not an allele call.

Why both arms matter: in release ``2026-08-07.1`` the NCBIfam-AMR HMM library
holds 784 models while the Reference Gene Catalog holds 10,078 proteins across
8,214 hierarchy nodes, only 410 of which have a matching HMM node. The HMM
library is family-level by design, so those counts are not a like-for-like
ratio, but whole allele series — most ``aac(3)`` variants and many ``bla``,
``tet``, ``erm``, ``sul``, ``dfr``, ``mcr`` and ``qnr`` alleles — are
detectable only through the catalog. Drakkar releases before 2.6 searched the
HMM library alone and recovered substantially fewer AMR genes as a result.

.. note::

   Point mutations are still out of scope for this target. AMRFinderPlus
   detects resistance-conferring mutations, such as in ``gyrA`` or ``rpoB``,
   only under ``--organism``, whose vocabulary is a short list of clinical
   taxa. Deriving it from GTDB-Tk taxonomy would be unreliable for most MAGs,
   so Drakkar omits the flag here and the ``amr`` source reports acquired
   genes only. For point mutations, run the :ref:`amr module <amr-workflow>`,
   whose manifest carries an explicit ``organism`` per assembly.

**CARD / RGI.** The ``card`` target runs RGI in protein mode over the same
Prodigal proteins, so its calls join the gene table on the same gene key. RGI
decides acceptance with per-model curated bit score cutoffs and reports the
tier it used in ``Cut_Off``: Perfect is an exact match to a curated reference
and Strict is above the model's cutoff. Drakkar does not pass
``--include_loose``, so sub-cutoff Loose calls are never emitted. The curated
cutoff itself is kept in the ``threshold`` column alongside the observed
``bitscore``, which makes the margin above cutoff directly auditable.

CARD and the NCBI Reference Gene Catalog are independently curated and do not
share an ontology, so the two sources are deliberately kept as separate rows
rather than reconciled. Agreement between them is evidence; disagreement is
informative, and collapsing them would destroy both signals. For coordinate-
reconciled loci across the two callers, use the :ref:`amr module
<amr-workflow>`.

**dbCAN / CAZy.** ``E-value < 1e-15`` and ``coverage > 0.35`` are
``run_dbcan``'s own defaults for the HMM method, established empirically in
the dbCAN papers and used by the dbCAN web server as the general annotation
standard. Note that dbCAN's organism-specific advice is stricter for bacteria
(``E < 1e-18``); Drakkar keeps the general default, which favours recall.

**DefenseFinder.** MacSyFinder v2 calls HMMER with ``--cut_ga`` when a profile
carries a GA score, which supersedes the i-evalue and profile-coverage
defaults, and defense systems must additionally satisfy the model's
co-localisation and quorum rules. Drakkar retains every gene DefenseFinder
reports and preserves the native record, including ``activity``, in
``details``.

**geNomad.** Drakkar keeps geNomad regions only at virus score >= 0.95 and
marker enrichment >= 5. For comparison, geNomad's own defaults are a score of
0.70 with no marker-enrichment requirement, and its ``--conservative`` preset
uses 0.80 and 1.50. Drakkar is therefore stricter than geNomad's most
stringent published preset. This is deliberate for MAG-level provirus calling,
where a false provirus call contaminates the host genome's functional profile,
but it is recall-limiting: expect Drakkar to miss proviruses that geNomad's
defaults would report.

Drakkar's own global thresholds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These apply only where no curated cutoff exists.

``--annotation-evalue``, default ``1e-10``
""""""""""""""""""""""""""""""""""""""""""

Applied to VFDB hits, to KOfam hits whose KO has no published cutoff, to
Foldseek structural hits, and to the HMMER hits that decide which genes are
"orphans" for the structural-annotation step.

Pearson (2013) reports that for protein-protein searches, expectation values
below 0.001 can reliably be used to infer homology. Drakkar's ``1e-10`` is
seven orders of magnitude more conservative than that boundary, which is
appropriate given that these are exactly the cases with no curated cutoff to
fall back on. It is also the value most commonly adopted in metagenomic
annotation practice.

The no-cutoff KO case is not marginal: in a recent ``ko_list`` release, 1,857
of 28,387 KOs (6.5%) carry no published threshold. anvi'o excludes such
"no-threshold KOs" from annotation entirely unless ``--include-nt-KOs`` is
given, and its relaxation heuristic for below-threshold hits uses an E-value
of ``1e-5``. Drakkar is therefore more inclusive than anvi'o for these KOs,
but gates them five orders of magnitude more strictly than anvi'o's heuristic.
These rows are identifiable and separable in the gene table: they carry
``score_type = evalue``, ``threshold`` equal to the configured e-value, and
``rank_score_type = negative_log10_evalue``, whereas cutoff-backed KO rows
carry ``full_bitscore`` or ``domain_bitscore`` and
``rank_score_type = bitscore_above_kofam_cutoff``.

``--annotation-identity``, default ``50`` percent
"""""""""""""""""""""""""""""""""""""""""""""""""

Applied to VFDB hits only. Rost (1999) found that above roughly 30% identity,
90% of aligned pairs are homologous, while below 25% fewer than 10% are, with
20-35% constituting the twilight zone. A 50% cutoff therefore sits clearly in
the safe zone for homology inference. For *function* transfer the requirement
is higher: Tian and Skolnick (2003) report that 40% identity supports
transferring the first three digits of an EC number, while all four digits
need above 60% identity for at least 90% accuracy.

Read the default accordingly. At 50% identity a ``vfdb`` row is well-supported
evidence that the gene is homologous to a known virulence factor and belongs
to that factor's class, and it is *not* sufficient evidence that the gene is
that specific virulence gene.

This is more permissive than common VFDB screening practice, which typically
uses 80% identity: ABRicate's default is 80%, EFSA guidance uses >80% identity
with >70% coverage, and many MAG studies use 80/80 or 90/80. Drakkar's
default is deliberately recall-oriented, because the long-form table retains
every accepted hit with its ``identity`` value, so tightening after the fact
costs nothing. If a result depends on a specific virulence gene being present,
either rerun with ``--annotation-identity 80`` or filter the table on
``identity >= 80``.

``--annotation-query-coverage`` and ``--annotation-target-coverage``, default ``0.5``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Applied to VFDB hits only, and both must be satisfied. Requiring coverage on
both sides is stricter than the single-sided coverage most screening tools
apply, and it targets the specific failure mode that a flat identity threshold
cannot catch. Rost (1999) notes that whether sequence identity implies
structural similarity depends crucially on alignment length: 10 identical
residues in an alignment of 16 is above 60% identity and implies nothing.
Requiring both coverages rejects a short conserved-motif match against a long
virulence protein (low target coverage) and a query whose virulence-like
region is a minority of its length (low query coverage), which is the usual
source of spurious hits in multi-domain proteins.

Note the two scales differ, and the gene table preserves both as reported:
``query_coverage`` and ``target_coverage`` are fractions from 0 to 1 (MMseqs2
``qcov``/``tcov``), while ``identity`` is a percentage from 0 to 100 (MMseqs2
``pident``). The ``coverage`` column holds the minimum of the two coverages
and is what ``vfdb`` rows are ranked on.

Choosing different thresholds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Goal
     - Setting
   * - Claim a specific virulence gene is present in a MAG
     - ``--annotation-identity 80 --annotation-query-coverage 0.8
       --annotation-target-coverage 0.8``, or filter the table equivalently.
   * - Survey virulence-factor classes across many MAGs
     - Keep the defaults; the class assignment in ``annotation_type`` is
       supported at 50% identity where the specific gene call is not.
   * - Exclude KOs that have no curated KOfam cutoff
     - Filter out ``source == "kegg" and score_type == "evalue"``, matching
       anvi'o's default behaviour.
   * - Restrict to high-confidence KO models
     - Keep rows where ``score - threshold`` is comfortably positive; the
       margin is already stored as ``rank_score`` for cutoff-backed KO rows.
   * - Make VFDB acceptance stricter than the safe zone of homology
     - Raise ``--annotation-evalue`` past ``1e-10`` only together with
       identity and coverage; e-value alone does not constrain partial hits.

Limitations to keep in mind
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **E-values are database-size dependent.** The same ``1e-10`` is not equally
  strict against VFDB and against AlphaFold/Swiss-Prot, because expectation
  values scale with the size of the searched database. Compare e-values within
  a source, not across sources.
- **Identity and coverage affect VFDB only.** Despite the general wording of
  ``--annotation-identity``, no other enabled source is filtered on identity
  or coverage by Drakkar. Foldseek rows carry an identity value but are gated
  on e-value alone, and their identity is a fraction rather than a percentage.
- **dbCAN and geNomad thresholds are not exposed.** They are fixed in the
  workflow rules. Changing them requires editing
  ``drakkar/workflow/rules/annotating_function.smk``, and the manifest will
  then no longer describe what was actually run.
- **AMR point mutations are not reported.** The ``amr`` and ``card`` targets
  report acquired resistance genes. A MAG with no AMR rows may still carry
  resistance-conferring mutations in core genes.
- **Per-KO model quality varies.** Roughly one KO in ten has an F-measure
  below 0.80 at its own optimal threshold. Passing the cutoff is a statement
  about that model, not a uniform confidence level.
- **MMseqs2 identity scale is version-dependent.** MMseqs2 changed ``pident``
  from a fraction to a percentage. If ``MMSEQS2_MODULE`` resolves to a legacy
  build that still emits fractions, every identity value would be at most 1.0
  and the 50% filter would silently reject every VFDB hit. The signature in
  ``annotation_qc.tsv`` is a ``vfdb`` row with a large ``reported_records``
  and ``retained_records`` of 0; check the module version before concluding
  that a MAG carries no virulence factors.
- **Thresholds are provenance, not defaults to forget.** Every value actually
  used, including ``kofam_acceptance: native_model_cutoffs``, is written to
  ``annotating/annotation_manifest.yaml``. Report the manifest values rather
  than the documented defaults.

References for the thresholds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Aramaki T. *et al.* (2020) KofamKOALA: KEGG Ortholog assignment based on
  profile HMM and adaptive score threshold. *Bioinformatics* 36:2251-2252.
  https://doi.org/10.1093/bioinformatics/btz859
- Mistry J. *et al.* (2021) Pfam: The protein families database in 2021.
  *Nucleic Acids Research* 49:D412-D419. https://doi.org/10.1093/nar/gkaa913
- Feldgarden M. *et al.* (2021) AMRFinderPlus and the Reference Gene Catalog
  facilitate examination of the genomic links among antimicrobial resistance,
  stress response, and virulence. *Scientific Reports* 11:12728.
  https://doi.org/10.1038/s41598-021-91456-0
- Zhang H. *et al.* (2018) dbCAN2: a meta server for automated
  carbohydrate-active enzyme annotation. *Nucleic Acids Research* 46:W95-W101.
  https://doi.org/10.1093/nar/gky418
- Abby S.S. *et al.* / Néron B. *et al.* (2023) MacSyFinder v2: Improved
  modelling and search engine to identify molecular systems in genomes.
  *Peer Community Journal* 3:e28. https://doi.org/10.24072/pcjournal.250
- Camargo A.P. *et al.* (2023) Identification of mobile genetic elements with
  geNomad. *Nature Biotechnology* 42:1303-1312.
  https://doi.org/10.1038/s41587-023-01953-y
- Pearson W.R. (2013) An introduction to sequence similarity ("homology")
  searching. *Current Protocols in Bioinformatics* 42:3.1.1-3.1.8.
  https://doi.org/10.1002/0471250953.bi0301s42
- Rost B. (1999) Twilight zone of protein sequence alignments. *Protein
  Engineering* 12:85-94. https://doi.org/10.1093/protein/12.2.85
- Tian W. and Skolnick J. (2003) How well is enzyme function conserved as a
  function of pairwise sequence identity? *Journal of Molecular Biology*
  333:863-882. https://doi.org/10.1016/j.jmb.2003.08.057
- Steinegger M. and Söding J. (2017) MMseqs2 enables sensitive protein
  sequence searching for the analysis of massive data sets. *Nature
  Biotechnology* 35:1026-1028. https://doi.org/10.1038/nbt.3988

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
     - The complete native AMRFinderPlus record, the detection ``method`` and
       its evidence tier, drug class and subclass, hierarchy node, closest
       reference accession and name, and HMM metadata when the HMM arm
       produced the call.
   * - ``card``
     - The complete native RGI record, the ``Cut_Off`` tier (Perfect or
       Strict), ARO accession, AMR gene family, resistance mechanism, drug
       class, model type, and any SNPs RGI reported.
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
     - Filter ``source == "ncbi_amrfinder"``. The drug class is in
       ``annotation_type`` and the gene symbol in ``annotation_id``.
   * - ``resistance_target``
     - Read ``details.drug_subclass`` on ``ncbi_amrfinder`` rows. Since 2.6
       these values come from AMRFinderPlus directly rather than from a
       mapping table keyed on HMM accession.
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
