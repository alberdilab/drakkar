# Changelog

This project tracks release notes here from this point forward.

## [Unreleased]

### Added

- No unreleased changes yet.

## [2.3.0] - 2026-08-27

### Added

- `drakkar database latest` queries each database source and reports whether the
  release wired into `config.yaml` is still the newest one available, printing
  the install command for every database that is behind. Nothing is downloaded
  and no output directory is touched. It covers the managed databases (`kegg`,
  `cazy`, `pfam`, `vfdb`, `amr`, `foldseek`) plus the GTDB reference data read
  from `GTDB_DB`. Sources are queried in parallel, and one that cannot be
  reached is reported as `unknown` for that database instead of failing the
  command.
- `drakkar database update` installs the newest release of every managed
  database that is behind, in one command instead of one invocation per
  database. It prints its plan and stops unless `--yes` is given, installs each
  release beside the configured one using the same per-database workflow, and
  repoints `config.yaml` at the releases that installed successfully
  (`--no-set-default` opts out). Installs run sequentially so each keeps its own
  Snakemake lock, log, failure report and `database_versions.yaml`; one that
  fails does not stop the others and never repoints `config.yaml`.
- `config.yaml` gained a `DATABASES_DIR` entry naming the base directory that
  holds the database releases Drakkar installs. `--directory` is now optional on
  every `drakkar database <name>` command and defaults to
  `DATABASES_DIR/<database>`, so a release can be installed with just
  `drakkar database pfam --version Pfam38.2`. Passing `--directory` still wins,
  and the command explains what to set when neither is available. Databases
  installed by other tools (GTDB, CheckM2, SingleM, geNomad, eggNOG) keep their
  own paths and are unaffected.

## [2.2.0] - 2026-08-26

### Added

- The Read fates chart in the preprocessing report carries a "Switch to bases"
  button beside its heading, redrawing the stack in gigabases instead of read
  counts. The button is offered only when every fate has base counts, and is
  hidden without JavaScript and when printing.
- The taxonomy section of the report opens with a circular phylogeny: the
  GTDB-Tk placement tree pruned to the catalogue genomes, drawn as a phylogram
  with one ring per genome attribute around it — phylum in the colours the rest
  of the section uses, then genome size, completeness and contamination. Each
  tip carries its figures as a tooltip, rings whose columns the database does
  not hold are left out, and bacteria and archaea get a figure each. The figure
  is inline SVG, so it adds no second plotting bundle to the page.
- The report database gained a `genome_tree` table, ingested from the pruned
  `annotating/bacteria.tree` and `annotating/archaea.tree` written by the
  `gtdbtk_pruned_trees` rule. Both files are optional, and a taxonomy section
  without them renders as before.
- Workflow commands now validate the databases they need before launching
  Snakemake. Every artifact of a managed release is checked for presence and
  content, including the pressed HMM indices, the KEGG hierarchy JSON, the
  KOfam `ko_list` cutoff table and the Pfam EC table, and the error names the
  missing files plus the `drakkar database` command that reinstalls the
  release. Only the databases the requested module and annotation types need
  are checked. `--skip-database-check` launches without the validation.
- Each run records the databases it used in its `drakkar_<run_id>.yaml` run
  metadata, and a later run in the same output directory is stopped when its
  databases differ from the recorded ones and outputs built with the earlier
  release are still present. Because Snakemake profiles rerun on file
  timestamps only, such a change would otherwise mix two database releases in
  one output directory without any trace. Directories written before this
  release fall back to the paths in `annotating/annotation_manifest.yaml`, a
  change with no affected outputs is reported as information only, and
  `--allow-database-change` proceeds deliberately.

### Changed

- The KOfam `ko_list` cutoff table is now checksummed in the
  `database_versions.yaml` written by `drakkar database kegg`, alongside the
  profile database and the hierarchy JSON.
- The report schema version moved to 4 for the new `genome_tree` table. An
  existing `drakkar.db` built by an earlier version is refused with the message
  it already gave for a schema change; rebuild it with `drakkar report --force`.

## [2.1.8] - 2026-08-26

### Changed

- Runtime requests were cut for the annotation and catalogue-indexing rules to
  match observed wall-time usage, which was far below the request. The
  size-scaled coefficient was reduced to a fraction of its previous value:
  `genomad` to 20%, `defensefinder` and `dbcan3` to 15%, `antismash`, `dbcan2`
  and `dbcan_summary` to 10%, and `dbcan`, `dbcan4`, `antismash_regions` and
  `index_catalogue` to 5%. The 10- and 15-minute floors and the per-attempt
  doubling on retries are unchanged, as are all memory requests.

### Added

- The per-rule resource table in the report names each rule's runtime floor —
  the smallest runtime it can request, read from the rule definitions —
  alongside the smallest runtime actually requested and the share of jobs
  sitting at it, so a rule priced by its floor rather than by its size
  coefficient can be told apart at a glance.

### Fixed

- Cluster annotation merging no longer aborts a MAG when a source reports the
  same cluster identifier twice. DefenseFinder does this for systems it lists
  in more than one equally scoring solution, which failed the whole
  `merge_cluster_annotations` job with a duplicate-key error. Records that are
  identical in every field are now collapsed into one, and records that share
  an identifier but differ in content — genuinely distinct clusters the source
  failed to name apart — each get an ordinal suffix (`<sys_id>#1`, `#2`). The
  native record is left untouched in the `details` column, the run log names
  what was collapsed and what was renamed, and the per-MAG QC JSON counts both.

## [2.1.7] - 2026-08-26

### Added

- The preprocessing section reports the estimated microbial fraction (SingleM)
  and Nonpareil coverage in a subsection of their own, with the sequencing
  effort spent against the effort projected for near-complete coverage, and a
  chart putting microbial fraction beside Nonpareil completeness per sample.
  The subsection appears only when `--fraction` or `--nonpareil` was run.
- Mean host reads are summarised beside the mean metagenomic reads, and every
  mean read count carries the same quantity in gigabases underneath.
- The Bins subsection of the cataloging section opens with the totals: how many
  bins the binners produced between them, how many survived as final bins after
  Binette reconciled them, and what share of the first the second is.
- Cataloging reports what Binette did, which the report previously never
  mentioned. Below the bins-per-binner chart, a Sankey diagram traces every bin
  each binner produced to its fate — kept as a final bin, kept but produced
  identically by another binner too, or replaced by a better-scoring bin — with
  the bins Binette assembled itself out of contigs no single binner had grouped
  that way shown as their own stream. The same figures are listed as a table
  underneath. The per-assembly Binette reports already written to
  `cataloging/final/<assembly>.tsv` are the source, so existing output
  directories gain the subsection without being re-run.
- The genome abundance table of the profiling section reports the metagenomic
  reads each sample handed the mapper and the share of them that landed on a
  catalogue genome, beside the mapped read count it already carried. All three
  are summarised as highlights above the table: mapped reads on their own say
  nothing about how well the catalogue represents a sample, and the share is
  what does.
- The job outcomes subsection draws the mix of final job states as a horizontal
  stacked bar across the page, so a run's failing slice is visible without
  being looked for in the table of percentages underneath.
- The per-rule resource table and chart report the 95th percentile of memory
  and runtime used alongside the median. The median says how much of a
  reservation a normal job leaves unused; the 95th percentile is the ceiling
  the rule's heaviest jobs reach, and is how far a request can actually be cut.
  In the chart each bar is the rule's median job and its whisker reaches that
  ceiling.
- Dereplication highlights the share of bins retained, and reports how the
  identity threshold actually acted. dRep's own data tables are read where they
  are present: the number of bins that had a MASH neighbour close enough to be
  compared at all, a histogram of the pairwise identities of those comparisons
  split by whether the pair was collapsed, with the threshold drawn on it, and
  a table counting the pairs in half-percent identity bands from 100% down.
  Reciprocal comparisons are averaged into one pair, as dRep clusters on.

### Changed

- The read fates chart runs the full width of the page, with its legend laid
  out in a single row above the plot instead of down its right side.
- Metagenomic bases are shown in gigabases rather than as a raw base count.
- Dereplication shows its yield as one horizontal stacked bar — retained beside
  collapsed within the same whole — instead of two bars whose ratio the reader
  had to work out.
- The phylum composition chart runs the full width of the page, with its key
  laid out as a single row spread across the width below the plot instead of
  stacked into several rows beside it. It names the eight most abundant phyla,
  down from twelve: past eight a long phylum name no longer fits its share of
  that row, and twelve named phyla also outran the palette, which drew two of
  them in the colour of another.
- Both phylum figures in the taxonomy section draw a phylum in the same colour.
  The colours are assigned once, in abundance order, and looked up by name, so
  the genomes-per-phylum chart — which ranks its bars by genome count instead —
  no longer gives a phylum a different colour from the composition chart below
  it. Phyla the composition chart pools as Other are drawn in its grey.
- The per-rule resource highlights separate the mean runtime of a single job
  from the total runtime of every job, which read as one figure before, and the
  runtime column names itself as per-job.
- The report database gains `assembly_bin`, `assembly_bin_origin`,
  `genome_cluster` and `genome_comparison`, and the schema version is bumped to
  3. A `drakkar.db` written by an earlier build is refused with the usual
  version message and has to be rebuilt with `drakkar report --force`.

### Fixed

- CPU efficiency is no longer 100% for every job. The benchmark asked sacct for
  `CPUTimeRAW`, which is `AllocCPUS x Elapsed` — the CPU time a job *reserved*,
  not the time it burned — so dividing it by that same product reported every
  job as perfectly efficient. `TotalCPU`, the user plus system time the job's
  steps actually consumed, is read instead, and with it the per-job and
  per-rule CPU efficiencies and the run's `used_cpu_sec` become real figures.
  Existing benchmark tables carry the old value and have to be regenerated.
- SingleM's read fraction is read whether it is written as `0.85` or as `85%`;
  the percent form was previously dropped and left the column empty.

## [2.1.6] - 2026-08-25

### Changed

- The HTML report opens with the Drakkar mark above the sidebar title, and the
  title itself now reads *Analysis Report* — the name is already in the image.
  The logo is base64-inlined like everything else on the page, so the report
  stays a single portable file.
- Timestamps are written for a reader rather than a parser. The ingest window
  in the sidebar, the *Ingested* column of the provenance table and the
  *Started*/*Finished* columns of the run table now read
  `25 Aug 2026 at 13:18 UTC` instead of
  `2026-08-25T13:18:19.205575+00:00`; a window that opens and closes within the
  same minute is shown once instead of as a span.
- The sidebar reaches the bottom of the window whatever the section it sits
  beside, instead of ending where its own content ends.

### Fixed

- Resource benchmarking no longer skips most of the workflow. Snakemake prints
  the `rule <name>:` block that carries the job's rule, threads and resources
  only for rules that define no `message:`; for every rule that does — which is
  nearly all of them, `assembly`, `assembly_map`, `metabat2`, `maxbin2`,
  `semibin2`, `comebin` and `binette` among them — it prints only
  `Job <n>: <message>`, so those launches were dropped and never reached the
  benchmark tables or the resource section of `drakkar_report.html`. The rule
  name and its wildcard values are now recovered from the per-job SLURM log
  path the executor reports, and the requested CPUs, memory and runtime — which
  the log does not print in this case — are read from SLURM accounting
  (`ReqCPUS`, `ReqMem`, `TimelimitRaw`).
- Retried jobs are counted again. Snakemake keeps the same internal job id
  across attempts and only the SLURM job id changes, so deduplicating launches
  by internal job id discarded every attempt after the first, undercounting
  relaunches, out-of-memory failures and total allocated CPU time.

## [2.1.5] - 2026-08-25

### Changed

- Report tables are no longer capped at 100 rows. Now that long tables are
  paged in the browser, every row a table has is written to
  `drakkar_report.html` — including the full per-job resource listing, which
  previously showed only the 100 longest-running jobs. The caps that remain
  apply to the figures, which cannot be paged.

### Added

- The Taxonomy section of `drakkar_report.html` now lists the lineage of every
  genome, one row per dereplicated bin and one column per GTDB rank, ordered by
  lineage so relatives sit together. Ranks GTDB-Tk could not place are shown as
  `Unclassified` rather than left blank.
- Every heading in `drakkar_report.html` — each section and each table or
  figure within it — is now followed by a short note in plain language saying
  what the numbers below it are and how to read them, so the page can be handed
  to a reader who did not run the workflow. The notes distinguish statistics
  that look interchangeable but are not, in particular the assembly table's
  read-weighted `Mapping rate %`, pooled across the samples of one assembly,
  against `Mean rate %` in the per-sample table, which averages one sample's
  rate over the assemblies it contributed to.

## [2.1.4] - 2026-08-25

### Added

- The "Runs and resources" section of `drakkar_report.html` now reports the
  computational benchmark of every run that has one. It states how many jobs
  were submitted and what share of them completed, failed, ran out of memory,
  hit their time limit or were relaunched, and it puts the CPUs, memory and
  runtime requested from the scheduler beside those actually used — per
  Snakemake rule, as medians against which a resource profile can be tuned,
  and per job for the longest-running launches. The figures come from the
  benchmark artifacts Drakkar already wrote after each SLURM run
  (`drakkar_<run_id>_resources.yaml` and `benchmark/drakkar_<run_id>.*.tsv`),
  which until now were only summarized in the terminal.
- The report database gained the `run_benchmark`, `benchmark_job` and
  `benchmark_rule` tables holding those figures, so requested-versus-used
  resources can be queried directly. The schema version is now 2; rebuild an
  existing `drakkar.db` with `drakkar report --force`.

### Changed

- The cataloging section is now two tables instead of one wide one. Assemblies
  carry contigs, total length, largest contig, N50, GC and mapping rate; bins
  carry the final, high- and medium-quality counts and mean completeness and
  contamination. The low-quality column is gone — those bins are discarded
  downstream and only the total still counts them. The per-sample mapping
  rates sit with the assemblies, which is what they measure: reads mapped back
  to the assembly, not to the bins recovered from it.
- Every table in the report can be sorted by any column, through the arrows in
  its headings. Numeric columns sort on their underlying values rather than on
  the rendered text, blanks stay at the bottom in both directions, and sorting
  a paged table returns it to the first page. The markup stays plain, so a
  page opened without scripting still shows every row in the renderer's order.

### Fixed

- `drakkar_<run_id>_resources.yaml` matches the same glob as the run metadata
  files and carries the same `run_id`, so the report ingested it as a run of
  its own and overwrote that run's provenance row — losing its version,
  modules, timestamps and status. Run metadata is now recognized by suffix.

## [2.1.3] - 2026-08-25

### Changed

- `drakkar_report.html` now opens as a left sidebar beside the report body.
  The sidebar carries the general information about the report — versions,
  runs, ingest window, and which sections were rendered, unavailable or
  excluded — together with the table of contents, which doubles as the
  navigation: one section is shown at a time instead of the whole report
  being laid out on a single page.
- Tables longer than twenty rows are paged in the browser. Every row is still
  written into the file, so the page remains complete when printed or opened
  without scripting; only twenty are on screen at a time.
- Every table is now preceded by the averages of its relevant numeric columns,
  shown as highlight cards, so per-sample and per-assembly tables state their
  central tendency before the reader scans the rows. The catalogue summary and
  the quantified-gene lengths, previously two-column metric tables, are shown
  as those same cards.

## [2.1.2] - 2026-08-25

### Fixed

- The `drakkar report` HTML tests used `unittest.TestCase.enterContext`, which
  only exists from Python 3.11 onwards, so they errored on the Python 3.10 job.
  They now register their temporary directories with `addCleanup` instead.

## [2.1.1] - 2026-08-25

### Added

- `drakkar report` now renders `drakkar_report.html` alongside `drakkar.db`.
  The page is a single self-contained file: the stylesheet is inlined and the
  Plotly bundle is embedded once, in the first figure, so the report opens
  offline and can be emailed or archived as-is.
- The report is rendered from the database and nothing else, so the source
  `.tsv`/`.tsv.xz` tables are still read exactly once, at ingest. Every figure
  and table comes from a SQL aggregate, so the tens of millions of rows that
  `gene_annotation`, `cluster_annotation` and `gene_expression` can hold are
  never pulled into memory.
- The page opens with a summary header naming the Drakkar version, the report
  schema version, the run identifiers, the ingest timestamps, and which
  sections were rendered, which were unavailable, and which were excluded by
  `--sections`. Sections absent from the database are named on the page rather
  than silently dropped, and a provenance table traces every number back to the
  file it was ingested from.
- New `--html-only` flag re-renders `drakkar_report.html` from an existing
  `drakkar.db` without re-ingesting the source tables. It is mutually
  exclusive with `--db-only`.

### Changed

- `--db-only` still stops after the database, and `--force` still governs
  rebuilding it. The HTML report is a derived artifact and is always
  overwritten.

## [2.1.0] - 2026-08-25

### Added

- New `drakkar report` command that builds `drakkar.db`, a lean SQLite
  projection of whatever tables a Drakkar output directory contains. The
  database is the single source of truth for the forthcoming HTML report, so
  the source tables are read exactly once, at ingest.
- The command screens the output directory and builds only the sections whose
  inputs are present, naming the missing files for anything it skips. Sections
  can also be selected explicitly with `--sections`
  (`preprocessing`, `cataloging`, `dereplication`, `profiling`, `taxonomy`,
  `function`, `expression`, `resources`, or `all`). `--list` reports what an
  output directory can support without building anything.
- `--db-only` skips HTML rendering, `--primary-hits-only` keeps just rank-1
  annotation hits for a smaller database, and `--force` rebuilds an existing
  one.
- The database records its own provenance: `schema_version` stamps the layout
  version and the Drakkar version that wrote it, and `ingest_log` records the
  source file, size, modification time and row count behind every table.

### Changed

- `annotating/genome_taxonomy.tsv` is parsed rather than copied when it reaches
  the report database. The file is a raw concatenation of the GTDB-Tk
  summaries, so its semicolon-packed `classification` string is split into
  seven rank columns and its free-text blob columns are dropped, leaving the
  ranked lineage alongside the placement evidence (`ani`, `af`, `red_value`,
  `msa_percent`, `warnings`).
- The long-form annotation tables are decluttered on the way into the database:
  gene coordinates move into a `gene` table instead of repeating on every hit
  row, annotation labels move into an `annotation_term` dimension table instead
  of repeating on every occurrence, the `details` JSON is dropped, and
  `source=prodigal` rows are omitted from `gene_annotation` because `gene`
  already records every predicted gene. The `.tsv.xz` tables remain the
  lossless archival record.
- Packed fields are normalized: `sample_mapping_rates`, `samples` and
  `coverage_samples` become `assembly_sample` rows, the four `<binner>_bins`
  columns become `assembly_binner` rows, and the wide `counts.tsv` / `bases.tsv`
  matrices are melted into long-form `genome_count` rows.

### Fixed

- The "structure annotation is work in progress" message no longer hardcodes a
  release number, so it reports the running Drakkar version instead of going
  stale at each version bump.

## [2.0.0] - 2026-08-25

### Breaking changes

- `gene_annotations.tsv.xz` is now a lossless long-form hit table instead of a
  one-row-per-gene wide summary. Every row includes explicit MAG and gene
  coordinates, annotation source, method and evidence provenance, hit rank,
  primary-hit status, native scores and alignment coordinates. All qualifying
  hits are retained; source-specific fields are stored as structured JSON in
  `details`. A Prodigal row keeps every gene in the table even when it has no
  functional annotation. Existing gene tables must be regenerated, and 1.x
  consumers must migrate from database-specific columns to the `source`,
  `annotation_id`, and `annotation` fields documented in the
  [annotation table reference](docs/source/annotation_tables.rst).
- `cluster_annotations.tsv.xz` now includes explicit `mag` and stable
  `cluster_id` fields plus method, evidence, and native-record provenance.
  Cluster inputs are schema-validated, and DefenseFinder system identifiers
  are no longer written into the `contig` column.

### Changed

- KOfam searches now retain HMMER's permissively reported candidates so the
  per-model full-sequence or domain bit-score thresholds from `ko_list` are
  authoritative; a global `1e-10` reporting filter no longer discards hits
  before those native cutoffs can be applied.
- KEGG annotation now fails if its native KOfam `ko_list` cutoff table is
  missing or empty instead of silently changing to a global e-value policy.
- VFDB/MMseqs searches now retain query and target lengths and coverage, with
  configurable minimum query and target coverage fractions (both default to
  `0.5`) in addition to identity and e-value filtering.
- Primary hits are ranked with source-specific evidence: KOfam cutoff margin,
  Pfam/AMR bit score, dbCAN coverage, VFDB alignment coverage then identity and
  bit score, or native predictor confidence. `rank_score` and
  `rank_score_type` record the ranking basis.
- Functional annotation writes `annotation_manifest.yaml` and
  `annotation_qc.tsv`, tying result-table checksums to enabled sources,
  filters, configured tool modules and environments, database releases, and
  per-MAG retained/rejected/unmapped record counts.
- Structure annotation with Foldseek/ProstT5 remains work in progress. It is no
  longer advertised or accepted as an annotation target in 2.0.0, and its
  experimental database installer is not exposed by the CLI.
- CAZy gene annotation now consumes dbCAN's coverage-aware HMM output instead
  of filtering raw `hmmscan` hits only by Drakkar's global e-value. Accepted
  calls use dbCAN's native `i-Evalue < 1e-15` and HMM coverage `> 0.35` rules;
  every accepted non-overlapping CAZyme domain is retained.
- The supported Python floor is now 3.10, matching syntax already used by the
  package. Tests run on Python 3.10, 3.11, and 3.12.
- The dbCAN environment pins Python 3.12 and pyhmmer 0.11.4. dbCAN 5.2.5 is
  incompatible with pyhmmer 0.12's string-returning alignment API and can exit
  successfully after its HMM search fails; the CAZy rule now also rejects a
  missing or empty coverage-filtered result explicitly.
- Tag-triggered releases now rerun the complete test suite on Python 3.10,
  3.11, and 3.12 before building or publishing artifacts, and verify that the
  tag matches the package version. Release instructions stage new files and
  wait for main-branch CI before creating the tag.
- Package metadata now exposes project links and classifiers, and installation
  documentation distinguishes the PyPI release from the GitHub development
  version while making the module-based HPC runtime requirement explicit.

### Fixed

- Partial functional reruns now merge only explicitly enabled sources, so
  stale files from a previously enabled gene or cluster annotation cannot leak
  into the new combined tables.
- Gene annotation merging now fails if the Prodigal GFF yields duplicate gene
  identifiers or if any functional hit does not resolve to a Prodigal gene in
  the same MAG, preventing invalid rows with missing coordinates.
- Repeated KO placements in the KEGG hierarchy are collapsed before joining
  hierarchy metadata to KOfam results, preventing one biological hit from
  being duplicated and ranked multiple times in the long-form table. Distinct
  EC associations remain preserved on the single hit row.
- VFDB Set B mapping now reads the VFC-bearing classification block instead of
  misidentifying the terminal organism block as `vf_type`. Corrected mappings
  carry `mapping_schema=drakkar-vfdb-v2`; Drakkar refuses legacy mappings with
  an actionable rebuild instruction. Install a fresh VFDB release with
  `drakkar database vfdb --directory /path/to/vfdb --set-default` before
  rerunning virulence annotation.

## [1.8.22] - 2026-08-20

### Added

- The `profiling_genomes` workflow now writes a third table to
  `profiling_genomes/final/`: `mags.tsv`, with one row per dereplicated genome
  and the columns `magid`, `completeness`, `contamination`, `size`, `contigs`,
  `longest_contig`, `n50`, `gc`, `cluster`, `cluster_members`, and `score`.
  The `magid` values match the genome identifiers used in `counts.tsv` and
  `bases.tsv`. Completeness and contamination come from CheckM2 (or from the
  file passed with `-q/--quality`), the contig statistics are computed from the
  dereplicated genome FASTA files, and the cluster columns come from dRep.
  `mags.tsv` is also included in `drakkar transfer --profile`.

## [1.8.21] - 2026-08-18

### Fixed

- `comebin` now puts its conda environment's `bin` first in `PATH` (and clears
  `PYTHONHOME`) before calling `run_comebin.sh`. COMEBin's launcher invokes a
  bare `python`, and when a conda environment is already active in the
  submitting shell, `conda activate` swaps the rule's environment into the
  *position* of the previously active one instead of prepending it, leaving the
  snakemake module's Python ahead of the environment's Python 3.7. COMEBin then
  ran under the wrong interpreter and failed with
  `ModuleNotFoundError: No module named 'torch'`, even though importing `torch`
  worked when the same environment was activated by hand.

## [1.8.20] - 2026-08-18

### Fixed

- `semibin2` and `comebin` now clear `PYTHONPATH` and set `PYTHONNOUSERSITE=1`
  before running, so the snakemake module's Python packages and `~/.local` can no
  longer shadow the packages inside their conda environments. SemiBin2 was
  crashing with `ModuleNotFoundError: No module named 'narwhals.stable.v2'`
  because an older `narwhals` from outside the environment was imported by the
  environment's `scikit-learn`, which resolves `narwhals >= 2.0.1` correctly
  inside the environment itself.

## [1.8.19] - 2026-08-18

### Fixed

- `semibin2` and `comebin` no longer size their memory and runtime requests from
  the assembly alone. Both binners read every coverage BAM to build the depth
  profiles used for training, so the requests now scale with the assembly size,
  the largest BAM and the number of BAMs, and the floor was raised from 16 GB to
  32 GB. A 186 MB assembly with a 1.8 GB BAM previously requested 16 GB and was
  killed by the SLURM out-of-memory handler while training. Retry scaling is now
  applied after the floor, so a job that was OOM-killed at the floor asks for
  more memory on the next attempt instead of repeating the same request.

## [1.8.18] - 2026-08-17

### Added

- `drakkar environments --list` reports every conda environment deployed in the
  environment directory (`-e/--env_path`, or `ENVIRONMENTS_DIR`) with its size
  and creation date, and classifies it as `in use`, `orphan` (built from a
  superseded definition), `incomplete` (a failed deployment) or `unknown` (a
  directory that is not a Snakemake environment). Environments are matched by
  comparing the `<hash>.yaml` copy Snakemake stores next to each environment
  against the definitions shipped in `workflow/envs/`, so the check does not
  depend on Snakemake's internal hashing.
- `drakkar environments --prune` removes orphaned and incomplete environments,
  along with definition files left without an environment. It is a dry run
  unless `--yes` is given, and never touches directories that are not named
  like a Snakemake environment. `--no-size` skips size computation on slow
  shared filesystems.

### Changed

- The `comebin` environment now carries a post-deployment script
  (`drakkar/workflow/envs/comebin.post-deploy.sh`) that replaces the
  conda-resolved PyTorch with the CUDA 11.3 build (`torch==1.10.2+cu113`, the
  last CUDA wheel published for Python 3.7) and verifies it, so COMEBin detects
  the GPU on the `gpuqueue` nodes. This mirrors what was already done for
  `semibin`. Because Snakemake hashes the post-deploy script together with the
  environment definition, the `comebin` environment is redeployed under a new
  hash the next time `drakkar environments` runs.

## [1.8.17] - 2026-08-17

### Fixed

- `drakkar environments` did nothing and exited with "Nothing to be done (all
  requested files are present and up to date)" without creating any
  environment. The `rule all` of the environments workflow lived in
  `rules/environments.smk`, and Snakemake's `Workflow.include()` restores the
  default target after parsing an included file, so the run had no target at
  all. The target rule is now declared in `workflow/Snakefile` after the
  include, using the `ENVIRONMENTS_TARGET` marker exported by the rules file.
  The command has been a silent no-op since the environments workflow was split
  into its own rules file; environments were only ever created on demand by
  actual workflow runs.

## [1.8.16] - 2026-08-17

### Changed

- `semibin2` now runs in its own Conda environment
  (`drakkar/workflow/envs/semibin.yaml`, SemiBin 2.3.0 on Python 3.11) instead
  of loading the `semibin`, `cuda`, `bedtools` and `hmmer` system modules. A
  post-deployment script (`drakkar/workflow/envs/semibin.post-deploy.sh`)
  replaces the conda-resolved PyTorch with the CUDA 12.4 build
  (`torch==2.5.1+cu124`) and verifies it, so the environment matches the setup
  validated on the A100 nodes. The `SEMIBIN2_MODULE`, `CUDA_MODULE` and
  `BEDTOOLS_MODULE` config keys have been removed, and the environment is
  pre-created by `drakkar environments`.

## [1.8.15] - 2026-08-17

### Changed

- Bumped the Binette module from `binette/1.0.5` to `binette/1.2.1`.
- `comebin` now runs in its own Conda environment
  (`drakkar/workflow/envs/comebin.yaml`, COMEBin 1.0.4 on Python 3.7.12)
  instead of loading the `comebin` and `cuda` system modules. The
  `COMEBIN_MODULE` config key has been removed, and the environment is
  pre-created by `drakkar environments`.

## [1.8.14] - 2026-08-16

### Added

- `-r/--reference` in `drakkar preprocessing` and `drakkar complete`, and the
  `reference_path` column of the sample info file, now accept NCBI genome
  assembly accessions (e.g. `GCF_000001405.40`, `GCA_000001635`) in addition to
  local paths and URLs. The accession is resolved against the NCBI genomes FTP
  site and the matching `*_genomic.fna.gz` file is downloaded into the run's
  reference cache. Unversioned accessions resolve to the latest available
  assembly version.
- Sample info tables can now be comma-separated or semicolon-separated as well
  as tab-separated. The delimiter is detected from the header line rather than
  from the file extension, so a `.csv` sample table no longer fails with a
  misleading `Missing value in column 'sample'` error.

### Changed

- All sample-table reads go through a single `read_input_table()` helper in the
  new `drakkar/input_tables.py`, which also trims whitespace and byte-order
  marks from column names and reports unreadable tables as input errors instead
  of raising a pandas traceback. `drakkar dereplicating`/`profiling` quality
  files (`-q/--quality`) now use the same delimiter detection, replacing their
  separate `csv.Sniffer`-based handling.

## [1.8.13] - 2026-08-15

### Added

- Bin filtering thresholds for cataloging are now forwarded to Binette in full:
  `MAX_CONTAMINATION` (default 10) is passed as `--max_contamination`, and the
  new `MIN_BIN_LENGTH` (default 200000) and `MAX_BIN_LENGTH` (default 10000000)
  config entries are passed as `--min_length` / `--max_length`. Previously only
  `--min_completeness` was applied, so bin contamination was not filtered at all
  during cataloging.
- `--min-completeness`, `--max-contamination`, `--min-bin-length`
  (`--min-length`) and `--max-bin-length` (`--max-length`) arguments for
  `drakkar cataloging` and `drakkar complete`. When omitted, the values in
  `workflow/config.yaml` are used.

### Changed

- Binette upgraded from 1.1.2 to 1.2.1 in `workflow/envs/cataloging.yaml`, which
  is the first release exposing the contamination and bin length filters.
  Binette 1.2 renames the `bin_id` column of `final_bins_quality_reports.tsv` to
  `name` (`binette_bin<n>`) and names bin FASTA files after it; the workflow now
  reads bin ids through `workflow/scripts/bin_report.py`, which accepts both the
  new and the legacy report layouts. Exported bin names
  (`<assembly>_bin_<n>.fa`) are unchanged.

## [1.8.12] - 2026-08-14

### Added

- `rules/genomics.smk`, a host-genomics branch that runs GATK duplicate marking
  and scattered joint germline variant calling on the host alignments produced
  by preprocessing. GATK variant outputs are indexed with `bcftools index -c`
  (CSI) rather than tabix, because `.tbi` shares BAI's ~512 Mbp ceiling. Not
  yet wired into the CLI.
- `reference_faidx` rule, producing `data/references/<reference>.fna.fai`.
  Kept separate from `prepare_reference` so references prepared by earlier runs
  gain the index without re-running `bowtie2-build`.
- `BCFTOOLS_MODULE` config entry.

### Changed

- **Reference mapping in preprocessing now produces CRAM instead of BAM.**
  `preprocessing/bowtie2/<sample>.bam` becomes `<sample>.cram`, and
  `preprocessing/final/<sample>.bam` becomes `<sample>.cram` with a `.crai`
  index alongside it. BAM's companion index, BAI, cannot address positions
  beyond 2^29-1 (~512 Mbp), so references carrying chromosomes above that size
  could not be indexed at all and `samtools_stats` failed with
  `cannot be stored in a bai index. Try using a csi index` /
  `Numerical result out of range`. The CRAM index has no such ceiling. CRAM is
  also roughly half the size of the equivalent BAM.
- `samtools_stats` now emits `<sample>.cram.crai` instead of `<sample>.bam.bai`.
- The failure report table is now written to
  `drakkar_<run_id>_failures.tsv` in the root of the output directory, next to
  `drakkar_<run_id>.yaml` and `drakkar_<run_id>_resources.yaml`, instead of
  `log/drakkar_<run_id>.failures.tsv`.

## [1.8.11] - 2026-08-14

### Added

- Tabular failure report printed automatically when Snakemake stops after
  failures. It lists one row per failed job (rule, target sample/assembly/MAG,
  failed attempts, failure category, and a short detail line), followed by a
  "what to do next" section with per-category recommendations and a verdict
  stating whether relaunching as is, relaunching with larger
  `--time-multiplier`/`--memory-multiplier`, or a manual fix is needed.
  Failure categories are inferred from the SLURM job state, the Snakemake
  error message, and the job logs of the failed rules.
- The same report is written to `log/drakkar_<run_id>.failures.tsv`, recorded
  in the run metadata, and can be reprinted with
  `drakkar logging -o <output_dir> --failures`, which also shows per-job log
  excerpts. Jobs that failed but succeeded on a retry are reported separately
  from unresolved failures.

### Changed

- Workflow failures now exit with the Snakemake exit code and a short pointer
  to the failure report instead of raising an uncaught `CalledProcessError`
  traceback.

## [1.8.10] - 2026-06-20

### Fixed

- KEGG KOfam database build now downloads the archived `ko_list.gz` (the
  archive no longer serves an uncompressed `ko_list`, which caused a 404) and
  decompresses it to `*_ko_list.tsv`.

## [1.8.9] - 2026-06-19

> Historical note: the structure features described below were experimental.
> They remain work in progress and their annotation target and database command
> are not available in Drakkar 2.0.0.

### Added

- New opt-in `structure` (alias `foldseek`) annotation type that structurally
  annotates homology-orphan genes with Foldseek/ProstT5 against the
  AlphaFold/Swiss-Prot database. ProstT5 runs inside the `foldseek` module
  (CPU, no GPU), gated on genes that KEGG/Pfam/CAZy did not annotate. Structural
  KO/EC/Pfam hits fill the existing columns only where sequence homology was
  empty, and a new `evidence` column (`sequence`/`structure`) records provenance.
- Added `FOLDSEEK_MODULE`, `FOLDSEEK_DB`, `PROSTT5_MODEL`, and `FOLDSEEK_MAP_DB`
  configuration entries.
- Added `build_foldseek_function_map.py` to generate `FOLDSEEK_MAP_DB` (UniProt
  accession -> KO/EC/Pfam) from UniProt's `uniprot_sprot.dat.gz`.
- New `drakkar database foldseek` command that installs the whole Foldseek bundle
  in one step: the AlphaFold/Swiss-Prot structure DB and ProstT5 weights (via
  `foldseek databases`) plus the UniProt function map. `--set-default` updates
  `FOLDSEEK_DB`, `PROSTT5_MODEL`, and `FOLDSEEK_MAP_DB` together. Like `vfdb`, the
  `--version` is optional and defaults to the UTC download date. Managed-database
  entries can now declare multiple config keys via `config_targets`.

### Changed

- KEGG and AMR annotations now use their databases' native per-model cutoffs
  instead of the flat `annotation_evalue`, mirroring the `--cut_ga` approach
  already used for Pfam:
  - KEGG: applies KOfam's per-KO adaptive bit-score thresholds from the `ko_list`
    file (full- or domain-score per `score_type`), falling back to the e-value
    cutoff only for KOs that have no curated threshold. The `ko_list` is now
    downloaded alongside the KOfam profiles by `drakkar database kegg`.
  - AMR: `hmmscan` now runs with `--cut_tc` (NCBIfam-AMRFinder trusted cutoffs,
    as used by AMRFinderPlus); the flat e-value filter is removed.

## [1.8.8] - 2026-06-18

### Changed

- Updated `SEMIBIN2_MODULE` from `semibin/2.1.0` to `semibin/2.3.0`, which includes a CUDA-enabled PyTorch build.
- Updated `CUDA_MODULE` from `cuda/12.9.1` to `cuda/12.8` to satisfy `semibin/2.3.0`'s prerequisite.
- Fixed module load order in `semibin2` and `comebin` rules: `CUDA_MODULE` is now loaded before the tool module, as required by the module system.

## [1.8.7] - 2026-06-18

### Changed

- Pfam annotation now uses HMMER's per-family curated gathering thresholds (`hmmscan --cut_ga`) instead of a flat e-value cutoff, matching how Pfam/InterPro/UniProt assign annotations. The flat `ANNOTATION_EVALUE` filter is no longer applied to Pfam hits in `merge_gene_annotations.py` (it would otherwise discard GA-passing hits whose full-sequence e-value is weaker than the threshold); the e-value is retained only to select the single best hit per gene. This is more sensitive for short/divergent families and may change existing Pfam results. Other HMM databases (KEGG, CAZy, AMR) are unaffected.

### Fixed

- The `genomad` rule failed with a numpy version error (module ships numpy 2.2.6, but 2.4.6 was imported) because Snakemake injects `PYTHONPATH` pointing at the drakkar conda env, which shadowed the genomad module's own `site-packages`. The rule now runs `unset PYTHONPATH` and `export PYTHONNOUSERSITE=1` after `module load` so the module's pinned numpy is used.

## [1.8.6] - 2026-06-13

### Fixed

- `query_sacct_for_jobs` now batches job IDs into chunks of 500 before calling `sacct`, preventing an `OSError: [Errno 7] Argument list too long` crash on runs with a large number of jobs.

### Changed

- Launch metadata now records `drakkar_version` so the exact package version used for a run is captured in the metadata file.

## [1.8.5] - 2026-06-13

### Fixed

- Virulence factor annotations were silently dropped because `mmseqs easy-search` outputs `fident` (fractional identity, 0–1) by default while `parse_vfdb` compared against a percentage threshold (`identity >= 50.0`). Every hit evaluated `0.99 >= 50.0 → False`, filtering the entire result set. Fixed by adding `--format-output query,target,pident,...` to the `vfdb` rule so the identity column is in percentage (0–100), consistent with the threshold. Added a regression test to catch any future reversion.
- Doubled the default runtime for `gtdbtk` to 360 minutes to accommodate larger datasets.

## [1.8.4] - 2026-06-09

### Fixed

- `nonpareil` rule now uses `shadow: "minimal"` so that the decompressed temporary FASTQ and nonpareil intermediate files (`.npl`, `.npa`, `.npc`) are written into Snakemake's per-job scratch directory and cleaned up automatically — even when the job is killed with SIGKILL, which bypasses shell `trap` handlers.

## [1.8.3] - 2026-06-08

### Changed

- `semibin2` and `comebin` rules now load `CUDA_MODULE` (configurable in `config.yaml`, default `cuda/12.9.1`) so that GPU hardware is properly detected at runtime.
- Minimum memory for `semibin2` and `comebin` raised from 8 GB to 16 GB.
- Minimum runtime for `semibin2` and `comebin` raised from 30 to 120 minutes; per-MB scaling increased 4× (from `size/2` to `size*2` minutes).

## [1.8.2] - 2026-06-05

### Fixed

- `MIN_COMPLETENESS` and `MAX_CONTAMINATION` were cast to `float` but Binette expects integers; corrected to `int`.

## [1.8.1] - 2026-06-03

### Added

- Added `--annotation-type clusters` shortcut that expands to `dbcan`, `antismash`, and `mobile`, producing only `cluster_annotations.tsv.xz` without triggering gene annotation. `defense` is intentionally excluded so that `RUN_GENE_ANNOTATIONS` remains false; users who want defense systems can combine `clusters,defense` explicitly.

## [1.8.0] - 2026-06-01

### Changed

- `semibin2` and `comebin` rules now submit to the `gpuqueue` SLURM partition and request one GPU via `--gres=gpu:1`, enabling GPU-accelerated binning for both tools.
- Fixed memory and runtime scaling for `semibin2` and `comebin`: resources are now derived from assembly file size alone instead of total input size (assembly + BAM files), which previously caused grossly inflated RAM and walltime requests.

## [1.7.10] - 2026-05-31

### Changed

- `--scratch_dir` in GTDB-Tk `classify_wf` is now opt-in via `--gtdb-scratch` on `drakkar annotating` and `drakkar complete`. By default pplacer loads the reference tree into memory (faster); pass `--gtdb-scratch` to write intermediates to disk and reduce peak RAM.

## [1.7.9] - 2026-05-31

### Changed

- Added `--scratch_dir` to the GTDB-Tk `classify_wf` call to reduce pplacer memory usage by writing intermediate files to disk.

## [1.7.8] - 2026-05-30

### Added

- Added `MIN_COMPLETENESS` (default 70) and `MAX_CONTAMINATION` (default 10) parameters to `workflow/config.yaml`. These are forwarded as `--min_completeness` to Binette in the cataloging workflow, and as `--completeness` / `--contamination` to dRep in the dereplicating and profiling-genomes workflows. Both values can be overridden in the project config file.

## [1.7.7] - 2026-05-26

### Added

- Added `--annotation-evalue` and `--annotation-identity` to `drakkar annotating` and `drakkar complete`, forwarding the thresholds to `merge_gene_annotations.py` for configurable merged annotation filtering. The default e-value threshold is `1e-10`, and the default identity threshold is `50`.

## [1.7.6] - 2026-05-21

### Fixed

- Ensured GTDB-Tk tree pruning runs in a dedicated Biopython conda environment and calls that environment's Python, preventing `ModuleNotFoundError: No module named 'Bio'` when the launcher environment lacks Biopython.

## [1.7.5] - 2026-05-20

### Added

- Added `--sanitize` flag to `drakkar preprocessing` and `drakkar complete` to run `seqkit sana` on each paired FASTQ file before quality filtering, followed by `seqkit pair` to re-synchronise the two mates. The sanitized reads feed directly into fastp, and a `seqkit_sana_reads` column (paired read count after sanitization) is appended to `preprocessing.tsv`; the column is `NA` when `--sanitize` is not used.
- Added `SEQKIT_MODULE: "seqkit/2.10.0"` to `workflow/config.yaml`.

## [1.7.4] - 2026-05-20

### Changed

- Replaced deprecated `--skip_ani_screen` flag with `--place_species` in the GTDB-Tk classify workflow, required by GTDB-Tk v2.7.0+.

## [1.7.3] - 2026-05-20

### Added

- Added `--slurm-qos` flag to all workflow subcommands for passing a SLURM Quality of Service string to `sbatch --qos`.

### Fixed

- Added `module purge` before every `module load` in all workflow rules to prevent module conflicts when a module version is already loaded in the environment.

## [1.7.2] - 2026-05-19

### Fixed

- `dbcan4` rule now passes `--input_raw_data` (the MAG FASTA file) and `--mode prok` to `run_dbcan substrate_prediction`, resolving a missing required input error.

## [1.7.1] - 2026-05-19

### Added

- Added a dedicated **Snakemake and SLURM management** section to the documentation covering resource caps (`SNAKEMAKE_MAX_GB`, `SNAKEMAKE_MAX_TIME`, `MEMORY_MULTIPLIER`, `TIME_MULTIPLIER` in `config.yaml`), resource multiplier flags (`--memory-multiplier`, `--time-multiplier`), all `--snakemake-*` override flags, all `--slurm-*` override flags, and SLURM benchmarking outputs.
- Documented `--overwrite` and `--skip-benchmark` in every workflow module's option list.
- Documented `benchmark/` and `drakkar_<run_id>_resources.yaml` in the outputs reference.
- Documented `drakkar update --skip-deps`.
- Expanded the CLI reference (`api.rst`) with the full set of common, Snakemake override, and SLURM override options.

## [1.7.0] - 2026-05-19

### Added

- Added `--snakemake-latency-wait`, `--snakemake-jobs`, `--snakemake-cores`, `--snakemake-executor`, `--snakemake-retries`, `--snakemake-rerun-incomplete`, and `--snakemake-keep-going` flags to all workflow subcommands for overriding Snakemake profile settings from the command line.
- Added `--slurm-partition`, `--slurm-account`, `--slurm-constraint`, `--slurm-nodes`, `--slurm-nodelist`, and `--slurm-extra` flags to all workflow subcommands for overriding SLURM resource defaults without editing profile files.
- Both override groups are surfaced as dedicated **Snakemake Overrides** and **SLURM Overrides** sections in `--help` for every workflow subcommand.

## [1.6.10] - 2026-05-15

### Added

- `drakkar status` now prompts for a run index when multiple runs are available and no run is specified.

### Fixed

- Final gene and cluster annotation tables now keep only the first per-genome TSV header when merging genome-specific annotation files.

## [1.6.9] - 2026-05-15

### Changed

- GTDB-Tk taxonomy annotation now resolves the configured GTDB database at rule execution time, reports the selected config key and path, and uses the configured GTDB-Tk module for classification.

### Fixed

- ENA/SRA accession FASTQ downloads now validate cached and downloaded files against ENA byte counts when available, retry truncated files, and reject clearly imbalanced paired FASTQ sizes before Snakemake starts when ENA byte counts are unavailable.

## [1.6.8] - 2026-05-10

### Added

- Added `drakkar status` to show latest or selected run progress by rule and sample, with optional run-directory/YAML selection, rule-only/sample-only views, and `--complete` helper-rule detail.

## [1.6.7] - 2026-05-10

### Changed

- Increased cataloging runtime scaling for SemiBin2 and ComeBin so their requested time is based on full assembly size instead of half assembly size.

### Fixed

- Fixed SLURM benchmark capture for current Snakemake executor logs such as `Job <id> has been submitted with SLURM jobid <id>`, and normalized embedded `sbatch` responses before querying `sacct`.

## [1.6.6] - 2026-05-10

### Added

- Added stricter input-transfer preflight checks before Snakemake starts: URL downloads now retry up to five times, empty cached/downloaded files are rejected, SFTP URLs are supported via `curl`, and missing or empty input files/directories are reported together with a clear stop message.

### Changed

- Split the large CLI and utility modules into focused implementation modules while keeping the existing `drakkar.cli` and `drakkar.utils` import surfaces compatible.

## [1.6.5] - 2026-05-09

### Added

- Added `comebin` as a cataloging binner and added `-b/--binners` to `drakkar cataloging` and `drakkar complete` to select comma-separated binners from `metabat`, `maxbin`, `semibin`, and `comebin`; the default remains all binners.

## [1.6.4] - 2026-05-07

### Fixed

- When no assembly column is present in the sample info file (or all assembly values are empty) and no `--mode` flag is given, DRAKKAR now defaults to individual assemblies instead of failing with an empty `assembly_to_samples.json`.

## [1.6.3] - 2026-05-07

### Changed

- Cataloging now prioritises preprocessed reads when resolving input paths from a sample info file: explicit `preprocessedreads1`/`preprocessedreads2` columns are used first, then reads found in `preprocessing/final/` from a prior preprocessing run, and raw `rawreads1`/`rawreads2` or `accession` paths are only used as a last resort. This allows the same input file to carry coverage/assembly grouping metadata while still assembling and mapping against quality-filtered reads.

## [1.6.2] - 2026-05-05

### Added

- Added pruned GTDB-Tk output trees for taxonomy annotation so ``annotating/bacteria.tree`` and optional ``annotating/archaea.tree`` keep only the input MAG genomes and exclude GTDB reference tips.

### Changed

- ENA/SRA and URL file downloads now retry up to three times with exponential backoff (5 s, 15 s, 45 s) before giving up on a single file.
- All file-downloading loops (samples, references, bins, MAGs, preprocessed reads, transcriptomes) now continue past individual download failures and report the full list of failed files at the end instead of stopping on the first error.

## [1.6.1] - 2026-05-04

### Changed

- Made SLURM resource benchmarking explicitly default-on across workflow commands, added a shared ``--skip-benchmark`` opt-out flag, and ensured empty benchmark tables are still written when no submitted jobs are detected so the benchmark output is visible by default.
- Broadened Snakemake SLURM submission parsing to recognize additional submission log styles such as ``Submitted batch job ...`` when building run benchmarks.
- Allowed sample tables for preprocessing, complete, profiling, and expressing to use an ``accession`` column with paired-end ENA/SRA run accessions, which DRAKKAR now resolves and downloads into the existing forward/reverse read caches automatically.
- Styled the screen-session startup hint so the inline ``screen -S mysession`` command is visually marked as code in Rich output and still remains clear in plain-text output.

## [1.6.0] - 2026-05-04

### Added

- Added a first SLURM-backed benchmark layer that parses submitted non-local Snakemake jobs, infers relaunch attempts, joins them to `sacct`, writes per-run benchmark tables under `benchmark/`, emits a root-level `drakkar_<run_id>_resources.yaml` summary for every workflow run, and exposes a compact resource-efficiency summary through `drakkar logging`.
- Added a richer root-level `profiling_genomes.tsv` summary with input reads, input bases, mapped reads, mapping percentage, and mapped bases per sample, and added a new root-level `dereplicating.tsv` summary for profiling and dereplicating runs with input/output bin counts plus average completeness and contamination when quality metadata are available.

## [1.5.4] - 2026-05-04

### Changed

- Animated the DRAKKAR banner sequence in interactive terminals so the ship appears first, the logo follows after a short pause, and the remaining intro or help content appears after a second pause.
- Allowed ``-B/--bins_file`` in profiling, dereplicating, annotating, inspecting, and expressing to point to a remote manifest URL, downloading the manifest and any referenced remote MAG/bin FASTA files into the local cache.
- Made ``--quality`` genome-name matching more flexible so MAG/bin names are accepted with or without FASTA extensions and still map back to the canonical genome filenames used internally.
- Added ``drakkar update --skip-deps`` to let users skip reinstalling Python dependencies during an update when they only want to refresh the package itself.

## [1.5.3] - 2026-05-03

### Changed

- Styled the update-success ASCII banner with the same ship color as the main DRAKKAR banner and rendered the ship version badge with a separate accent color when Rich output is available.

### Fixed

- Made the cataloging Binette checkpoint tolerate empty ``temporary_files/diamond_result.tsv`` or ``diamond_result.tsv.gz`` files by exporting a header-only ``final_bins_quality_reports.tsv`` instead of failing the workflow.

## [1.5.2] - 2026-05-03

### Changed

- Stopped ``drakkar logging`` from dumping the Snakemake log by default, and replaced that section with a short guide showing how to inspect summaries, excerpts, paths, full logs, and available runs.
- Added an explicit ``--excerpt`` mode to ``drakkar logging`` so failure excerpts or fallback tails are shown only on request.

### Fixed

- Ensured cataloging still writes ``cataloging/final/all_bin_metadata.csv`` when Binette yields zero bins, by materializing an empty metadata table with the expected columns.

## [1.5.1] - 2026-05-03

### Added

- Added a plain-help fallback note explaining when Rich-styled help is unavailable because the ``rich`` dependency is missing.

### Changed

- Rendered the DRAKKAR ship and logo above CLI help output, consolidated the top-level workflow table into ``Data Generation and Analysis``, and expanded workflow descriptions in the top-level help menu.

### Fixed

- Changed ``drakkar update`` to reinstall package dependencies instead of skipping them, so upgraded environments can pick up required runtime packages such as ``rich``.
- Normalized release-changelog generation so published versions always remain separated by a blank line.

## [1.5.0] - 2026-05-03

### Added

- Added parsed execution summaries to ``drakkar logging``, including planned jobs, observed rule executions, workflow progress, and detected error types, plus a ``--summary`` mode for summary-only output.
- Added a version-aware success banner at the end of ``drakkar update`` showing the version that was just installed.

### Changed

- Moved DRAKKAR ASCII art assets into ``drakkar/ascii.py`` so display helpers in ``utils.py`` only handle formatting and printing.
- Reworked ``drakkar --help`` and subcommand help into grouped Rich layouts with workflow categories, command families, option sections, and concrete examples.

## [1.4.0] - 2026-05-03

### Added

- Accepted ``assembly`` as the preferred sample-table column for cataloging assembly groups, while keeping legacy ``coassembly`` support.
- Added ``drakkar logging`` plus persistent per-run Snakemake logs under ``log/drakkar_<run_id>.snakemake.log`` to help troubleshoot failed or locked workflow directories.

## [1.3.5] - 2026-05-02

### Added

- Added `--memory-multiplier` and `--time-multiplier` to Snakemake-running Drakkar commands.
- Added `SNAKEMAKE_MAX_TIME` to cap runtime requests, defaulting to 20160 minutes (14 days).

### Changed

- Routed explicit workflow memory and runtime resources through capped scaling helpers, and scaled default Snakemake resources for rules without explicit resources.

## [1.3.4] - 2026-05-02

### Fixed

- Changed the SemiBin2 cataloging table to read reclustered `output_bins/*.fa` files instead of pre-reclustering `contig_bins.tsv`.
- Counted MetaBAT2, MaxBin2, and SemiBin2 bins in `cataloging.tsv` from Binette input quality reports instead of contig-to-bin table rows.

## [1.3.3] - 2026-05-01

### Changed

- Rendered argparse help with Rich panels and tables for command and subcommand help.

### Fixed

- Forced Drakkar CLI truecolor output for normal terminal streams so inherited non-color environment settings do not hide the Rich theme.
- Merged SingleM microbial fraction rows whose sample names include read-mate suffixes such as `_1` back into the main preprocessing summary rows.

## [1.3.2] - 2026-05-01

### Changed

- Moved the startup version badge to the ship banner, centered the intro against the DRAKKAR logo, and enriched CLI status/help styling through the shared Rich output layer.

## [1.3.1] - 2026-05-01

### Changed

- Styled the startup banner with a logo-integrated version badge, centered intro frame, and Rich colors inspired by Viking sea, iron, and gold tones.

## [1.3.0] - 2026-05-01

### Added

- Cataloging now runs QUAST assembly statistics and samtools flagstat mapping summaries, then writes a root-level `cataloging.tsv` summary table.

### Changed

- Host/metagenomic read and base count sidecar files are now temporary preprocessing intermediates; the retained summary is `preprocessing.tsv`.

## [1.2.3] - 2026-05-01

### Added

- All SingleM `pipe` and microbial fraction steps now use the configured `SINGLEM_DB` metapackage.

## [1.2.2] - 2026-04-30

### Added

- Added `--fraction` to `drakkar preprocessing` so standalone preprocessing can run SingleM and SingleM microbial fraction outputs.
- Added `--nonpareil` to `drakkar preprocessing` and `drakkar complete` to estimate Nonpareil coverage/diversity from preprocessed metagenomic reads.
- Expanded `preprocessing.tsv` with fastp, host/metagenomic, SingleM fraction, and Nonpareil summary columns.

## [1.2.1] - 2026-04-30

### Fixed

- Allowed preprocessing `-r/--reference` and `-x/--reference-index` to accept `http`, `https`, and `ftp` URLs.

## [1.2.0] - 2026-04-30

### Added

- Added `-x/--reference-index` for preprocessing and complete workflows to use tarballs containing a host FASTA plus prebuilt Bowtie2 index files, including tarballs listed in sample table `reference_path` values.

### Fixed

- Added a preflight check for current/output directory access before writing Drakkar run metadata, so protected directories produce a clear CLI error instead of a Python `PermissionError`.

## [1.1.3] - 2026-04-27

### Added

- Added a `.gitignore` for Python caches, test artifacts, build outputs, documentation builds, Snakemake state, and local OS/editor files.
- Added `--gtdb-version` to `drakkar annotating` and `drakkar complete` to select numbered `GTDB_DB_<version>` config entries for taxonomy annotation.

### Fixed

- Capped dynamic Snakemake `mem_mb` resource requests using configurable `SNAKEMAKE_MAX_GB`, defaulting to 1024 GB.

## [1.1.2] - 2026-04-25

### Added

- Added Rich-backed CLI output, prompts, help rendering, and `drakkar --version` for the 1.1.2 release.

## [1.1.1] - 2026-04-24

### Added

- Added `drakkar config --view` and `drakkar config --edit` to inspect or modify the installed `workflow/config.yaml` from the CLI.

### Changed

- Added `--overwrite` to output-writing workflows so a locked output directory can be deleted and recreated after a broken Snakemake session, with an interactive prompt when the flag is not provided.
- Changed `drakkar database amr` to download versioned NCBI AMRFinder releases such as `2025-07-16.1` instead of the floating `latest` endpoint.
- Changed `drakkar database kegg` to download KOfam profiles from versioned KEGG monthly archive directories using `YYYY-MM-DD` archive dates, then merge and `hmmpress` the extracted HMMs.
- Changed the default database download/preparation runtime to `120` minutes for the KEGG archive rule, and added `--download-runtime` to `drakkar database` so this limit can be overridden from the CLI.
- Changed managed annotation database config entries to point to the installed release directory, with the workflow resolving the expected internal files automatically.
- Changed `drakkar database pfam` to download `Pfam-A.hmm.gz` from versioned Pfam release directories such as `Pfam37.4`, then unzip and `hmmpress` the requested release.
- Changed `drakkar database vfdb` so `--version` can be omitted and then defaults to the UTC download date, which is used as the release folder and logged version.
- Fixed `drakkar database cazy` to download the requested upstream dbCAN release from the versioned `Databases/<version>/` endpoint instead of silently saving an HTML landing page.

## [1.1.0] - 2026-04-24

### Added

- Added a new `drakkar database <name>` maintenance command for installing one managed annotation database release at a time (`kegg`, `cazy`, `pfam`, `vfdb`, `amr`).
- Added per-run database version logging in `database_versions.yaml`, including source URLs, requested versions, release directories, timestamps, and SHA256 checksums of installed assets.

### Changed

- Managed database installs now use explicit `--directory` and `--version` arguments, and `--set-default` updates the concrete path stored in `config.yaml` only after a successful installation.
- Restored direct-path config usage for non-managed and externally sourced databases such as `GTDB_DB`.

### Documentation

- Documented the new `database` workflow and its version logging behavior in the README and usage guide.

## [1.0.1] - 2026-04-21

### Changed

- Added release automation scaffolding, including this changelog, package-version wiring, and release helper tooling.
