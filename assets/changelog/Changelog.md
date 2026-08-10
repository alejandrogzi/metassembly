<p align="center">
  <p align="center">
    <img width=100 align="center" src="../figures/logo.png" >
  </p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="../figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="../figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="../figures/hillerlab-light.png"
    >
  </picture>
</p>

  <span>
    <h1 align="center">
        xasm
    </h1>
  </span>

  <span>
    <h2 align="center">
        CHANGELOG
    </h2>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/xasm" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/xasm?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://nbisweden.github.io/workshop-RNAseq/2011/lab_assembly.html">metassembly</a> .
        <a href="https://github.com/hillerlab/xasm/blob/main/assets/pipeline/xasm.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

--- 

All notable changes to xasm will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.1.1] - 2026-08-10

### Added

- ORF calling with XORF, gated by `xorf_call_orfs` (default `true`). The pinned `xorf` submodule is wired in as a gitlink and its ORF subworkflow is called through a new `XORF_RUN` wrapper (`subworkflows/xorf/main.nf`) that adapts channel shapes and prepares the protein database: a custom `xorf_custom_database` (`.dmnd`/`.dmnd.gz` used as-is, `.fa`/`.fasta` appended to SwissProt and reindexed with diamond) or the default zenodo download. The full ORF chain runs with per-step parameters (`xorf_chunk_size`, `xorf_skip_netstart`, `xorf_predict_min_score_max_predictions`, `xorf_predict_max_predictions`, `xorf_predict_threshold`, `xorf_predict_keep_raw`, `xorf_do_polishing`, `xorf_rename_*`, `xorf_samba_weights`/`xorf_samba_local_weights`), publishing under `09_polish/xorf/{00_concat,01_renamed,02_merged,03_duplicates,04_results}`; the `rename_predictions.py` helper is vendored into `bin/` and `assets/scripts/xorf.sh` manages the submodule.
- Two-pass polishing now operates on ORF-validated transcripts: when `do_twopass_polish` is set, `xorf_call_orfs` is required (validated at startup) and the twopass HQ input is XORF's ORF predictions rather than the raw first-pass HQ, so intron reclassification and retention rescuing run against transcripts with ORF evidence.
- Nonsense-mediated decay filtering: new `ISOTOOLS_NMD` module runs `iso-nmd` (premature-termination-codon detection) on the final HQ set; NMD candidates publish to `10_final/nmd` and the published final transcript set is the NMD-passing reads.

### Changed

- Output layout: single-exon scraps are published under `10_final/scraps` (previously `09_polish`) and fusion calls under `10_final/fusions` (previously `09_polish/fusions`), grouping every final transcript category under `10_final`.
- XORF input transcripts are named `<prefix>_flnc`; the deterministic id flows through XORF outputs and two-pass artifact names (e.g. `09_polish/twopass/classify/<prefix>_flnc@<prefix>_flnc.reference_introns.tsv`).
- Repository layout: the end-to-end fixture moved from `test_data/` to `assets/test/test_data/` and the test harness from `tests/` to `assets/ci/`; all `nextflow.config` profiles, CI workflows, `run.sh`/`verify.py`, the fixture builder, and docs were updated accordingly. The e2e golden suite now also asserts the XORF chain (fixture protein database, ORF task counts, `09_polish/xorf/` artifacts).
- The `INTRONIC` process label was raised from `process_medium` to `process_high`.
- Version bumped to `0.1.1` in the pipeline manifest.

### Removed

- The `publish.yml` GitHub Actions workflow (automatic release publishing) — no longer needed.

---

## [0.1.0] - 2026-08-07

### Added

- ruSTAR as an alternative aligner (`aligner = "ruSTAR"`): a Rust reimplementation of STAR with the same CLI flags and STAR-compatible log output, reading both FASTQ and CBQ natively. Because the aligner image ships no `samtools`, BAM indexing is delegated to the new `SAMTOOLS_INDEX` module. The two-pass scheme and the pre-curated-junction shortcut apply as in the STAR path.
- CBQ (columnar BINSEQ) read support. Native `.cbq` files in `input_dir` are detected automatically; `bqtools_encode_fastqs` encodes FASTQ inputs before QC (QC then runs via `bqc` instead of `fastp`), and `bqtools_encode_before_alignment` keeps the `fastp`+deacon path and encodes to `.cbq` only for the aligner (the two are mutually exclusive, configurable in `params.json` and `main.nf` help).
- `BQTOOLS_ENCODE` module to collapse a paired FASTQ set into a single `.cbq`.
- `BQC` module: CBQ-native all-in-one QC (adapter, trimming, filtering) with a structured JSON report and `bqc_*` parameters. Paired counts are doubled so `fastp_min_trimmed_reads` means the same thing on both paths.
- deacon upgraded to 0.16.0, which emits `.cbq` on the CBQ decontamination path.
- `test-rustar` end-to-end profile exercising the whole CBQ path (`bqtools_encode_fastqs` → `bqc` → deacon → ruSTAR → `samtools` index → cleanup) on the existing fixture, wired into the nightly matrix and CI profile validation, with a smoke checker that asserts the CBQ processes run and the STAR/FASTP path does not.

### Changed

- New `aligner` parameter (`STAR` | `ruSTAR`, default `STAR`) separating alignment configuration; version bumped to `0.1.0` in the pipeline manifest.
- Coverage tracks are refused up front under ruSTAR (`rustar-aligner` rejects `--outWigStrand Unstranded`), so `aligner = ruSTAR` requires `star_make_coverage = false`.
- Alignment under ruSTAR rebuilds the `(meta, bam, bai)` tuple via `SAMTOOLS_INDEX` instead of treating STAR-style inline indexing as the norm.

### Fixed

- Aletsch per-chromosome cleanup race: with `assembly_by_chr = true` the shared full BAM is read by every per-chromosome local-assembly task, and the `!aletsch_keep_bam && !star_make_coverage` cleanup deleted it after whichever chromosome finished first. The in-task deletion was removed and `REMOVE_BAMS` now waits for the last local assembly before deleting, so `--star_make_coverage false` is safe without `--aletsch_keep_bam`.
- `REMOVE_BAMS` referenced a `biocontainers/bash` image that does not exist on the configured registry; it now uses a valid ubuntu base.
- CHROMSIZE stub failed in `-stub-run`: the stub fell back to the always-empty `meta.id` and created a file where the real `chromsize -o` creates a directory. It now mirrors the script (`genome.baseName`, `mkdir -p`).

---

## [0.0.20] - 2026-08-04

### Added

- Per-chromosome assembly (`assembly_by_chr`, enabled by default). Aligned BAMs are now split by reference sequence with the new `BAMSPLIT_CHROM` module, local assembly and metassembly run per chromosome, and the results are merged into a single metassembly GTF with globally unique transcript IDs. The `assembly_exclude_chromosomes` parameter lets you keep chromosomes in the alignment while omitting them from the assembly fan-out.
- StringTie3 as an alternative local assembler (`assembler = "stringtie3"`), including nascent RNA options (`stringtie3_enable_nascent_assembly`, `stringtie3_include_nascent_rna`).
- TransMeta as an alternative meta-assembler (`metassembler = "transmeta"`), with annotation-guided merging controlled by `transmeta_use_annotation`.
- Two-pass polishing (`do_twopass_polish`). Retention discards from the first pass are reclassified against the first-pass intronIC evidence, UTR retentions are ignored in the second pass, and rescued transcripts are merged back into the clean HQ set. Output lands under `09_polish/twopass/`.
- `star_twopass_junctions_file` parameter to provide a pre-curated junction file (STAR sjdb 4-column format), which skips the STAR first pass and the junction-merging step entirely.
- GFF support: annotations in `.gff` format are now accepted alongside GTF and BED.
- Aletsch insert-size fallbacks (`aletsch_fallback_insert_size`, `aletsch_fallback_insert_std`) that replace a zero profiled insert-size mean when compact per-chromosome inputs fall below Aletsch's internal sample minimum.
- `--intron-track` output for the iso-classify process, emitted during the 0.0.19b milestone.
- End-to-end test suite: a committed 256 KB synthetic fixture (`test_data/e2e/`) with reads, genome, annotation, splice scores, and repeats; `tests/e2e/run.sh` and `tests/e2e/verify.py` with golden counts; new `test`, `test-sb`, and `test-tm` profiles (the old small fixture is now `test-polish`); and a nightly end-to-end workflow (`e2e-nightly.yml`) alongside the expanded CI smoke test.
- Pinned container image digests for the beaver, bed2gtf, chromsize, genepred, isox-py, xloci, and bamsplit modules to make runs reproducible.
- New `REMOVE_BAMS` and `RENAME_BAM` housekeeping modules, and a README logo via `.gitattributes`.

### Changed

- Version bumped to `0.0.20` in the pipeline manifest, which now also lists the authors.
- Beaver and TransMeta outputs are prefixed per chromosome (`<prefix>_<chr>` when `assembly_by_chr` is on) so per-chromosome results can be traced back to their source sample.
- Aletsch now derives its library type from `single_end`/`paired_end` metadata rather than the strandedness field, and supports a per-chromosome `-l` flag. Transcript counts are computed with `awk` instead of `grep -w`.
- The polished GTF/BED conversion path now also handles GFF input in `prepare_indexes` and the `--from polish` checkpoint.
- STAR alignment flows straight to the second pass with a pre-curated junction file when `star_twopass_junctions_file` is set.
- `params.json` was reorganized with a dedicated assembly section and a new polishing section (`from`, `assembly_by_chr`, `assembler`, `metassembler`, `do_twopass_polish`, etc.).
- CI pins Nextflow `25.10.0`, validates every test profile, and runs the default end-to-end smoke test on each push.

### Fixed

- Aletsch runs on single-chromosome chunked input: when a BAM contains only one populated chromosome, the process now selects that chromosome from the original full BAM instead of consuming the split chunk, which Aletsch 1.1.x would otherwise fail to flush.
- RENAME_GTF no longer deletes the input GTF (or its realpath target) after renaming.
- Stub behavior of the Aletsch module now produces the expected GTF output shape.
- Default strandedness metadata is reported as `unstranded` rather than `paired_end`.

---

## [0.0.19] - 2026-07-23

### Added

- New `--from bigbed` checkpoint that allows resuming the pipeline directly from the BigBed conversion step, skipping both meta-assembly and polishing entirely. This is useful when BED files are already available from a prior run and only the BigBed conversion (or re-conversion with different settings) is needed.
- `SORT_BED` module (`modules/custom/sort/main.nf`) that performs coordinate-based sorting (`sort -k1,1 -k2,2n -k3,3n`) on BED files prior to BigBed conversion. Sorting is now applied to full-length transcripts, scraps, and fusions within the polishing subworkflow, ensuring BED input compatibility with `bedToBigBed`.
- BigBed conversion step added to the `--from polish` checkpoint. Previously, BigBed conversion only ran during the full pipeline; now users resuming from polishing also get BigBed output for all categories (HQ, retentions, strong RTs, weak RTs, artifacts, fusions, scraps).
- `prefix` parameter (`params.prefix`) to override the output filename prefix used by the Beaver process. When set, Beaver output files will use this prefix instead of the default.
- `test` profile in `nextflow.config` for running xasm against a minimal synthetic dataset (`test_data/`). The profile configures a single-chromosome test genome, a test annotation GTF, and a test metassembly GTF, with BigBed conversion disabled by default.
- Test data files: `test_data/test_genome.fa`, `test_data/test_annotation.gtf`, and `test_data/test_metassembly.gtf`.
- `.gitignore` entries to preserve `test_data/` and its contents while ignoring other test-related files.

### Changed

- Disabled automatic work directory cleanup (`cleanup = false`) to preserve intermediate files for debugging and pipeline resumption.
- Updated `workflow.onComplete` handler to correctly distinguish between BED and BigBed output directories (`10_final` vs `11_bbs`) when reporting pipeline results, and to handle the case where BigBed conversion is active.
- Updated `params.json` and the task handler script (`assets/scripts/xasm.sh`) to expose the new `all_bed_path` parameter for the `--from bigbed` checkpoint and the `prefix` parameter.
- PublishDir pattern for BigTools processes corrected from `*.bed` to `*.bb` in `nextflow.config`.
- Added `FROM_BIGBED` publishDir configuration to write BigBed output to `11_bbs`.

### Fixed

- Groovy syntax errors in the `onComplete` handler where `and`/`or` keywords were used instead of the correct `&&`/`||` operators, which would cause runtime failures on pipeline completion.
- Fixed the `onComplete` handler referencing a non-existent `10_results` directory instead of the correct `10_final` directory.

---

## [0.0.18a] - 2026-06-16

### Added

- BigBed conversion (`bedToBigBed`) now runs by default for all transcript categories (HQ, retentions, strong RTs, weak RTs, artifacts, fusions, scraps) during the full pipeline run.
- New `11_bbs` output directory for storing BigBed (`.bb`) files alongside the existing `10_final` BED output.
- New `SORT_BED` instances for full-length transcripts, scraps, and fusions in the polishing subworkflow to ensure sorted BED input for downstream tools.

### Changed

- Updated the shell task handler (`assets/scripts/xasm.sh`) to assert correct input format.
- Changed `process_medium` resource label for the intron process.

---

## [0.0.18] - 2026-06-10

### Added

- Retained intronic transcript (RT) and artifact detection and classification implementation in the polishing subworkflow.
- Explicit publish step for classified transcript categories.

### Changed

- Adjusted resource allocations for multi-threaded processes to improve stability under concurrent execution.

### Fixed

- Retention transcripts now emit sorted (striped) output correctly.
- xloci intron extraction now applies `--unmask` as intended.
- Array handling fixes for GTF/BED input channels in the polishing subworkflow.
- Versioning channel propagation to `--from` checkpoints.
- Error strategy updated to `retry` for transient failures in external tool processes.

---

## [0.0.17] - 2026-05-20

### Added

- `--from polish` checkpoint for resuming the pipeline from the polishing step, skipping meta-assembly.
- JSON-based parameter input support.

### Changed

- Major refactor of the pipeline workflow structure and subworkflow organization.

---

## [0.0.16] - 2026-05-10

### Added

- Detach and strip operations for transcript classification.
- isotoools UTR processing integration.
- Array-based input support for select processes.

### Changed

- Reverted isotoools UTR behavior to previous state while retaining the detach/strip infrastructure.

### Fixed

- Compliance with `nf-core` linting standards.

---

## [0.0.15] - 2026-04-28

### Added

- GTF and BED file acceptance as input formats.
- Fusion transcript output as a default pipeline output.
- Pipeline DAG visualization output.
- Gene prediction linting step.
- Copyright headers.
