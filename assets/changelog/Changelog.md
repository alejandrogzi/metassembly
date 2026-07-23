# Changelog

All notable changes to xasm will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
