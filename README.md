
<p align="center">
  <p align="center">
    <img width=200 align="center" src="./assets/figures/hillerlab.png" >
  </p>

  <span>
    <h1 align="center">
        xasm
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/xasm" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/xasm?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> meta-assemble transcriptomes from bulk RNA-seq datasets at scale  </span>
        <br>
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

<div align="center">

<pre style="font-size: 18px;">

[ QC + trim ] ──► [ decontaminate ] ──► [ STAR-2pass align ]

S1  ───────●●──────────●●●──────────────────●●────────────
S2  ──────────────●●●────────●●─────────────●●●───────────
S3  ────●●──────────────────●●●●────────────────●─────────
S4  ─────────●●●────────────●●───────────●●●──────────────
S5  ─────────────●●────────────────●●●──────────●●────────

▼

═══════════════●●════════════●●●●════════════●●═══════════

└────────────────── metassembled transcriptome ──────────────┘
</pre>

</div>

---

## Usage

> [!NOTE]
> Requirements: Nextflow ≥ 25.04.6, Docker or Apptainer, Java.

```bash
git clone https://github.com/hillerlab/xasm.git
cd xasm
```

Edit `params.json` (set `genome`, `assembly_prefix`), then:

```bash
# Docker
nextflow run main.nf -params-file params.json -profile docker

# Apptainer / Singularity
nextflow run main.nf -params-file params.json -profile apptainer
```

Smoke test:
```bash
nextflow run main.nf -profile test,apptainer
```

> [!NOTE]
> You can also specify these options directly in `params.json`.

A helper sh script is provided to run the pipeline on a SLURM cluster. See details below.

<details>
<summary>Click to expand</summary>


Edit the path variables at the top of `assets/hpc/do_xasm.sh` (cache dir, container image, manifest path), then submit:

```bash
sbatch --array=1-<N> do_xasm.sh
```

Each array task spawns one Nextflow head job that submits all compute as child SLURM jobs.

STAR_ALIGN_1PASS, STAR_ALIGN_2PASS, FASTP, DEACON_FILTER, and BEDGRAPHTOBIGWIG run as SLURM job arrays. Partition routing, array sizes, and resource tiers are documented inline in `nextflow.config` — edit there to match your cluster.

</details>

---

## Output

```
results/
├── 00_prepare/
│   ├── deacon_index/    *fasta
│   └── wget/            *2bit
├── 01_deacon_filter/    *fastq
├── 02_star_index/       star/
├── 03_star_2pass/       *bam/*bai
├── 04_junctions/        *bed
├── 05_coverage/         *bigwig
├── 06_beaver/           *gtf
├── 07_remove_dirt/      *gtf
├── 08_gxf2bed/          *bed  <- --from polish checkpoint
├── 09_polish/
│   ├── fusions/          *bed
│   ├── orphans/          *bed
│   ├── intron_sequences/ *tsv
│   ├── iic/              *tsv
│   └── classify/         *tsv
├── 10_final/             *bed  <- --from bigbed checkpoint
├── 11_bbs/               *bb
└── pipeline_info/    timeline, trace, DAG
```

---

## Where to edit

| File | What |
|------|------|
| `params.json` | Genome paths, alignment settings, checkpoints — per run |
| `nextflow.config` | Compute resources, profiles, container, SLURM — rarely |
