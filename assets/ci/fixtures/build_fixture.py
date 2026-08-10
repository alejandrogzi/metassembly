#!/usr/bin/env python3
"""Build the compact, deterministic mm39-derived xasm E2E fixture."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import random
import shutil
import subprocess
import tempfile
from pathlib import Path


LOCI = (
    {
        "source_chrom": "chr7",
        "source_start": 16_693_603,
        "source_end": 16_698_532,
        "window_start": 16_693_103,
        "window_end": 16_699_032,
        "chrom": "chrTestA",
        "transcript_id": "ENSMUST00000038163.8",
        "gene_id": "test_gene_plus",
        "strand": "+",
        "exon_sizes": (213, 1289, 2588),
        "exon_offsets": (0, 530, 2341),
    },
    {
        "source_chrom": "chr1",
        "source_start": 74_324_088,
        "source_end": 74_326_613,
        "window_start": 74_323_588,
        "window_end": 74_327_113,
        "chrom": "chrTestB",
        "transcript_id": "ENSMUST00000113805.8",
        "gene_id": "test_gene_plus_b",
        "strand": "+",
        "exon_sizes": (208, 169, 363),
        "exon_offsets": (0, 912, 2162),
    },
)

READ_LENGTH = 100
# Overlapping mates provide the compact insert-size observations required by
# Aletsch's profiler while remaining representative of short-fragment RNA-seq.
INSERT_SIZE = 180
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def read_fasta(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines() if not line.startswith(">"))


def exons_for(locus: dict) -> list[tuple[int, int]]:
    offset = locus["source_start"] - locus["window_start"]
    return [
        (offset + block_start, offset + block_start + block_size)
        for block_start, block_size in zip(locus["exon_offsets"], locus["exon_sizes"])
    ]


def transcript_sequence(genome: dict[str, str], locus: dict, skip_middle: bool = False) -> str:
    exons = exons_for(locus)
    if skip_middle:
        exons = [exon for index, exon in enumerate(exons) if index != len(exons) // 2]
    sequence = "".join(genome[locus["chrom"]][start:end] for start, end in exons)
    return sequence if locus["strand"] == "+" else reverse_complement(sequence)


def fragments(sequence: str, count: int, phase: int) -> list[str]:
    span = len(sequence) - INSERT_SIZE
    if span <= 0:
        raise ValueError("Transcript is shorter than the configured insert size")
    return [sequence[(phase + index * 37) % span :][:INSERT_SIZE] for index in range(count)]


def add_pair(records: list[tuple[str, str, str, str]], name: str, fragment: str, quality: str = "I") -> None:
    read1 = fragment[:READ_LENGTH]
    read2 = reverse_complement(fragment[-READ_LENGTH:])
    # Keep identical query names in both mates; Aletsch pairs BAM records by QNAME.
    records.append((name, read1, quality * len(read1), "1"))
    records.append((name, read2, quality * len(read2), "2"))


def write_fastqs(output_dir: Path, genome: dict[str, str]) -> dict[str, dict[str, int]]:
    randomizer = random.Random(7122026)
    manifest: dict[str, dict[str, int]] = {}

    for sample_index, sample in enumerate(("sampleA", "sampleB")):
        records: list[tuple[str, str, str, str]] = []
        counts = {"canonical_pairs": 0, "alternate_pairs": 0, "flush_pairs": 0, "intron_pairs": 0,
                  "deacon_negative_pairs": 0, "alignment_negative_pairs": 0,
                  "low_quality_pairs": 0}

        for locus_index, locus in enumerate(LOCI):
            sequence = transcript_sequence(genome, locus)
            count = ((2500, 1800) if sample == "sampleA" else (2200, 1600))[locus_index]
            for index, fragment in enumerate(fragments(sequence, count, 11 * (sample_index + locus_index + 1))):
                add_pair(records, f"{sample}:canonical:{locus['chrom']}:{index}", fragment)
            counts["canonical_pairs"] += count

        if sample == "sampleB":
            alternate = transcript_sequence(genome, LOCI[1], skip_middle=True)
            for index, fragment in enumerate(fragments(alternate, 180, 19)):
                add_pair(records, f"{sample}:alternate:chrTestB:{index}", fragment)
            counts["alternate_pairs"] = 180

        # Aletsch 1.1.x flushes a reference when the next populated reference
        # starts. This test-only chromosome follows chrTestB in the BAM and is
        # excluded from assembly, ensuring the final real test chromosome is
        # exercised without adding a third assembly branch.
        for index, fragment in enumerate(fragments(genome["chrTestFlush"], 120, sample_index * 17)):
            add_pair(records, f"{sample}:flush:{index}", fragment)
        counts["flush_pairs"] = 120

        # Genomic fragments spanning the first intron create assembly/polishing negatives.
        locus = LOCI[0]
        first_exon_end = exons_for(locus)[0][1]
        intron_template = genome[locus["chrom"]][first_exon_end - 120:first_exon_end + 430]
        for index, fragment in enumerate(fragments(intron_template, 90, sample_index * 7)):
            add_pair(records, f"{sample}:intron:{index}", fragment)
        counts["intron_pairs"] = 90

        # These contain no target 31-mers and should be rejected by the target-genome Deacon index.
        contaminant = "".join(randomizer.choice("ACGT") for _ in range(900))
        for index, fragment in enumerate(fragments(contaminant, 60, sample_index * 13)):
            add_pair(records, f"{sample}:deacon_negative:{index}", fragment)
        counts["deacon_negative_pairs"] = 60

        # One target k-mer plus random sequence: retained by Deacon but below STAR's alignment threshold.
        target_seed = genome["chrTestA"][800:831]
        for index in range(35):
            random_tail = "".join(randomizer.choice("ACGT") for _ in range(INSERT_SIZE - len(target_seed)))
            add_pair(records, f"{sample}:alignment_negative:{index}", target_seed + random_tail)
        counts["alignment_negative_pairs"] = 35

        for index in range(25):
            low_quality = ("ACGT" * 63)[:INSERT_SIZE - len(ADAPTER)] + ADAPTER
            add_pair(records, f"{sample}:low_quality:{index}", low_quality, quality="!")
        counts["low_quality_pairs"] = 25

        for mate in ("1", "2"):
            path = output_dir / "reads" / f"{sample}_{mate}.fastq.gz"
            with path.open("wb") as raw_handle, gzip.GzipFile(
                filename="", mode="wb", fileobj=raw_handle, mtime=0
            ) as gzip_handle, io.TextIOWrapper(gzip_handle, newline="\n") as handle:
                for name, sequence, quality, record_mate in records:
                    if record_mate == mate:
                        handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")
        manifest[sample] = counts

    return manifest


def write_gtf(path: Path) -> None:
    with path.open("w") as handle:
        for locus in LOCI:
            exons = exons_for(locus)
            start, end = exons[0][0] + 1, exons[-1][1]
            attrs = f'gene_id "{locus["gene_id"]}"; transcript_id "{locus["transcript_id"]}";'
            handle.write(f'{locus["chrom"]}\tmm39\tgene\t{start}\t{end}\t.\t{locus["strand"]}\t.\tgene_id "{locus["gene_id"]}";\n')
            handle.write(f'{locus["chrom"]}\tmm39\ttranscript\t{start}\t{end}\t.\t{locus["strand"]}\t.\t{attrs}\n')
            ordered = exons if locus["strand"] == "+" else list(reversed(exons))
            for exon_number, (exon_start, exon_end) in enumerate(ordered, 1):
                exon_attrs = attrs + f' exon_number "{exon_number}";'
                handle.write(f'{locus["chrom"]}\tmm39\texon\t{exon_start + 1}\t{exon_end}\t.\t{locus["strand"]}\t.\t{exon_attrs}\n')
                handle.write(f'{locus["chrom"]}\tmm39\tCDS\t{exon_start + 1}\t{exon_end}\t.\t{locus["strand"]}\t0\t{exon_attrs}\n')


def write_score_tracks(output_dir: Path, bedgraph_to_bigwig: Path, chrom_sizes: Path) -> None:
    tracks: dict[str, list[tuple[str, int, int, float]]] = {
        "spliceAiAcceptorPlus": [], "spliceAiDonorPlus": [],
        "spliceAiAcceptorMinus": [], "spliceAiDonorMinus": [],
    }
    # Keep every strand-specific BigWig structurally valid even when this tiny
    # fixture has no true junction on that strand.
    for entries in tracks.values():
        entries.extend((locus["chrom"], 0, 1, 0.0) for locus in LOCI)
    for locus in LOCI:
        exons = exons_for(locus)
        for left, right in zip(exons, exons[1:]):
            if locus["strand"] == "+":
                tracks["spliceAiDonorPlus"].append((locus["chrom"], left[1] - 1, left[1], 0.99))
                tracks["spliceAiAcceptorPlus"].append((locus["chrom"], right[0], right[0] + 1, 0.99))
            else:
                tracks["spliceAiAcceptorMinus"].append((locus["chrom"], left[1] - 1, left[1], 0.99))
                tracks["spliceAiDonorMinus"].append((locus["chrom"], right[0], right[0] + 1, 0.99))

    for name, entries in tracks.items():
        bedgraph = output_dir / "spliceai" / f"{name}.bedGraph"
        entries.sort()
        bedgraph.write_text("".join(f"{chrom}\t{start}\t{end}\t{score}\n" for chrom, start, end, score in entries))
        run([str(bedgraph_to_bigwig), str(bedgraph), str(chrom_sizes), str(bedgraph.with_suffix(".bw"))])
        bedgraph.unlink()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-2bit", type=Path, required=True)
    parser.add_argument("--two-bit-to-fa", type=Path, required=True)
    parser.add_argument("--fa-to-two-bit", type=Path, required=True)
    parser.add_argument("--bedgraph-to-bigwig", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("assets/test/test_data/e2e"))
    args = parser.parse_args()

    output_dir = args.output.resolve()
    for directory in (output_dir, output_dir / "reads", output_dir / "spliceai"):
        directory.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="xasm-fixture-") as temporary:
        temporary_dir = Path(temporary)
        genome: dict[str, str] = {}
        combined_fasta = temporary_dir / "genome.fa"
        with combined_fasta.open("w") as output:
            for locus in LOCI:
                extracted = temporary_dir / f'{locus["chrom"]}.fa'
                run([
                    str(args.two_bit_to_fa), str(args.source_2bit), str(extracted),
                    f'-seq={locus["source_chrom"]}', f'-start={locus["window_start"]}',
                    f'-end={locus["window_end"]}', "-noMask",
                ])
                sequence = read_fasta(extracted).upper()
                genome[locus["chrom"]] = sequence
                output.write(f'>{locus["chrom"]}\n')
                output.write("\n".join(sequence[index:index + 80] for index in range(0, len(sequence), 80)) + "\n")
            flush_randomizer = random.Random(38163)
            flush_sequence = "".join(flush_randomizer.choice("ACGT") for _ in range(1200))
            genome["chrTestFlush"] = flush_sequence
            output.write(">chrTestFlush\n")
            output.write("\n".join(flush_sequence[index:index + 80] for index in range(0, len(flush_sequence), 80)) + "\n")
        run([str(args.fa_to_two_bit), str(combined_fasta), str(output_dir / "genome.2bit")])

    chrom_sizes = output_dir / "chrom.sizes"
    chrom_sizes.write_text(
        "".join(f'{chrom}\t{len(sequence)}\n' for chrom, sequence in genome.items())
    )
    write_gtf(output_dir / "annotation.gtf")
    (output_dir / "repeats.bed").write_text("chrTestA\t720\t900\ttest_repeat\t0\t+\n")
    read_counts = write_fastqs(output_dir, genome)
    write_score_tracks(output_dir, args.bedgraph_to_bigwig, chrom_sizes)

    produced = sorted(path for path in output_dir.rglob("*") if path.is_file())
    manifest = {
        "seed": 7122026,
        "read_length": READ_LENGTH,
        "insert_size": INSERT_SIZE,
        "source": "mm39",
        "loci": list(LOCI),
        "reads": read_counts,
        "sha256": {str(path.relative_to(output_dir)): sha256(path) for path in produced if path.name != "manifest.json"},
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
