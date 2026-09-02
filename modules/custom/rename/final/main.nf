/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RENAME_FINAL_TRANSCRIPTS — Replace UNKNOWN protein tags in final transcript ids.

    Transcripts are packed into CDS-overlap components with py-packbed. Within a
    component carrying UNKNOWN tags, every UNKNOWN takes the protein name supported
    by the most members (ties: alphabetical). Components without any UNKNOWN, or
    where every protein tag is UNKNOWN, are flushed unchanged. The record count is
    preserved; only the name column may change.

    NOTE: container-only (py-packbed ships as a maturin wheel, not on bioconda).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RENAME_FINAL_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/py-packbed:latest' }"
    // The image entrypoints `python3`, which would swallow Nextflow's
    // `/bin/bash -ue` wrapper (bash then tries to run /bin/bash as a script).
    // `/usr/bin/env` re-execs the wrapper; apptainer/singularity ignore the OCI
    // entrypoint on exec, so they need no override.
    containerOptions "${ workflow.containerEngine == 'singularity' ? '' : '--entrypoint /usr/bin/env' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.renamed.bed"), emit: reads
    tuple val(meta), path("*.renamed.tsv"), emit: tsv
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // python escapes are doubled (\\t etc.) because the GString would eat a single backslash
    """
    cat > rename_final_transcripts.py << 'PYEOF'
#!/usr/bin/env python3
'''
Replace UNKNOWN protein tags in final transcript ids using CDS components.

Transcripts are packed into CDS-overlap components with py-packbed. Within a
component that carries UNKNOWN tags, every UNKNOWN takes the protein name
supported by the most members (ties: alphabetical). Components without any
UNKNOWN, or where every protein tag is UNKNOWN, are flushed unchanged. The
output BED holds exactly the input records in input order; only column 4 may
change. The output TSV holds one old/new pair per changed record.

--selftest builds a synthetic BED exercising every branch and asserts on it.
'''

import argparse
import os
import sys
import tempfile
from collections import Counter

import packbed

UNKNOWN = "UNKNOWN"


def protein_tag(name):
    # The protein tag is the last @-separated component; untagged names never
    # carry UNKNOWN and are left untouched.
    return name.rsplit("@", 1)[1] if "@" in name else None


def rename_map(bed):
    changes = {}
    packed = packbed.pack([bed], ["query"], overlap_type="cds")
    for group in packed.values():
        for component in group:
            members = [gene.name for gene in component]
            tags = [protein_tag(name) for name in members]
            if UNKNOWN not in tags:
                continue
            supported = [tag for tag in tags if tag is not None and tag != UNKNOWN]
            if not supported:
                continue
            counts = Counter(supported)
            winner = min(counts, key=lambda tag: (-counts[tag], tag))
            for name in members:
                if protein_tag(name) == UNKNOWN:
                    changes[name] = name[: name.rfind("@")] + "@" + winner
    return changes


def apply_map(bed, changes, output_bed, output_tsv):
    total = 0
    renamed = 0
    ids = set()
    with open(bed) as source, open(output_bed, "w") as sink, open(output_tsv, "w") as log:
        for line in source:
            if not line.strip():
                continue
            fields = line.rstrip("\\n").split("\\t")
            old = fields[3]
            new = changes.get(old, old)
            if new in ids:
                raise ValueError("renaming produced a duplicate transcript id: " + new)
            ids.add(new)
            if new != old:
                log.write(f"{old}\\t{new}\\n")
                renamed += 1
                fields[3] = new
            sink.write("\\t".join(fields) + "\\n")
            total += 1
    print(f"renamed {renamed} of {total} transcripts")
    return total, renamed


def selftest():
    fixture = {
        # mixed component: UNKNOWN takes the majority protein
        "chr1": ["A#OR1@UNKNOWN", "B#OR2@UNKNOWN", "C#OR3@PROT2", "D#OR4@PROT2"],
        # every tag UNKNOWN: left as is
        "chr2": ["E#OR1@UNKNOWN", "F#OR2@UNKNOWN"],
        # no UNKNOWN: flushed unchanged
        "chr3": ["G#OR1@PROT4", "H#OR2@PROT4"],
        # untagged name votes for nothing and stays untouched
        "chr4": ["I#OR1@UNKNOWN", "J#OR2@PROT5", "K-untagged"],
        # tied support: alphabetical wins
        "chr5": ["L#OR1@UNKNOWN", "M#OR2@PA", "N#OR3@PB"],
    }
    with tempfile.TemporaryDirectory() as work:
        bed = os.path.join(work, "fixture.bed")
        with open(bed, "w") as handle:
            for chrom, names in fixture.items():
                for name in names:
                    fields = [chrom, "0", "200", name, "0", "+", "50", "150", "0,0,0", "2", "100,100,", "0,0,"]
                    handle.write("\\t".join(fields) + "\\n")
        expected = {
            "A#OR1@UNKNOWN": "A#OR1@PROT2",
            "B#OR2@UNKNOWN": "B#OR2@PROT2",
            "I#OR1@UNKNOWN": "I#OR1@PROT5",
            "L#OR1@UNKNOWN": "L#OR1@PA",
        }
        assert rename_map(bed) == expected, rename_map(bed)

        out_bed = os.path.join(work, "out.bed")
        out_tsv = os.path.join(work, "out.tsv")
        total, renamed = apply_map(bed, expected, out_bed, out_tsv)
        assert (total, renamed) == (14, 4), (total, renamed)
        with open(out_bed) as handle:
            out_lines = handle.read().splitlines()
        assert len(out_lines) == 14
        assert out_lines[0].split("\\t")[3] == "A#OR1@PROT2"
        assert out_lines[4].split("\\t")[3] == "E#OR1@UNKNOWN"
        assert out_lines[10].split("\\t")[3] == "K-untagged"
        assert out_lines[11].split("\\t")[3] == "L#OR1@PA"
        with open(out_tsv) as handle:
            assert len(handle.read().splitlines()) == 4
    print("selftest ok")


def main():
    parser = argparse.ArgumentParser(
        description="Rename UNKNOWN protein tags in a final transcripts BED12."
    )
    parser.add_argument("--bed")
    parser.add_argument("--output-bed", default="renamed.bed")
    parser.add_argument("--output-tsv", default="renamed.tsv")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        selftest()
        return
    if not args.bed:
        parser.error("--bed is required")
    apply_map(args.bed, rename_map(args.bed), args.output_bed, args.output_tsv)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.exit(f"ERROR: {error}")
PYEOF
    python3 rename_final_transcripts.py --selftest

    python3 rename_final_transcripts.py \\
        --bed ${bed} \\
        --output-bed ${prefix}.renamed.bed \\
        --output-tsv ${prefix}.renamed.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        py-packbed: \$(python3 -c "import importlib.metadata as m; print(m.version('py-packbed'))" 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.renamed.bed ${prefix}.renamed.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        py-packbed: \$(python3 -c "import importlib.metadata as m; print(m.version('py-packbed'))" 2>/dev/null || echo unknown)
    END_VERSIONS
    """
}
