# End-to-end fixture

This directory contains a deterministic, committed fixture derived from two
high-confidence mm39 coding transcripts on different chromosomes. Coordinates
are remapped to `chrTestA` and `chrTestB`; no full source genome is needed to
run the tests. A small synthetic `chrTestFlush` carries reads solely to work
around Aletsch 1.1.x's final-reference flush behavior and is excluded from
assembly by the test profiles.

The fixture contains two paired-end samples, a compact 2bit genome, GTF
annotation, repeats, four SpliceAI-compatible BigWig tracks, and a
selenocysteine-sites BED (TGA codons, including two intronic sites on the test
loci) that drives XORF's selenocysteine masking. Read names
encode their intended class. `manifest.json` records source coordinates,
designed counts, seed, and checksums.

Regeneration requires UCSC `twoBitToFa`, `faToTwoBit`, and
`bedGraphToBigWig`:

```bash
python assets/ci/fixtures/build_fixture.py \
  --source-2bit /path/to/mm39.2bit \
  --two-bit-to-fa /path/to/twoBitToFa \
  --fa-to-two-bit /path/to/faToTwoBit \
  --bedgraph-to-bigwig /path/to/bedGraphToBigWig
```

Regeneration is intentionally separate from normal tests. Review every
manifest and golden-count change.
