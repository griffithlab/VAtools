#!/usr/bin/env bash
# Builds docs/_static/vatools-examples.tar.gz, the small example dataset
# that docs/*.rst example commands are written against. Not run by CI or
# the Sphinx build -- rerun by hand (from the repo root or anywhere) when
# tests/test_data/ changes in a way that should be reflected in the docs.
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root

WORKDIR=$(mktemp -d)
OUT="$WORKDIR/vatools-examples"
mkdir -p "$OUT"
TD=tests/test_data

# Files copied as-is
cp "$TD/snvs.bam_readcount"                                 "$OUT/sample.snv.bam_readcount"
cp "$TD/indels.bam_readcount"                                "$OUT/sample.indel.bam_readcount"
cp "$TD/kallisto.genes"                                      "$OUT/sample.kallisto.genes.tsv"
cp "$TD/kallisto.transcripts"                                 "$OUT/sample.kallisto.transcripts.tsv"
cp "$TD/kallisto.transcript_version.transcripts"               "$OUT/sample.kallisto.transcripts_versioned.tsv"
cp "$TD/transcripts.gtf"                                      "$OUT/sample.stringtie.transcripts.gtf"
cp "$TD/info.tsv"                                              "$OUT/sample.info_single.tsv"
cp "$TD/multi_col.tsv"                                         "$OUT/sample.info_multi.tsv"
# Trimmed to just the row matching sample.vcf's variant -- the source file's
# second row (22:18644674) has no matching variant in sample.vcf, and
# vep-annotation-reporter's -t merge crashes (uncaught KeyError) on TSV rows
# it can't find in the VCF, unlike transform-split-values -t which handles
# that gracefully.
head -2 "$TD/transform_split_values/variants.tsv"               > "$OUT/sample.variant_report.tsv"
cp "$TD/ref_transcript_mismatch_reporter/csq_mismatch.vcf"      "$OUT/sample.mismatch.vcf"

# Header-trimmed VCFs (see scripts/trim_example_vcf.py)
python3 scripts/trim_example_vcf.py "$TD/input.vcf"                   > "$OUT/sample.vcf"
python3 scripts/trim_example_vcf.py "$TD/multiple_samples.vcf"        > "$OUT/sample.multi_sample.vcf"
python3 scripts/trim_example_vcf.py "$TD/input.snvs_and_indels.vcf"   > "$OUT/sample.snvs_and_indels.vcf"
python3 scripts/trim_example_vcf.py "$TD/input.no_gt_in_format.vcf"   > "$OUT/sample.no_gt.vcf"
python3 scripts/trim_example_vcf.py "$TD/multiple_transcripts.vcf"    > "$OUT/sample.multi_transcript.vcf"

# Whole-transcriptome tables, trimmed to a handful of rows
{ head -1 "$TD/genes.tsv"; sed -n '2,4p' "$TD/genes.tsv"; grep 'ENSG00000184979' "$TD/genes.tsv"; } \
  > "$OUT/sample.stringtie.genes.tsv"
{ head -1 "$TD/genes.fpkm_tracking"; grep 'ENSG00000184979' "$TD/genes.fpkm_tracking"; } \
  > "$OUT/sample.cufflinks.genes.fpkm_tracking"
{ head -1 "$TD/isoforms.fpkm_tracking"; grep 'ENST00000215794' "$TD/isoforms.fpkm_tracking"; } \
  > "$OUT/sample.cufflinks.isoforms.fpkm_tracking"

# Hand-authored preferred-transcript TSV (no real fixture matches our variant)
printf 'transcript_id\nENST00000423297\n' > "$OUT/sample.preferred_transcripts.list.tsv"

tar czf docs/_static/vatools-examples.tar.gz -C "$WORKDIR" vatools-examples
rm -rf "$WORKDIR"
echo "Wrote docs/_static/vatools-examples.tar.gz"
