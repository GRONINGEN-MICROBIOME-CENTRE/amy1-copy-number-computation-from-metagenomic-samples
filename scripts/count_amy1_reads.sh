#!/usr/bin/env bash
#
# count_amy1_reads.sh
#
# Run AFTER all align_sample.sh jobs finish (check with: squeue -u $USER).
# Counts reads in AMY1A/B/C vs. flanking control regions for every BAM
# produced by the alignment step, and writes one summary row per sample.
#
# Method (Kamitaki et al., Nature 2026): sum reads across the three
# AMY1A/AMY1B/AMY1C intervals, normalize against reads in flanking
# single-copy regions (0.5 Mb upstream of AMY2B + 0.5 Mb downstream of
# AMY1C). The resulting ratio is intended to scale with AMY1 diploid
# copy number. Coordinates below are GRCh38 and are fixed by the assay
# itself - they are not something you should need to edit.
#
# USAGE:
#   bash count_amy1_reads.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
source "${REPO_ROOT}/config/config.sh"

BAM_DIR="${WORKDIR}/02_aligned_nomapqfilter"
OUT="${WORKDIR}/amy1_readcounts.tsv"

# AMY1A/B/C coordinates (GRCh38, from the paper's Methods)
AMY1A="chr1:103638545-103666411"
AMY1B="chr1:103685558-103713427"
AMY1C="chr1:103732687-103760549"

# Flanking / normalizer regions: 0.5 Mb upstream of AMY2B, 0.5 Mb downstream of AMY1C
FLANK_UP="chr1:103053815-103553815"
FLANK_DOWN="chr1:103760549-104260549"

if [ ! -f "$OUT" ]; then
    echo -e "sample\tamy1a_reads\tamy1b_reads\tamy1c_reads\tamy1_total\tflank_up_reads\tflank_down_reads\tflank_total\tamy1_flank_ratio\ttotal_mapped_reads\tamy1_zero_flag" > "$OUT"
fi

for bam in "$BAM_DIR"/*.nomapqfilter.bam; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" .nomapqfilter.bam)

    if grep -qP "^${sample}\t" "$OUT" 2>/dev/null; then
        echo "Skipping ${sample} - already in output" >&2
        continue
    fi

    echo "Processing ${sample}..." >&2

    a=$(samtools view -c "$bam" "$AMY1A")
    b=$(samtools view -c "$bam" "$AMY1B")
    c=$(samtools view -c "$bam" "$AMY1C")
    amy1_total=$((a + b + c))

    fu=$(samtools view -c "$bam" "$FLANK_UP")
    fd=$(samtools view -c "$bam" "$FLANK_DOWN")
    flank_total=$((fu + fd))

    total_mapped=$(samtools idxstats "$bam" | awk '{sum+=$3} END {print sum}')

    ratio=$(awk -v a="$amy1_total" -v f="$flank_total" 'BEGIN { if (f>0) printf "%.4f", a/f; else print "NA" }')

    zero_flag="no"
    if [ "$amy1_total" -eq 0 ]; then
        zero_flag="YES"
    fi

    echo -e "${sample}\t${a}\t${b}\t${c}\t${amy1_total}\t${fu}\t${fd}\t${flank_total}\t${ratio}\t${total_mapped}\t${zero_flag}" >> "$OUT"
done

echo "" >&2
echo "Done. Results written to $OUT" >&2
echo "" >&2
column -t "$OUT"
