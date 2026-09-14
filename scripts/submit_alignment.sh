#!/usr/bin/env bash
#
# submit_alignment.sh
#
# Auto-discovers every sample folder under $FASTQ_DIR (config.sh) and
# submits align_sample.sh for each one via sbatch, applying the SLURM
# resource settings from config.sh as command-line overrides.
#
# USAGE:
#   bash submit_alignment.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
source "${REPO_ROOT}/config/config.sh"

mkdir -p "${WORKDIR}/logs"

# Auto-discover sample folders. Filenames are matched by glob
# (*_1.fq.gz / *_2.fq.gz) rather than an exact pattern, so naming
# conventions beyond that suffix don't matter.
SAMPLES=()
for dir in "${FASTQ_DIR}"/*/; do
    sample=$(basename "${dir}")
    r1=$(find "${dir}" -maxdepth 1 -name "*_1.fq.gz" | sort | head -1)
    r2=$(find "${dir}" -maxdepth 1 -name "*_2.fq.gz" | sort | head -1)
    if [ -n "${r1}" ] && [ -n "${r2}" ]; then
        SAMPLES+=("${sample}")
    else
        echo "WARNING: skipping ${sample} - no *_1.fq.gz/*_2.fq.gz found" >&2
    fi
done

echo "Found ${#SAMPLES[@]} sample(s) to submit."

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Submitting ${SAMPLE}..."
    sbatch \
        --partition="${SLURM_PARTITION}" \
        --cpus-per-task="${SLURM_CPUS}" \
        --mem="${SLURM_MEM}" \
        --time="${SLURM_TIME}" \
        --output="${WORKDIR}/logs/%x_%A_%a.out" \
        "${SCRIPT_DIR}/align_sample.sh" "${SAMPLE}"
done

echo ""
echo "All ${#SAMPLES[@]} job(s) submitted. Check status with: squeue -u \$USER"
echo "Logs will appear in: ${WORKDIR}/logs/"
