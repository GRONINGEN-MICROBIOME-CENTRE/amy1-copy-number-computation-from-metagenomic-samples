#!/usr/bin/env bash
#SBATCH --job-name=amy1_align
#SBATCH --output=logs/%x_%A_%a.out
#
# align_sample.sh
#
# Trims adapters and aligns one sample's paired FASTQs to the human
# reference, WITHOUT filtering on mapping quality (MAPQ).
#
# Why no MAPQ filter: AMY1A, AMY1B and AMY1C sit inside a segmental
# duplication with >99% sequence identity between the three copies.
# Reads originating there cannot be mapped with high confidence and
# receive low MAPQ scores - a standard MAPQ >= 20 filter (the default
# in many pipelines, including whichever one produced your existing
# BAMs) silently discards nearly all reads from this locus before they
# ever reach a BAM file, making downstream copy-number estimation
# impossible. This script produces a BAM suitable specifically for
# AMY1 copy-number analysis - it is not a general-purpose replacement
# for your standard variant-calling/imputation alignment.
#
# Resource requests (partition/cpus/mem/time) are supplied by
# submit_alignment.sh via sbatch command-line flags, which override the
# placeholder #SBATCH lines above - edit config/config.sh, not this file,
# to change them.
#
# USAGE:
#   sbatch align_sample.sh <SAMPLE_ID>
#   (SAMPLE_ID = a folder name under $FASTQ_DIR, e.g. sample01)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
source "${REPO_ROOT}/config/config.sh"

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

if [ $# -lt 1 ]; then
    echo "ERROR: No sample ID provided"
    echo "Usage: $0 <SAMPLE_ID>"
    exit 1
fi
SAMPLE="$1"

echo "========================================"
echo " align_sample.sh - sample: ${SAMPLE}"
echo " Started: $(date)"
echo " Host: $(hostname)"
echo "========================================"

TRIM_DIR="${WORKDIR}/01_trimmed"
ALIGN_DIR="${WORKDIR}/02_aligned_nomapqfilter"
LOG_DIR="${WORKDIR}/logs"

TRIM_THREADS=8
BOWTIE2_THREADS="${SLURM_CPUS}"
SAMTOOLS_THREADS=8

mkdir -p "${TRIM_DIR}" "${ALIGN_DIR}" "${LOG_DIR}"

R1_SRC=$(find "${FASTQ_DIR}/${SAMPLE}" -maxdepth 1 -name "*_1.fq.gz" | sort | head -1)
R2_SRC=$(find "${FASTQ_DIR}/${SAMPLE}" -maxdepth 1 -name "*_2.fq.gz" | sort | head -1)
[ -z "${R1_SRC}" ] && { echo "ERROR: no *_1.fq.gz found under ${FASTQ_DIR}/${SAMPLE}"; exit 1; }
[ -z "${R2_SRC}" ] && { echo "ERROR: no *_2.fq.gz found under ${FASTQ_DIR}/${SAMPLE}"; exit 1; }
echo "[${SAMPLE}] R1: ${R1_SRC}"
echo "[${SAMPLE}] R2: ${R2_SRC}"

# =============================================================================
# STAGE 1: ADAPTER TRIMMING
# =============================================================================
TRIM_R1="${TRIM_DIR}/${SAMPLE}_R1_paired.fastq.gz"
TRIM_R2="${TRIM_DIR}/${SAMPLE}_R2_paired.fastq.gz"
TRIM_R1_U="${TRIM_DIR}/${SAMPLE}_R1_unpaired.fastq.gz"
TRIM_R2_U="${TRIM_DIR}/${SAMPLE}_R2_unpaired.fastq.gz"

if [ -f "${TRIM_R1}" ] && [ -s "${TRIM_R1}" ] && [ -f "${TRIM_R2}" ] && [ -s "${TRIM_R2}" ]; then
    echo "[${SAMPLE}] Already trimmed - skipping"
else
    echo "[${SAMPLE}] Trimming (Trimmomatic)..."
    trimmomatic PE \
        -threads "${TRIM_THREADS}" \
        -phred33 \
        "${R1_SRC}" "${R2_SRC}" \
        "${TRIM_R1}" "${TRIM_R1_U}" \
        "${TRIM_R2}" "${TRIM_R2_U}" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
    rm -f "${TRIM_R1_U}" "${TRIM_R2_U}"
fi

# =============================================================================
# STAGE 2: ALIGNMENT - Bowtie2, no MAPQ filter (see header note)
# =============================================================================
BAM="${ALIGN_DIR}/${SAMPLE}.nomapqfilter.bam"

if [ -f "${BAM}" ] && [ -s "${BAM}" ] && [ -f "${BAM}.bai" ]; then
    echo "[${SAMPLE}] BAM already exists - skipping"
else
    echo "[${SAMPLE}] Bowtie2 mapping to human genome (no MAPQ filter)..."

    bowtie2 \
        -x "${BOWTIE2_INDEX}" \
        -1 "${TRIM_R1}" \
        -2 "${TRIM_R2}" \
        -p "${BOWTIE2_THREADS}" \
        --no-unal \
        --very-sensitive \
        --rg-id "${SAMPLE}" \
        --rg "SM:${SAMPLE}" \
        --rg "PL:ILLUMINA" \
        --rg "LB:lib1" \
    | samtools view -bS -F 4 -f 2 \
    | samtools sort -@ "${SAMTOOLS_THREADS}" -o "${BAM}"

    samtools index "${BAM}"

    FLAGSTAT="${ALIGN_DIR}/${SAMPLE}_flagstat.txt"
    samtools flagstat "${BAM}" > "${FLAGSTAT}"
    MAPPED=$(grep "mapped (" "${FLAGSTAT}" | head -1 | awk '{print $1}')
    echo "[${SAMPLE}]   Mapped reads (no MAPQ filter): ${MAPPED}"
fi

echo "[${SAMPLE}] Done: $(date)"
echo "[${SAMPLE}] BAM: ${BAM}"
