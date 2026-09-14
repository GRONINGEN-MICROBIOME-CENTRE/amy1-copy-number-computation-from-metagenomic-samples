# AMY1 Copy Number Estimation Pipeline

Estimates AMY1 diploid copy number from short-read WGS data using the
read-depth-ratio method described in Kamitaki et al., *Nature* 2026
(SPARK/All of Us cohorts). Originally built for low-depth
(0.5-5x) host-extracted metagenomic WGS data, but works on any
short-read WGS BAM.

## Background

AMY1A, AMY1B and AMY1C sit inside a segmental duplication with >99%
sequence identity between the three gene copies. Reads originating from
this locus cannot be mapped with high confidence and receive low
mapping-quality (MAPQ) scores. **Most standard alignment pipelines
apply a MAPQ filter (commonly MAPQ >= 20) that silently discards nearly
all reads from this locus** - this is invisible in normal use, but it
makes AMY1 copy-number estimation from the resulting BAMs impossible
(in our testing: essentially zero informative reads survive).

This pipeline re-aligns from FASTQ **without** a MAPQ filter, then
counts reads across the three AMY1A/B/C intervals and normalizes
against flanking single-copy control regions (0.5 Mb upstream of
AMY2B + 0.5 Mb downstream of AMY1C). The resulting ratio is intended to
scale with an individual's total AMY1 diploid copy number.

## Requirements

- A SLURM cluster
- Conda environment with: `trimmomatic`, `bowtie2`, `samtools`
- A GRCh38 reference (no ALT contigs) with a prebuilt Bowtie2 index
- Paired-end short-read WGS FASTQs

See `environment.yml` for a conda environment that satisfies the tool
requirements.

## Setup

```bash
git clone <this-repo-url>
cd amy1-copy-number-pipeline
cp config/config.sh.example config/config.sh
nano config/config.sh   # fill in your paths, conda env, SLURM partition, etc.
```

Input layout expected under `$FASTQ_DIR`: one subfolder per sample,
each containing paired FASTQs matched by glob (`*_1.fq.gz` /
`*_2.fq.gz`) - exact filenames beyond that suffix don't matter.

```
$FASTQ_DIR/
├── sample01/
│   ├── sample01_1.fq.gz
│   └── sample01_2.fq.gz
├── sample02/
│   ├── sample02_run1_1.fq.gz
│   └── sample02_run1_2.fq.gz
...
```

## Usage

```bash
# 1. Align every sample (submits one SLURM job per sample)
bash scripts/submit_alignment.sh
squeue -u $USER   # wait for all jobs to finish

# 2. Count AMY1 vs. flanking reads for every finished BAM
bash scripts/count_amy1_reads.sh
```

Both `submit_alignment.sh`/`align_sample.sh` and `count_amy1_reads.sh`
are resume-safe: re-running after an interruption skips samples already
completed rather than redoing them.

Output: `$WORKDIR/amy1_readcounts.tsv`

| Column | Meaning |
|---|---|
| `sample` | Sample ID (folder name under `$FASTQ_DIR`) |
| `amy1a_reads`, `amy1b_reads`, `amy1c_reads` | Reads in each of the three AMY1 gene intervals |
| `amy1_total` | Sum of the three above |
| `flank_up_reads`, `flank_down_reads` | Reads in the two flanking control windows |
| `flank_total` | Sum of the two above |
| `amy1_flank_ratio` | `amy1_total / flank_total` - the core signal |
| `total_mapped_reads` | Total mapped reads in the BAM (context/QC) |
| `amy1_zero_flag` | `YES` if `amy1_total == 0` (likely insufficient depth for this sample) |

## Known limitations

- **This ratio is not a calibrated copy number.** Converting the ratio
  into an integer copy number requires GC-bias correction and an
  integer-centering calibration step (as used in the paper's UKB
  pipeline) - not included here. Treat `amy1_flank_ratio` as a relative
  signal, not an absolute copy-number call, unless you've added that
  calibration yourself.
- **Depth matters a lot below ~1x.** In our validation (downsampling
  real high-depth WGS samples), ratio accuracy stayed within ~2-4% of
  the full-depth value from 1x up to 10x, but degraded sharply below
  1x (average error ~17% at 0.5x, up to 27% in one sample). Below ~1x,
  treat individual estimates with real caution.
- **Extraction method may matter, separately from depth.** In samples
  from a metagenomic host-extraction pipeline, we observed the flanking
  control region's coverage running below genome-wide average more
  often than in matched clean WGS samples at the same depth - this is
  a preliminary finding (not yet confirmed across a large sample set)
  but suggests extraction-method-specific effects beyond simple
  depth/counting noise are possible and worth checking for your own
  data.
- **Biological replicates should agree.** AMY1 copy number is a fixed
  germline trait - if you have multiple samples from the same
  individual (e.g. different timepoints), their ratios should match
  within counting noise. Large disagreement is a useful QC signal.

## Citation

Kamitaki, N. et al. Diversification and prevalence of a *AMY1* copy
number-driven digestive advantage in human populations. *Nature*
(2026).
