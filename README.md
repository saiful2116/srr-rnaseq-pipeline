# srr-rnaseq-pipeline

Analyze public mouse RNA-seq data **starting from nothing but SRR numbers**. No need to download FASTQ files yourself: give the pipeline a list of SRR accessions and it does the rest.

**What you get:**

1. **Aligned BAMs from SRR numbers.** Each sample is fetched from SRA, trimmed, aligned with STAR (mm10 / GRCm38, GENCODE vM25) and turned into a sorted, indexed BAM plus a bigWig track. Reads are streamed through a temporary folder and deleted after each sample, so raw FASTQ files never pile up on your disk.
2. **A count matrix ready for DESeq2.** featureCounts runs on all samples and writes `featurecounts/count_matrix.tsv`: one row per gene, one column per SRR.
3. **BAMs for your genes of interest.** Each BAM is filtered to the regions in a BED file (`GENE_BED`), giving small BAMs that contain only the genes you care about, ideal for IGV or sharing.

```
SRR list → download → trimming (fastp / cutadapt for SMART-seq) → STAR → BAM + bigWig
                                                    ↓ all samples
                         featureCounts matrix (DESeq2) · gene-list BAMs · correlation + PCA
```

## Features

- **Automatic trimming choice:** checks the first 200,000 reads for the SMART oligo. Above 0.5 % it uses cutadapt with SMART-seq adapters, otherwise fastp with adapter auto-detection.
- **Paired-end or single-end** detected automatically.
- **Strandedness** inferred from STAR `ReadsPerGene` and passed to featureCounts.
- **Resumable:** samples that already have an indexed BAM are skipped.
- **Saves disk space:** FASTQ, `.sra` and STAR temp files are deleted after each sample. BAMs and QC reports are kept.

## Requirements

Tools on your `PATH`:

| Tool | Used for |
|---|---|
| [sra-tools](https://github.com/ncbi/sra-tools) (`prefetch`, `fasterq-dump`) | download |
| [fastp](https://github.com/OpenGene/fastp) | trimming (standard libraries) |
| [cutadapt](https://cutadapt.readthedocs.io) | trimming (SMART-seq) |
| [STAR](https://github.com/alexdobin/STAR) | alignment |
| [samtools](http://www.htslib.org) | indexing, flagstat, gene-list BAMs |
| [deepTools](https://deeptools.readthedocs.io) | bigWig, correlation, PCA |
| [subread](https://subread.sourceforge.net) (`featureCounts`) | gene counts |
| `wget`, `awk`, `gunzip` | reference download |

All are available from bioconda:

```bash
conda create -n rnaseq -c bioconda -c conda-forge \
  sra-tools fastp cutadapt star samtools deeptools subread
conda activate rnaseq
```

**Hardware:** building the STAR index needs ~32 GB RAM (~1 h, done once). Alignment uses up to 16 GB for BAM sorting. Default: 12 threads.

## Setup

Edit the settings at the top of `srr_pipeline.sh`:

| Variable | Meaning |
|---|---|
| `THREADS` | CPU threads (default 12) |
| `REF_DIR` | shared reference folder (default `~/ref`). Genome, GTF and STAR index are downloaded/built here if missing |
| `GENE_BED` | BED file (chr, start, end) of mm10 regions for the gene-list BAMs |
| `STRAND` | `auto`, or force `0` unstranded / `1` stranded / `2` reverse |
| `TRIM` | `auto`, or force `fastp` / `smart` |

## Usage

```bash
TMPDIR=~/tmp bash srr_pipeline.sh SRR_Acc_List.txt [output_folder]
```

- `SRR_Acc_List.txt`: one SRR per line (extra columns are ignored, `#` lines are skipped). You can download it from the SRA Run Selector.
- `output_folder`: defaults to `results`.

Run once per dataset, e.g.

```bash
bash srr_pipeline.sh keratinocytes.txt results_keratinocytes
```

## Using the counts in DESeq2

```r
library(DESeq2)
counts  <- read.delim("results/featurecounts/count_matrix.tsv", row.names = 1, check.names = FALSE)
coldata <- data.frame(row.names = colnames(counts),
                      condition = c("ctrl", "ctrl", "treat", "treat"))  # one entry per SRR, same order
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
```

## Output

```
results/
├── bam/            sorted + indexed BAMs
├── bigwig/         CPM-normalized coverage tracks
├── star/           STAR logs, ReadsPerGene, splice junctions
├── fastp/          fastp HTML/JSON reports
├── cutadapt/       cutadapt reports (SMART-seq samples)
├── logs/           trimming choice, flagstat, strandedness
├── deeptools/      correlation heatmap, PCA (≥ 2 samples)
├── featurecounts/  counts.txt and count_matrix.tsv (ready for DESeq2/edgeR)
└── genelist_bam/   BAMs containing only the GENE_BED regions
```
