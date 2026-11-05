#!/usr/bin/env bash
# Per-sample RNA-seq pipeline from an SRR list (mouse mm10 / GRCm38, GENCODE vM25):
#   download -> trimming (fastp, or cutadapt for SMART-seq) -> STAR -> sorted BAM -> bigWig
# then, across all samples: deepTools correlation/PCA, featureCounts, gene-list BAMs.
# FASTQ and .sra files are deleted after each sample; QC reports and BAMs are kept.
# Gene list (step 9): GENE_BED below, a BED file of mm10 regions (chr, start, end), e.g. genes.bed.
# Usage: TMPDIR=~/amarRNAseq/tmp bash srr_pipeline.sh SRR_Acc_List.txt [output_folder]
#   SRR_Acc_List.txt: one SRR per line (extra columns are ignored)
#   Run once per dataset, e.g.  bash srr_pipeline.sh SRR_Acc_List.txt results
set -euo pipefail

SRR_LIST="$(readlink -f "${1:?Usage: $0 SRR_Acc_List.txt [output_folder]}")"
THREADS=12
REF_DIR="$HOME/ref"         # shared reference folder (all projects)
FASTA="$REF_DIR/GRCm38.primary_assembly.genome.fa"   # mm10 genome, only needed to build an index
GTF="$REF_DIR/gencode.vM25.annotation.gtf"           # annotation matching the index
STAR_INDEX="$REF_DIR/star_mm10_vM25"                      # existing STAR index
STRAND=auto   # featureCounts -s: auto (from STAR ReadsPerGene), or 0 = unstranded, 1 = stranded, 2 = reverse
TRIM=auto     # auto = cutadapt if SMART primer found in the reads, else fastp; or force: fastp / smart
SMART_THRESHOLD=0.5                                   # % of reads containing the SMART oligo -> SMART trimming
GENE_BED="$HOME/Path/to/the/genes.bed"   # gene regions for step 9 (BED: chr, start, end)
OUT="${2:-results}"
TMP="${TMPDIR:-/tmp}/srr_tmp"

mkdir -p "$OUT"/{fastp,cutadapt,bam,star,bigwig,logs} "$TMP" "$REF_DIR"

# SMART-seq library sequences, removed with cutadapt together with everything after them
SMART_ADAPTERS=( "A{12}"                          # poly-A stretch
                 ACTCTGCGTTGATACCACTGCTT          # SMART oligo, reverse complement
                 GTGTAGATCTCGGTGGTCGCCGTATCATT    # Illumina P5, reverse complement
                 CAAGCAGAAGACGGCATACGAGAT         # Illumina P7
                 AGATCGGAAGAGC                    # TruSeq universal
                 CTGTCTCTTATACACATCT )            # Nextera

# 0. Reference files (downloaded once)
GENCODE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25
if [[ ! -s "$GTF" ]]; then
  echo "Downloading GENCODE vM25 annotation..."
  wget -c -O "$GTF.gz" "$GENCODE/gencode.vM25.annotation.gtf.gz"
  gunzip -f "$GTF.gz"
fi

# 0b. Build the STAR index once (~32 GB RAM, ~1 h)
if [[ ! -s "$STAR_INDEX/SA" ]]; then
  echo "No STAR index found at $STAR_INDEX - building it (one time, ~1 hour)..."
  if [[ ! -s "$FASTA" ]]; then
    echo "Downloading mm10 (GRCm38) genome..."
    wget -c -O "$FASTA.gz" "$GENCODE/GRCm38.primary_assembly.genome.fa.gz"
    gunzip -f "$FASTA.gz"
  fi
  mkdir -p "$STAR_INDEX"
  STAR --runMode genomeGenerate --runThreadN "$THREADS" \
       --genomeDir "$STAR_INDEX" --genomeFastaFiles "$FASTA" \
       --sjdbGTFfile "$GTF" --sjdbOverhang 100 \
       --outFileNamePrefix "$STAR_INDEX/"
fi

while read -r SRR _; do
  SRR="${SRR//[$'\r\t ']/}"          # strip whitespace / Windows line endings
  [[ -z "$SRR" || "$SRR" == \#* ]] && continue

  if [[ -s "$OUT/bam/${SRR}.sorted.bam.bai" ]]; then
    echo "[$SRR] already done, skipping"; continue
  fi
  echo "[$SRR] start $(date)"
  WD="$TMP/$SRR"; rm -rf "$WD"; mkdir -p "$WD"

  # 1. Download
  prefetch "$SRR" -O "$WD" --max-size 200G
  fasterq-dump "$WD/$SRR" -O "$WD" -e "$THREADS" -t "$WD" --split-3

  # 2. Paired or single?
  if [[ -f "$WD/${SRR}_1.fastq" && -f "$WD/${SRR}_2.fastq" ]]; then
    R1="$WD/${SRR}_1.fastq"; R2="$WD/${SRR}_2.fastq"; PE=1
  else
    R1="$WD/${SRR}.fastq"; R2=""; PE=0
  fi

  # 4. Trimming
  #    4a. SMART-seq check: % of the first 200,000 reads that contain the SMART oligo
  PCT=$(head -n 800000 "$R1" | awk 'NR%4==2{n++; if($0 ~ /AAGCAGTGGTATCAACGCAGAGT|ACTCTGCGTTGATACCACTGCTT/) h++}
                                    END{printf "%.2f", (n? 100*h/n : 0)}')
  MODE=$TRIM
  [[ $MODE == auto ]] && MODE=$(awk -v p="$PCT" -v t="$SMART_THRESHOLD" 'BEGIN{print (p>=t ? "smart" : "fastp")}')
  echo "trimming=$MODE smart_oligo_pct=$PCT layout=$([[ $PE -eq 1 ]] && echo PE || echo SE)" | tee "$OUT/logs/${SRR}.trimming.txt"

  T1="$WD/${SRR}_1.trim.fq.gz"; T2=""; [[ $PE -eq 1 ]] && T2="$WD/${SRR}_2.trim.fq.gz"
  if [[ $MODE == smart ]]; then
    #  4b. cutadapt: remove poly-A, SMART oligo and Illumina adapters + everything after them
    A1=(); A2=(); for a in "${SMART_ADAPTERS[@]}"; do A1+=(-a "$a"); A2+=(-A "$a"); done
    if [[ $PE -eq 1 ]]; then
      cutadapt -j "$THREADS" -n 3 -m 20 -q 15 -e 0.1 "${A1[@]}" "${A2[@]}" \
               -o "$T1" -p "$T2" "$R1" "$R2" > "$OUT/cutadapt/${SRR}.cutadapt.txt"
    else
      cutadapt -j "$THREADS" -n 3 -m 20 -q 15 -e 0.1 "${A1[@]}" \
               -o "$T1" "$R1" > "$OUT/cutadapt/${SRR}.cutadapt.txt"
    fi
  else
    #  4c. fastp, adapter auto-detection (standard libraries)
    if [[ $PE -eq 1 ]]; then
      fastp -i "$R1" -I "$R2" -o "$T1" -O "$T2" \
            --detect_adapter_for_pe -w "$THREADS" \
            -j "$OUT/fastp/${SRR}.fastp.json" -h "$OUT/fastp/${SRR}.fastp.html"
    else
      fastp -i "$R1" -o "$T1" -w "$THREADS" \
            -j "$OUT/fastp/${SRR}.fastp.json" -h "$OUT/fastp/${SRR}.fastp.html"
    fi
  fi
  rm -f "$R1" ${R2:+"$R2"}           # raw FASTQ no longer needed

  # 6. STAR alignment (default filters) -> coordinate-sorted BAM (+ per-gene counts for the strandedness check)
  STAR --runThreadN "$THREADS" --genomeDir "$STAR_INDEX" \
       --readFilesIn $T1 $T2 --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 \
       --quantMode GeneCounts --outSAMattributes NH HI AS nM \
       --outTmpDir "$WD/STARtmp" --outFileNamePrefix "$WD/${SRR}."
  mv "$WD/${SRR}.Aligned.sortedByCoord.out.bam" "$OUT/bam/${SRR}.sorted.bam"
  mv "$WD/${SRR}.Log.final.out"        "$OUT/star/${SRR}.Log.final.out"
  mv "$WD/${SRR}.Log.out"              "$OUT/star/${SRR}.Log.out"
  mv "$WD/${SRR}.ReadsPerGene.out.tab" "$OUT/star/${SRR}.ReadsPerGene.out.tab"
  mv "$WD/${SRR}.SJ.out.tab"           "$OUT/star/${SRR}.SJ.out.tab"
  samtools index "$OUT/bam/${SRR}.sorted.bam"
  samtools flagstat "$OUT/bam/${SRR}.sorted.bam" > "$OUT/logs/${SRR}.flagstat.txt"

  # 6b. deepTools: CPM-normalized coverage track (bigWig) for genome browsers
  bamCoverage -b "$OUT/bam/${SRR}.sorted.bam" -o "$OUT/bigwig/${SRR}.CPM.bw" \
              --normalizeUsing CPM --binSize 10 -p "$THREADS"

  # 7. Remove FASTQ + .sra + STAR temp files for this sample
  rm -rf "$WD"
  echo "[$SRR] done $(date)"
done < "$SRR_LIST"

BAMS=( "$OUT"/bam/*.sorted.bam )

# 7b. deepTools across samples: correlation heatmap + PCA (needs >= 2 BAMs)
if (( ${#BAMS[@]} >= 2 )); then
  mkdir -p "$OUT/deeptools"
  LABELS=( $(for b in "${BAMS[@]}"; do basename "$b" .sorted.bam; done) )
  multiBamSummary bins -b "${BAMS[@]}" --labels "${LABELS[@]}" \
                  -o "$OUT/deeptools/bam_summary.npz" -p "$THREADS"
  plotCorrelation -in "$OUT/deeptools/bam_summary.npz" --corMethod spearman \
                  --whatToPlot heatmap --plotNumbers --skipZeros \
                  -o "$OUT/deeptools/correlation_heatmap.pdf" \
                  --outFileCorMatrix "$OUT/deeptools/correlation_matrix.tab"
  plotPCA -in "$OUT/deeptools/bam_summary.npz" \
          -o "$OUT/deeptools/PCA.pdf" --outFileNameData "$OUT/deeptools/PCA.tab"
fi

# 7c. strandedness from STAR ReadsPerGene (column 3 = forward, column 4 = reverse)
if [[ $STRAND == auto ]]; then
  STRAND=$(awk 'FNR>4{f+=$3; r+=$4} END{x=(r>0?f/r:99); print (x>3 ? 1 : (x<1/3 ? 2 : 0))}' "$OUT"/star/*.ReadsPerGene.out.tab)
  awk 'FNR>4{f+=$3; r+=$4} END{printf "forward/reverse ratio = %.2f\n", (r>0?f/r:99)}' "$OUT"/star/*.ReadsPerGene.out.tab \
      > "$OUT/logs/strandedness.txt"
  echo "featureCounts -s $STRAND" >> "$OUT/logs/strandedness.txt"
fi
echo "Strandedness: -s $STRAND"

# 7d. featureCounts: gene-level counts for all samples in one matrix
#     paired-end: count fragments (-p --countReadPairs), both mates mapped (-B), no chimeric pairs (-C)
mkdir -p "$OUT/featurecounts"
PE_FLAGS=""
if [[ -n "$(samtools view -f 1 "${BAMS[0]}" 2>/dev/null | head -1 || true)" ]]; then
  PE_FLAGS="-p --countReadPairs -B -C"
fi
featureCounts -T "$THREADS" -a "$GTF" -g gene_name -s "$STRAND" $PE_FLAGS \
              -o "$OUT/featurecounts/counts.txt" "${BAMS[@]}"
# clean matrix for DESeq2/edgeR: gene + one column per SRR
awk 'BEGIN{OFS="\t"} /^#/ {next}
     $1=="Geneid" { printf "gene"; for(i=7;i<=NF;i++){ n=$i; sub(/.*\//,"",n); sub(/\.sorted\.bam$/,"",n); printf "\t%s", n }; print ""; next }
     { printf "%s", $1; for(i=7;i<=NF;i++) printf "\t%s", $i; print "" }' \
    "$OUT/featurecounts/counts.txt" > "$OUT/featurecounts/count_matrix.tsv"

# 8. FastQC-only report: raw vs trimmed reads 
# multiqc "$OUT/fastqc_raw" "$OUT/fastqc_trimmed" -o "$OUT/multiqc_fastqc" -f -d -dd 1

# 9. BAMs containing only the genes in the gene list
mkdir -p "$OUT/genelist_bam"
for b in "${BAMS[@]}"; do
  s=$(basename "$b" .sorted.bam)
  samtools view -b -M -L "$GENE_BED" "$b" > "$OUT/genelist_bam/${s}.genes.bam"
  samtools index "$OUT/genelist_bam/${s}.genes.bam"
done

echo "All done. Gene-list BAMs: $OUT/genelist_bam/"
