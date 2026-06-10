
# scRNA-seq Pipeline: FASTQ to Depth Files

This pipeline processes SMART-seq data (with UMI) from raw FASTQ files to
per-cell 3′ UTR depth files suitable for APA analysis with APAcatcher/WAPI-KDE.

---

## Dependencies

| Tool | Version tested | Install |
|------|---------------|---------|
| Trim Galore | ≥ 0.6.7 | `conda install -c bioconda trim-galore` |
| STAR | ≥ 2.7.10a | `conda install -c bioconda star` |
| UMI-tools | ≥ 1.1.2 | `conda install -c bioconda umi_tools` |
| samtools | ≥ 1.15 | `conda install -c bioconda samtools` |
| featureCounts | ≥ 2.0.3 | `conda install -c bioconda subread` |

---

## Input

```
data/
  ├── cell_001_R1.fastq.gz      # Read 1 (contains UMI)
  ├── cell_001_R2.fastq.gz      # Read 2 (cDNA read)
  ├── cell_002_R1.fastq.gz
  ├── cell_002_R2.fastq.gz
  └── ...
```

Each cell is a separate pair of FASTQ files (one library per cell).

---

## Reference preparation (run once)

```bash
# Download human reference genome and annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v38.annotation.gtf.gz

# Build STAR index
STAR \
    --runMode genomeGenerate \
    --genomeDir ref/star_index \
    --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile gencode.v38.annotation.gtf \
    --runThreadN 8
```

---

## Pipeline (per cell)

Replace `CELL` with the cell identifier (e.g., `cell_001`).

### Step 1: Adapter trimming

```bash
CELL=cell_001

trim_galore \
    --paired \
    --cores 4 \
    --output_dir trimmed/ \
    data/${CELL}_R1.fastq.gz \
    data/${CELL}_R2.fastq.gz
```

Output: `trimmed/${CELL}_R1_val_1.fq.gz`, `trimmed/${CELL}_R2_val_2.fq.gz`

---

### Step 2: Extract UMI

UMI is located at the 5′ end of Read 1. Adjust `--bc-pattern` to match your
library structure (e.g., `NNNNNNNN` for 8-bp UMI).

```bash
umi_tools extract \
    --bc-pattern=NNNNNNNN \
    --stdin  trimmed/${CELL}_R1_val_1.fq.gz \
    --stdout trimmed/${CELL}_R1_umi.fq.gz \
    --read2-in  trimmed/${CELL}_R2_val_2.fq.gz \
    --read2-out trimmed/${CELL}_R2_umi.fq.gz \
    --log logs/${CELL}_umi_extract.log
```

Output: `trimmed/${CELL}_R1_umi.fq.gz`, `trimmed/${CELL}_R2_umi.fq.gz`

---

### Step 3: Align to reference genome

```bash
STAR \
    --runThreadN 8 \
    --genomeDir ref/star_index \
    --readFilesIn \
        trimmed/${CELL}_R2_umi.fq.gz \
        trimmed/${CELL}_R1_umi.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outFileNamePrefix bam/${CELL}_ \
    --outSAMunmapped Within \
    --runMode alignReads

# Index BAM
samtools index bam/${CELL}_Aligned.sortedByCoord.out.bam
```

Output: `bam/${CELL}_Aligned.sortedByCoord.out.bam`

---

### Step 4: Deduplicate by UMI

```bash
umi_tools dedup \
    --stdin  bam/${CELL}_Aligned.sortedByCoord.out.bam \
    --stdout bam/${CELL}_dedup.bam \
    --log    logs/${CELL}_umi_dedup.log \
    --paired

samtools index bam/${CELL}_dedup.bam
```

Output: `bam/${CELL}_dedup.bam`

---

### Step 5: Generate 3′ UTR depth file

Extract per-base coverage within annotated 3′ UTR regions.
`UTR_Human.bed` is the UTR annotation file used by APAcatcher.

```bash
# Extract coverage over 3' UTR regions only
samtools view -b -L UTR_Human.bed bam/${CELL}_dedup.bam \
    | samtools depth -a -b UTR_Human.bed - \
    | awk -v cell=${CELL} 'BEGIN{OFS="\t"} {print $1, $2, $3, cell, "."}' \
    > depth/${CELL}_UTR_depth.txt
```

Output format (tab-delimited):
```
chr1    1001    5    cell_001    .
chr1    1002    8    cell_001    .
chr1    1003    3    cell_001    .
```

Columns: `chromosome`, `position`, `depth`, `cell_id`, `strand`

> **Note**: strand information should be added based on the UTR annotation.

---
---
