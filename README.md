
# APAcatcher
## System Requirements
- Operating System: Linux / macOS / Windows
- Hardware: No non-standard hardware required (runs on a standard laptop/desktop)

## Installation
Typical install time on a "normal" desktop computer: about 10 minutes

We recommend creating a dedicated Conda environment:

```bash
conda create -n apacatcher_env python=3.8
conda activate apacatcher_env
```

Then install the required Python packages:

```bash
pip install \
  numpy==1.24.3 \
  pandas==2.0.3 \
  scipy==1.10.1 \
  torch==2.2.1 \
  ruptures==1.1.9 \
  biopython==1.83 \
  einops==0.8.0 \
  tqdm==4.66.5
```

> **Note:** The script has been tested with Python 3.8.18.

---
## Preparation: Bam to depth file
APAcatcher using depth file as input
```
# extract reads in 3’UTR regions
samtools view -hb -L RefSeq_UTR_final.bed sample1.bam > sample1_3UTR.bam


# convert bam to depth.txt
samtools depth sample1_3UTR.bam -b RefSeq_UTR_final.bed > sample1_3UTR_read_coverage.txt
```

## Usage Overview

APAcatcher operates in four main stages:

1. **Generate input files**
2. **Identify high-confidence APA sites**
3. **Quantify with Salmon**
4. **Assemble the final data matrix**

Estimated Runtime per Stage (~30M reads, 8 kernel)

| Stage                                | Typical Runtime*             |
|--------------------------------------|------------------------------|
| Generate input files                 | ~5 min                       |
| Identify high-confidence APA sites   | ~12 min                      |
| Quantify with Salmon                 | ~20-30 min                   |
| Assemble the final data matrix       | ~5-10 min                    |

Runtime estimates vary depending on dataset size, sequencing depth, and available computing resources.

---

### 1. Generate Input Files

#### Original depth file format

```
chr1    70009   0
chr1    70010   0
chr1    70011   0
…
```

First, annotate each depth file with gene%%transcript and strand information (based on `RefSeq_UTR_final.bed`):

```bash
python add_geneinfo.py \
  -g RefSeq_UTR_final.bed \
  -d input_depth_file_dir
```

#### Annotated depth file example

```
chr1    70009   0   OR4F5%%1    +
chr1    70010   0   OR4F5%%1    +
chr1    70011   0   OR4F5%%1    +
…
```

---

### 2. Identify High-Confidence APA Sites

#### 2.1 PAS detection with PELT + Deep Learning

Options for `main.py`:

| Flag                 | Description                        | Default |
| -------------------- | ---------------------------------- | ------- |
| `--input_folder`     | Path to annotated depth files      |         |
| `--genome_file`      | Reference genome FASTA             |         |
| `--output_folder`    | Directory to save predicted PAS    |         |
| `--tpm_threshold`    | Minimum TPM to consider            | `1`     |
| `--length_threshold` | Minimum 3' UTR length (bp)         | `100`   |
| `--penalty`          | Penalty for change-point detection | `50`    |
| `--min_size`         | Minimum segment size for PELT      | `30`    |
| `--num_processes`    | Number of parallel workers         | `4`     |

**Example:**

```bash
python main.py \
  --input_folder depth_file_dir \
  --genome_file hg38.fa \
  --output_folder high_confidence_pas_folder
```

---

#### 2.2 Intra-group clustering

Options for `cluster_bed_files.sh`:

| Flag | Description                         | Default |
| ---- | ----------------------------------- | ------- |
| `-i` | Directory with per-sample BED files |         |
| `-o` | Output directory                    |         |
| `-d` | Merge distance (bp)                 | `70`    |
| `-c` | Minimum replicate count             | `2`     |

**Example:**

```bash
./cluster_bed_files.sh \
  -i high_confidence_pas_folder \
  -o cluster_high_confidence_pas_folder
```

---

#### 2.3 Inter-group merging
if you only have one group，pass this command

Options for `combind.sh`:

| Flag                | Description                    | Default |
| ------------------- | ------------------------------ | ------- |
| `-i`                | Directory with group BED files |         |
| `-o`                | Output directory               |         |
| `-d MERGE_DISTANCE` | Merge distance (bp)            | `70`    |

**Example:**

```bash
./combind.sh \
  -i cluster_high_confidence_pas_folder \
  -o cluster_high_confidence_pas_folder/combind
```

---

#### 2.4 Remove annotated-proximal sites

Options for `process_last.py`:
| Flag                | Description                              | Default |
| ------------------- | -----------------------------------------| ------- |
| `--input`      | single group or mulitiple group BED file |         |
| `--output`     | Output directory                         |         |
```bash
python process_last.py \
  --input <PAS_BED> \
  --output final_site_for_quantification.bed
```

* Use `pas_site.bed` for a single group.
* Use `combind_pas_site.bed` for multiple groups.

---

**High-confidence APA sites (BED) example:**

```
# chr   start     end       cluster_size  gene%%transcript  strand
chr1   944201    944201    4             NOC2L%%1          -
chr1   965721    965721    2             KLHL17%%1         +
…
```

---

### 3. Quantify with Salmon

#### 3.1 Build Salmon index

Options for `get_salmon_index.sh`:

| Flag | Description                                               |
| ---- | --------------------------------------------------------- |
| `-f` | final\_site\_for\_quantification.bed                      |
| `-r` | RefSeq\_UTR\_final.bed                                    |
| `-l` | RefSeq\_UTR\_lastexon\_final.bed                          |
| `-g` | UCSC\_hg38.fa                                             |
| `-o` | Output FASTA for 3′ UTR isoforms (e.g. `3UTRisoforms.fa`) |
| `-i` | Salmon index directory (e.g. `3UTRisoforms_library`)      |

**Example:**

```bash
./get_salmon_index.sh \
  -f final_site_for_quantification.bed \
  -r RefSeq_UTR_final.bed \
  -l RefSeq_UTR_lastexon_final.bed \
  -g hg38.fa \
  -o 3UTRisoforms.fa \
  -i 3UTRisoforms_library
```

---

#### 3.2 Quantification

Options for `get_quant.sh`:

| Flag | Description                                 |
| ---- | ------------------------------------------- |
| `-i` | Salmon index directory                      |
| `-d` | Directory of input FASTQ files              |
| `-o` | Output directory for quantification results |

**Example:**

```bash
./get_quant.sh \
  -i 3UTRisoforms_library \
  -d fastq_directory \
  -o quant_results
```

---

#### 3.3 Merge quantification across samples

Options for `merge_quant.sh`:

```bash
merge_quant.sh \
  -l sample_list.txt \
  -b /path/to/quant_results \
  -o merged_tpm.txt
```

* `sample_list.txt`: one sample directory per line
* `-b`: base directory containing each sample’s results

---

### 4. Assemble Final Data Matrix

Options for `get_final_result.py`:

| Flag            | Description                          |
| --------------- | ------------------------------------ |
| `--group_files` | Paths to sample-group text files     |
| `--merge_file`  | Path to `merged_tpm.txt`             |
| `--output_dir`  | Directory to save the final matrices |
| `--length`      | Minimum 3′ UTR length to include     |

**Example:**

```bash
python get_final_result.py \
  --group_files group_A.txt group_B.txt \
  --merge_file merged_tpm.txt \
  --output_dir final_result \
  --length 100
```

---

#### Final 3UTR index example

```
# Name                                    Length  Transcript   start       end         strand    sample1_indexUTR  sample2_indexUTR  sample3_indexUTR
AACS%%3:+::chr12:125142091-125142380      289     AACS%%3      125142091   125142380   +         0.6331            0.8141            0.6945
…
```

---

*All scripts and parameters are fully customizable. For questions or issues, please open an issue on the project repository.*
