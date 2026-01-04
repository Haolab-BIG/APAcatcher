
# APAcatcher
## Introduction

**Alternative polyadenylation (APA)**, a key post-transcriptional regulator affecting > 70% of human genes,
generates mRNA isoforms with distinct 3’ ends. Despite its widespread biological importance, APA detection remains
challenging in full-length bulk and single-cell RNA-seq, especially for datasets derived from rare specimens.
Here, we developed **APAcatcher**, a deep learning embedded approach integrating both **in vivo** cleavage patterns
from sequencing data and **in vitro** DNA sequence features, achieving accurate and robust APA detection and quantification. 

## Supported Species

- **Human**
- **Mouse**
- **Drosophila**
- **Arabidopsis**
- **Yeast**
- **Zebrafish**
- **Caenorhabditis elegans (C.elegans)**
## System Requirements
- Operating System: Linux / macOS 
- Hardware: No non-standard hardware required (runs on a standard laptop/desktop)
---

## Installation (Singularity Version)
You can download Singularity on https://docs.sylabs.io/guides/3.2/user-guide/installation.html
## 1. Download APAcatcher Singularity Image

Pull the pre-built image directly from the Singularity Library:

```bash
singularity pull --arch amd64 library://c10969/default/apacatcher:latest
```
## Notes on Running Commands

All commands are identical to the Conda version **except**:

- Each command must be prefixed with  
  ```bash
  singularity exec --bind /path/to/data:/data apacatcher.sif
  ```
- Always use **absolute paths**

Example:
```bash
singularity exec --bind /path/to/data:/data apacatcher.sif main.py --help
```

---
## Installation (Conda Version)
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
> **Note:** The script has been tested with Python 3.8.18, bedtools(v2.31.0), salmon(v1.10.1).

---
## Preparation: Bam to depth file

APAcatcher using depth file as input

**💡 Recommended Genome Annotations**
To achieve the best results, we recommend using the following GTF versions to generate your BAM files for each species:

| Species          | Recommended GTF File                         |
|------------------|----------------------------------------------|
|Human             | gencode.v46.primary_assembly.annotation.gtf  |
| Mouse            | mm10.ncbiRefSeq.gtf                          |
| Zebrafish        | Danio_rerio.GRCz11.115.gtf                   |
| Drosophila       | Drosophila_melanogaster.BDGP6.54.115.gtf     |
| Arabidopsis      | Arabidopsis_thaliana.TAIR10.62.gtf           |
| Yeast            | Saccharomyces_cerevisiae.R64-1-1.115.gtf     |
| C. elegans       | Caenorhabditis_elegans.WBcel235.115.gtf      |


### Singularity version
```
# extract reads in 3’UTR regions
singularity exec --bind /path/to/data:/data apacatcher.sif \
  samtools view -hb -L /data/UTR_Human.bed /data/sample1.bam > sample1_3UTR.bam

# convert bam to depth.txt
singularity exec --bind /path/to/data:/data apacatcher.sif \
  samtools depth sample1_3UTR.bam -b /data/UTR_Human.bed > sample1_3UTR_read_coverage.txt
```

### Conda version
```
# extract reads in 3’UTR regions
samtools view -hb -L UTR_Human.bed sample1.bam > sample1_3UTR.bam


# convert bam to depth.txt
samtools depth sample1_3UTR.bam -b UTR_Human.bed > sample1_3UTR_read_coverage.txt
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
| Generate input files                 | ~5 min per sample            |
| Identify high-confidence APA sites   | ~10 min per samlpe           |
| Quantify with Salmon                 | ~10-20 min per samlpe        |
| Assemble the final data matrix       | ~1-2 min                     |

Runtime estimates vary depending on dataset size, sequencing depth, and available computing resources.

---
**Notes**
⚠️ Please run the script using absolute paths 

### 1. Generate Input Files

#### Original depth file format

```
chr1    70009   0
chr1    70010   0
chr1    70011   0
…
```

First, annotate each depth file with gene%%transcript_number and strand information (based on RefSeq_UTR_final.bed, please see the detailed information from the README.md in the “demo” file):
#### Singularity version
```
singularity exec --bind /path/to/data:/data apacatcher.sif add_geneinfo.py \
-g /data/UTR_Human.bed \
-d /data/input_depth_file_dir
```
#### Conda version
```bash
python add_geneinfo.py \
  -g UTR_Human.bed \
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
🚨 Runtime Warning
PAS detection is the most time-consuming step in the pipeline.
Depending on the average 3′ UTR length, processing each sample may require 30 minutes.(In Human H1 cellline,40M reads,3 replicates)


Options for `main.py`:

| Flag                  | Description                        | Default |
| --------------------  | ---------------------------------- | ------- |
| `--input_folder`      | Path to annotated depth files      |         |
| `--genome_file`       | Reference genome FASTA             |         |
| `--output_folder`     | Directory to save predicted PAS    |         |
| `--tpm_threshold`     | Minimum TPM to consider            | `1`     |
| `--length_threshold`  | Minimum 3' UTR length (bp)         | `100`   |
| `--penalty`           | Penalty for change-point detection | `50`    |
| `--min_size`          | Minimum segment size for PELT      | `30`    |
| `--num_processes`     | Number of parallel workers         | `4`     |
| `--model_path`        | Path to the pre-trained model for different species. The species (e.g., Human) is automatically inferred from the filename.         | `APAcatcher_Human.pth`     |

**Example:**
#### Singularity version
```
singularity exec --bind /path/to/data:/data apacatcher.sif \
  main.py \
  --input_folder /data/depth_file_dir \
  --genome_file /data/hg38.fa \
  --model_path /data/muti_species/APAcatcher.Human \
  --num_processes 8 \
  --output_folder /data/high_confidence_pas_folder
```

#### Conda version
```bash
python main.py \
  --input_folder depth_file_dir \
  --genome_file hg38.fa \
  --model_path ./muti_species/APAcatcher.Human \
  --num_processes 8 \
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
#### Singularity version
```bash
singularity exec --bind /path/to/data:/data apacatcher.sif \
  cluster_bed_files.sh \
  -i /data/high_confidence_pas_folder \
  -o /data/cluster_high_confidence_pas_folder
```


#### Conda version
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
#### Singularity version
```
singularity exec --bind /path/to/data:/data apacatcher.sif \
  combind.sh \
  -i /data/cluster_high_confidence_pas_folder \
  -o /data/cluster_high_confidence_pas_folder/combind
```

#### Conda version
```bash
./combind.sh \
  -i cluster_high_confidence_pas_folder \
  -o cluster_high_confidence_pas_folder/combind
```

---

#### 2.4 Remove annotated-proximal sites

Options for `process_last.py`:
| Flag                | Description                              |
| ------------------- | -----------------------------------------|
| `--input`      | single group or mulitiple group BED file |
| `--output`     | Output directory                         |

#### Singularity version
```
singularity exec --bind /path/to/data:/data apacatcher.sif \
  process_last.py \
  --input /data/<PAS_BED> \
  --output /data/final_site_for_quantification.bed
```


#### Conda version
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
# chr   start     end       cluster_size  gene%%transcript_number  strand
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
| `-r` | UTR\_Human.bed                                            |
| `-l` | UTR\_lastexon\_Human.bed                                  |
| `-g` | UCSC\_hg38.fa                                             |
| `-o` | Output FASTA for 3′ UTR isoforms (e.g. `3UTRisoforms.fa`) |
| `-i` | Salmon index directory (e.g. `3UTRisoforms_library`)      |

**Example:**

#### Singularity version
```bash
singularity exec --bind /path/to/data:/data apacatcher.sif \
  get_salmon_index.sh \
  -f /data/final_site.bed \
  -r /data/UTR_Human.bed \
  -l /data/UTR_lastExon_Human.bed \
  -g /data/hg38.fa \
  -o /data/3UTRisoforms.fa \
  -i /data/3UTRisoforms_library
```

#### Conda version
```bash
./get_salmon_index.sh \
  -f final_site_for_quantification.bed \
  -r UTR_Human.bed \
  -l UTR_lastExon_Human.bed \
  -g hg38.fa \
  -o 3UTRisoforms.fa \
  -i 3UTRisoforms_library
```

---

#### 3.2 Quantification
⚠️ **`FASTQ Naming Requirement`** 
Input FASTQ files must follow a paired-end naming convention so that the script can correctly match read pairs.
Supported formats include, but are not limited to:
- **`sample1_R1.fastq.gz and sample1_R2.fastq.gz`** 
- **`sample1_1.fastq.gz and sample1_2.fastq.gz`** 


Ensure that read pairs share the same sample prefix and differ only in the read identifier (R1/R2 or 1/2).
Options for `get_quant.sh`:

| Flag | Description                                 |
| ---- | ------------------------------------------- |
| `-i` | Salmon index directory                      |
| `-d` | Directory of input FASTQ files              |
| `-o` | Output directory for quantification results |


#### Singularity version
**Example for paired-end:**
```bash
singularity exec --bind /path/to/data:/data apacatcher.sif \
  get_quant.sh \
  -i /data/3UTRisoforms_library \
  -d /data/fastq_directory \
  -o /data/quant_results
```

**Example for single-end:**
```bash
  singularity exec --bind /path/to/data:/data apacatcher.sif \
  get_quant_single_end.sh \
  -i /data/3UTRisoforms_library \
  -d /data/fastq_directory \
  -o /data/quant_results
```

#### Conda version
**Example for paired-end:**
```bash
./get_quant.sh \
  -i 3UTRisoforms_library \
  -d fastq_directory \
  -o quant_results
```
**Example for single-end:**
```bash
  ./get_quant_single_end.sh \
  -i 3UTRisoforms_library \
  -d fastq_directory \
  -o quant_results
```
---

#### 3.3 Merge Quantification Across Samples

#### Singularity version
```bash
  singularity exec --bind /path/to/data:/data apacatcher.sif \
  merge_quant.sh \
  -l /data/3UTRisoforms_library \
  -b /data/quant_results \
  -o /data/merged_tpm.txt
```

#### Conda version
```bash
./merge_quant.sh \
  -l sample_list.txt \
  -b /path/to/quant_results \
  -o merged_tpm.txt
```

* `sample_list.txt`: one sample directory per line  
* `-b`: base directory containing each sample’s quantification results  

📝 **Sample list format**  
The `sample_list.txt` file should contain **one sample name per line**, corresponding to the sample prefix used in FASTQ file names and quantification output directories.

For example, if your FASTQ files are named `sample1_R1.fastq.gz` / `sample1_R2.fastq.gz` (or `sample1_1.fastq.gz` / `sample1_2.fastq.gz`), the sample name should be:

```text
sample1
sample2
sample3
```

---

### 4. Assemble Final Data Matrix

Options for `get_final_result.py`:

| Flag            | Description                                   |
| --------------- | --------------------------------------------- |
| `--group_files` | Paths to sample group definition files        |
| `--merge_file`  | Path to `merged_tpm.txt`                      |
| `--output_dir`  | Directory to save the final matrices          |
| `--length`      | Minimum 3′ UTR length to include              |

🧩 **Group file format**  
Each group file (e.g., `group_A.txt`, `group_B.txt`) should contain **one sample name per line**, matching the sample names defined in `sample_list.txt`.

- `group_A.txt`: treatment / KO / case group samples  
- `group_B.txt`: control / WT / reference group samples  

Example:

```text
# group_A.txt (e.g., KO or treatment group)
sample2
sample4
sample6
```

```text
# group_B.txt (e.g., WT or control group)
sample1
sample3
sample5
```

**Example:**
#### Singularity version
```bash
  singularity exec --bind /path/to/data:/data apacatcher.sif \
  get_final_result.py \
  --group_files /data/group_A.txt /data/group_B.txt \
  --merge_file /data/merged_tpm.txt \
  --output_dir /data/final_result \
  --length 100
```


#### Conda version
```bash
python get_final_result.py \
  --group_files group_A.txt group_B.txt \
  --merge_file merged_tpm.txt \
  --output_dir final_result \
  --length 100
```

---

#### Final 3′ UTR index example

```text
# Name                                    Length  Transcript   start       end         strand    sample1_indexUTR  sample2_indexUTR  sample3_indexUTR
AACS%%3:+::chr12:125142091-125142380      289     AACS%%3      125142091   125142380   +         0.6331            0.8141            0.6945
…
```

---

#### Output files

The `get_final_result.py` script generates the following output files in the specified output directory:

- **`3UTR_averageLength.txt`**  
  Average 3′ UTR length for each gene across samples, calculated based on isoform usage.

- **`3UTR_usage.txt`**  
  Relative usage of alternative 3′ UTR isoforms for each gene in each sample.

- **`3UTR_index.txt`**  
  Sample-wise 3′ UTR index values, representing the weighted average 3′ UTR length and serving as the primary input for downstream APA and differential analysis.

- **`PDUI.txt`**  
  Percentage of distal poly(A) site usage (PDUI) for each gene across samples.

- **`PPUI.txt`**  
  Percentage of proximal poly(A) site usage (PPUI), complementary to PDUI.

- **`TPM.txt`**  
  Transcript-level or isoform-level TPM matrix aggregated across all samples, used for downstream differential expression or APA-related analyses.

All output files are tab-delimited.

⚠️ **Data Preprocessing Note**:

In the generated .txt files, any original values recorded as **0** (which represent null or missing values)
have been automatically converted to **NA**(Display as empty in output files).This conversion is implemented to ensure data integrity
and compatibility for subsequent statistical testing and differential expression analysis.


### 5. DE analysis
Options for `DE_analysis.py`:
| Flag | Description                                            |
| ---- | ------------------------------------------------------ |
| `-t` | Path to the TPM file (e.g., `TPM.txt`)                 |
| `-i` | Path to the 3′ UTR index file (e.g., `3UTR_index.txt`) |
| `-a` | Sample list for group A (e.g., KO group)               |
| `-b` | Sample list for group B (e.g., WT group)               |
| `-m` | Mode for DE analysis (e.g., `index`)                   |
| `-o` | Output file name                                       |

**Example:**
#### Singularity version
```bash
  singularity exec --bind /path/to/data:/data apacatcher.sif \
  DE_analysis.py \
  -t /data/final_result/TPM.txt \
  -i /data/final_result/3UTR_index.txt \
  -a /data/groupA.txt \
  -b /data/groupB.txt \
  -m index \
  -o /data/KO_WT_index.txt
```


#### Conda version
```bash
python DE_analysis.py \
  -t /final_result/TPM.txt \
  -i /final_result/3UTR_index.txt \
  -a groupA.txt \
  -b groupB.txt \
  -m index \
  -o KO_WT_index.txt
```
**Notes**
---
Ensure that sample names in Group_A and Group_B files are consistent with the column names in the TPM file.

---

*All scripts and parameters are fully customizable. For questions or issues, please open an issue on the project repository.*
