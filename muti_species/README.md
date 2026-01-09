# muti_species

This directory contains **species-specific resources** for running **APAcatcher**, including  
pre-trained deep learning models (PyTorch) and UTR annotation files required for RNA-seq
quantification and downstream APA analysis.

## Supported Species

- *Arabidopsis*
- *Caenorhabditis elegans (C.elegans)*
- *Drosophila*
- **Mouse**
- **Yeast**
- **Zebrafish**
- **Human**
---

## Directory Structure

Each species is organized in an independent subdirectory with the following structure:

```text
other_species/
├── Arabidopsis/
│   ├── APAcatcher_Arabidopsis_301.pth
│   ├── UTR_Arabidopsis.bed
│   └── UTR_lastExon_Arabidopsis.bed
├── C.elegans/
│   ├── APAcatcher_Celegans.pth
│   ├── UTR_Celegans.bed
│   └── UTR_lastExon_Celegans.bed
├── Drosophila/
├── Mouse/
├── Yeast/
└── Zebrafish/
└── Human/
```


### File Descriptions

Each species-specific folder contains three files:

- **`APAcatcher_[Species].pth`**  
  Pre-trained PyTorch model checkpoint (e.g., `APAcatcher_Celegans.pth`).  
  This file is loaded via the `--model_path` argument when running the main inference script (typically `main.py`).  
  **Important**: Each species uses an independently trained model optimized for its sequence characteristics.

- **`UTR_[Species].bed`**  
  BED file defining 3' UTR regions for the species.  
  Primary use:  
  - Filtering or computing read depth when converting BAM files to coverage txt using tools like `samtools depth` or custom scripts.

- **`UTR_lastExon_[Species].bed`**  
  BED file defining the last exon / terminal UTR regions.  
  Primary use:  
  - Required for Salmon quantification with custom indexing of APA isoforms or terminal exons.

## Usage

### Loading the Model
When running APAcatcher inference:

```bash
python main.py --model_path muti_species/C.elegans/APAcatcher_Celegans.pth [other arguments...]
```
Replace the path with the appropriate species folder and model file.


## Using BED Files
For read depth extraction from BAM files (e.g., with samtools):

```bash
samtools view -hb -L UTR_[Species].bed sample1.bam > sample_3UTR.bam
samtools depth sample_3UTR.bam -b UTR_[Species].bed > sample_3UTR_read_coverage.txt
```

For Salmon index involving terminal exons:
```bash
./get_salmon_index.sh \
  -f final_site_for_quantification.bed \
  -r UTR_[Species].bed \
  -l UTR_lastExon_[Species].bed \
  -g species.fa \
  -o 3UTRisoforms.fa \
  -i 3UTRisoforms_library
```
