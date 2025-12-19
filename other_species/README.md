# other_species

This directory contains species-specific resources for running APAcatcher, including
pre-trained deep learning models (PyTorch) and UTR annotation files required for RNA-seq
quantification and downstream APA analysis.

Supported species:
- Arabidopsis
- Caenorhabditis elegans (Celegans)
- Drosophila
- Mouse
- Yeast
- Zebrafish


Directory Structure

other_species/
├── Arabidopsis/
│   ├── APAcatcher_Arabidopsis.pth
│   ├── UTR_Arabidopsis.bed
│   └── UTR_lastExon_Arabidopsis.bed
├── Celegans/
│   ├── APAcatcher_Celegans.pth
│   ├── UTR_Celegans.bed
│   └── UTR_lastExon_Celegans.bed
├── Drosophila/
├── Mouse/
├── Yeast/
└── Zebrafish/


File Descriptions

1. Deep Learning Model (*.pth)

File name format:
APAcatcher_<Species>.pth

Example:
APAcatcher_Celegans.pth

Description:
Species-specific pre-trained deep learning model implemented in PyTorch.

Usage:
This file is required when running main.py and should be specified via the
--model_path argument.

Example command:
python main.py \
  --model_path other_species/Celegans/APAcatcher_Celegans.pth \
  [other arguments]

Important:
Each species must use its corresponding model parameters.
Models are not interchangeable across species.


2. UTR_<Species>.bed

Example:
UTR_Arabidopsis.bed

Description:
BED file containing UTR region annotations for the corresponding species.

Purpose:
- Used when converting BAM files to read depth text files via samtools
- Required for extracting coverage information in UTR regions

Typical usage:
samtools depth -b UTR_Arabidopsis.bed aligned.bam > utr_depth.txt


3. UTR_lastExon_<Species>.bed

Example:
UTR_lastExon_Arabidopsis.bed

Description:
BED file containing annotations of last exon UTR regions.

Purpose:
- Used during Salmon quantification
- Ensures correct transcript-level abundance estimation for APA analysis


Workflow Summary

For a given species, the typical workflow is:

1. RNA-seq quantification
   Use Salmon with the provided UTR-related BED files

2. Read depth extraction
   Use samtools depth with UTR_<Species>.bed

3. APAcatcher prediction
   Run main.py with the corresponding APAcatcher_<Species>.pth model


Notes

- All models are trained independently per species
- BED files are species-specific and required for correct quantification
- Mixing models or BED files across species may lead to incorrect results
