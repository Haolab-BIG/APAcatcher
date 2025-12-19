# other_species

This directory contains **species-specific resources** for running **APAcatcher**, including  
pre-trained deep learning models (PyTorch) and UTR annotation files required for RNA-seq
quantification and downstream APA analysis.

## Supported Species

- **Arabidopsis**
- **Caenorhabditis elegans (Celegans)**
- **Drosophila**
- **Mouse**
- **Yeast**
- **Zebrafish**

---

## Directory Structure

Each species is organized in an independent subdirectory with the following structure:

```text
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
