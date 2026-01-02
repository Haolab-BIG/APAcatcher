## Preparation: Download FASTQ Files

Before generating BAM and depth files, raw sequencing data should be downloaded in FASTQ format.

We provide a shell script (`download.sh`) containing all download commands. Please make sure the script has execution permission and run it as follows:

```bash
chmod +x download.sh
./download.sh

## Preparation: Bam to Depth File

```bash
cd /Sample_data/bam
for bam in *.bam; do
    base=$(basename "$bam" _sorted.bam)

    # Extract reads in 3’UTR regions
    samtools view -hb -L UTR_Human.bed "$bam" > "${base}_3UTR.bam"

    # Convert BAM to depth.txt
    samtools depth "${base}_3UTR.bam" -b RefSeq_UTR_final.bed > "${base}_3UTR_read_coverage.txt"
done
```

## 1. Generate Input Files

```bash
python add_geneinfo.py -g ./UTR_Human.bed -d ./Sample_data/bam
```

## 2. Identification of High-Confidence APA Sites

### 2.1 PAS Detection Using PELT and Deep Learning

```bash
python main.py \
  --input_folder ./Sample_data/bam \
  --genome_file ~/hg38.fa \
  --num_processes 16 \
  --tpm_threshold 1 \
  --length_threshold 100 \
  --penalty 50 \
  --min_size 30 \
  --model_path ./APAcatcher_Human.pth \
  --output_folder ./Sample_output/high_confidence_pas_folder
```

### 2.2 Intra-Group Clustering

> Note: Samples must be organized into separate folders by experimental group.  
> SRR5076687 belong to pre_IGF

> SRR5076712 belong to post_IGF

> SRR5076688 belong to pre_DGF

> SRR5076713 belong to post_DGF

```bash
for group in pre_IGF post_IGF pre_DGF post_DGF; do
    ./cluster_bed_files.sh \
        -i ./Sample_output/high_confidence_pas_folder/$group \
        -o ./Sample_output/cluster_high_confidence_pas_folder
done
```

### 2.3 Inter-Group Merging

```bash
./combind.sh
-i absolute_path/Sample_output/cluster_high_confidence_pas_folder
-o absolute_path/Sample_output/cluster_high_confidence_pas_folder/combind
```

### 2.4 Removal of Annotated-Proximal Sites

```bash
python process_last.py \
  --input ./Sample_output/cluster_high_confidence_pas_folder/combind/combind_pas_site.bed \
  --output ./Sample_output/final_site_for_quantification.bed
```

## 3. Quantification with Salmon

### 3.1 Build Salmon Index

```bash
./get_salmon_index.sh \
  -f ./Sample_output/final_site_for_quantification.bed \
  -r ./UTR_Human.bed \
  -l ./UTR_lastExon_Human.bed \
  -g hg38.fa \
  -o ./Sample_output/3UTRisoforms.fa \
  -i ./Sample_output/3UTRisoforms_library
```

### 3.2 Quantification

```bash
./get_quant.sh \
  -i ./Sample_output/3UTRisoforms_library \
  -d ./Sample_data/fastq \
  -o ./Sample_output/quant_results
```

### 3.3 Merge Quantification Across Samples

```bash
./merge_quant.sh \
  -l sample_list.txt \
  -b ./Sample_output/quant_results \
  -o ./Sample_output/merged_tpm.txt
```

## 4. Assemble Final Data Matrix

```bash
python get_final_result.py \
  --group_files pre_IGF.txt post_IGF.txt pre_DGF.txt post_DGF.txt \
  --merge_file ./Sample_output/merged_tpm.txt \
  --output_dir ./Sample_output/final_result \
  --length 100
```

### 5. DE analysis
```bash
python DE_analysis.py \
  -t /final_result/TPM.txt \
  -i /final_result/3UTR_index.txt \
  -a post_IGF.txt \
  -b pre_IGF.txt \
  -m index \
  -o post_pre_IGF_index.txt



python DE_analysis.py \
  -t /final_result/TPM.txt \
  -i /final_result/3UTR_index.txt \
  -a post_DGF.txt \
  -b pre_DGF.txt \
  -m index \
  -o post_pre_DGF_index.txt 
```
