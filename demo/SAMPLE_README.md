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
> SRR5076687 and SRR5076688 belong to Group_Pre  
> SRR5076712 and SRR5076713 belong to Group_Post

```bash
./cluster_bed_files.sh -i ./Sample_output/high_confidence_pas_folder/Group_Pre -o ./Sample_output/cluster_high_confidence_pas_folder
./cluster_bed_files.sh -i ./Sample_output/high_confidence_pas_folder/Group_Post -o ./Sample_output/cluster_high_confidence_pas_folder
```

### 2.3 Inter-Group Merging

```bash
./combind.sh -i ./Sample_output/cluster_high_confidence_pas_folder -o ./Sample_output/cluster_high_confidence_pas_folder/combind
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
