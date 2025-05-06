# APAcatcher
## Install APAcatcher
## run APAcatcher in 4 step
### 1.generate input files
```bash
python add_geneinfo.py -g RefSeq_UTR_final.bed -d depth_file_dir -p 8
```
### 2.get pesudo High confidence APA sites
```bash
#2.1 Using PELT and DL model get PloyA sites
python main.py --input_folder depth_file_dir --genome_file hg38.fa --output_folder high_confidence_pas_folder --tpm_threshold 1 --length_threshold 100  --penalty 50 --min_size 30 --num_processes 8

#2.2 each group cluster
./cluster.sh <BED_FILES_DIRECTORY> <OUTPUT_DIRECTORY>
./cluster.sh high_confidence_pas_folder high_confidence_pas_folder/pas.bed


#2.3 merge different grop
#if you only have one group pass this command
./combind.sh high_confidence_pas_folder high_confidence_pas_folder


#2.4 
python process_last.py --input_file s_site.bed --output_file final_site_for_quantification.bed

```
### 3.using salmon to quantification

```bash
#3.1 build salmon index

./get_salmon_index.sh /final_site_for_quantification.bed /RefSeq_UTR_final.bed /RefSeq_UTR_lastexon_final.bed /hg38.fa /quant_result/3UTRisoforms_sequences.fa /mnt/pengc/APA_project/alogrithm/single_cell/quant_result/3UTRisoforms_library

#3.2 quantification
./get_quant.sh 3UTRisoforms_library clean_fastq_dir quant_result_dir
#merge quant result from different sample
./merge_quant.sh sample.txt quant_result_dir/final_quant_result.txt

```

### 4.get data matrix for downstream analysis
```bash
python get_final_result.py --group_files group_A.txt group_B.txt ... --merge_file final_quant_result.txt --output_dir final_result 
```

