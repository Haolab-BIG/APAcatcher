# APAcatcher
## Install APAcatcher
Running APAcatcher requires the following packages be installed
```bash
python==3.8.18
numpy==1.24.3
pandas==2.0.3
scipy==1.10.1
torch==2.2.1
ruptures==1.1.9
biopython==1.83
tqdm==4.66.5
```
## run APAcatcher in 4 steps
### 1.generate input files
```bash
python add_geneinfo.py -g RefSeq_UTR_final.bed -d depth_file_dir -p 8
```
### 2.get high confidence APA sites
The options for running APAcatcher for PAS identification
```bash
--input_folder                  'Path to the input folder containing .txt files.'
--genome_file                   'Path to the genome fasta file.'
--output_folder                 'Path to the output folder where results will be saved.'
--tpm_threshold     default=1 'Threshold for the tpm in gene data.'
--length_threshold  default=100 'Threshold for the length of gene data.'
--penalty           default=50  'Penalty value for change point detection.'
--min_size          default=30  'Minimum size for change point detection.'
--num_processes     default=4   'Number of parallel processes to use'
```
Example prediction from depth file
```bash
#2.1 Using PELT and DL model get PloyA sites
python main.py --input_folder depth_file_dir --genome_file hg38.fa --output_folder high_confidence_pas_folder --tpm_threshold 1 --length_threshold 100  --penalty 50 --min_size 30 --num_processes 8

#2.2 each group cluster
./cluster_bed_files.sh <BED_FILES_DIRECTORY> <OUTPUT_DIRECTORY>
./cluster_bed_files.sh high_confidence_pas_folder high_confidence_pas_folder/pas.bed


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

