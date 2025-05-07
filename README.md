# APAcatcher

## Install APAcatcher
We recommend to creat a conda enviroment firstly using conda create -n apacatcher_env python=3.8 and then conda activate apacatcher_env
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
Example of original depth_file
```bash
chr1    70009   0
chr1    70010   0
chr1    70011   0
chr1    70012   0
chr1    70013   0
chr1    70014   0
...
```
First you need add gene%%transcript and strand based on RefSeq_UTR_final.bed
```bash
python add_geneinfo.py -g RefSeq_UTR_final.bed -d depth_file_dir
```
An example of the depth file after adding the information of gene%%transcript and strand
```bash
chr1    70009   0       OR4F5%%1        +
chr1    70010   0       OR4F5%%1        +
chr1    70011   0       OR4F5%%1        +
chr1    70012   0       OR4F5%%1        +
chr1    70013   0       OR4F5%%1        +
chr1    70014   0       OR4F5%%1        +
chr1    70015   0       OR4F5%%1        +
chr1    70016   0       OR4F5%%1        +
...
```
### 2.get high confidence APA sites
The options for running APAcatcher for PAS identification
```bash
--input_folder                  'Path to the input folder containing .txt files.'
--genome_file                   'Path to the genome fasta file.'
--output_folder                 'Path to the output folder where results will be saved.'
--tpm_threshold     default=1   'Threshold for the tpm.'
--length_threshold  default=100 'Threshold for the length of 3'UTR.'
--penalty           default=50  'Penalty value for change point detection.'
--min_size          default=30  'Minimum size for change point detection.'
--num_processes     default=4   'Number of parallel processes to use'
```
Example prediction from depth file
```bash
#2.1 Using PELT and DL model get PloyA sites
python main.py --input_folder depth_file_dir --genome_file hg38.fa --output_folder high_confidence_pas_folder

#2.2 Cluster the sites obtained within each group.
./cluster_bed_files.sh <BED_FILES_DIRECTORY> <OUTPUT_DIRECTORY>
./cluster_bed_files.sh high_confidence_pas_folder high_confidence_pas_folder/pas.bed


#2.3 combind the sites obtained from different groups.
#if you only have one group pass this command
./combind.sh -i high_confidence_pas_folder high_confidence_pas_folder/combind_pas_site.bed

#2.4 Remove sites located within 100 bp of annotated sites.
#if you only have one group use this bed file
python process_last.py --input_file pas_site.bed --output_file final_site_for_quantification.bed
#if you have more than one group use this bed file
python process_last.py --input_file combind_pas_site.bed --output_file final_site_for_quantification.bed

```
Example of high confidence APA sites bed file
```bash
#chr    #start    #end     #cluster size         #gene%%transcript     #strand
chr1    944201    944201    4                     NOC2L%%1              -
chr1    965721    965721    2                     KLHL17%%1             +
chr1    1014542   1014542   2                     ISG15%%1              +
chr1    1056118   1056118   2                     AGRN%%1               +
chr1    1081822   1081822   2                     C1orf159%%1           -
chr1    1216930   1216930   2                     SDF4%%1               -
...
```
### 3.using salmon to quantification

```bash
#3.1 build salmon index

./get_salmon_index.sh /final_site_for_quantification.bed /RefSeq_UTR_final.bed /RefSeq_UTR_lastexon_final.bed /hg38.fa \
/quant_result/3UTRisoforms_sequences.fa /quant_result/3UTRisoforms_library

#3.2 quantification
./get_quant.sh 3UTRisoforms_library clean_fastq_dir quant_result_dir
#merge quant result from different sample
./merge_quant.sh sample.txt quant_result_dir/final_quant_result.txt

```

### 4.get data matrix for downstream analysis
```bash
python get_final_result.py --group_files group_A.txt group_B.txt ... --merge_file final_quant_result.txt --output_dir final_result 
```
Example of 3UTR_index file
```bash
#Name                                   Length  Transcript      start           end           strand    sample1_indexUTR        sample2_indexUTR          sample3_indexUTR
AACS%%3:+::chr12:125142091-125142380    289     AACS%%3         125142091       125142380     +         0.6331098039907349      0.8141000731552296        0.6945178987151257
AAMP%%1:-::chr2:218264128-218264608     480     AAMP%%1         218264128       218264608     -         0.9140432121116838      0.9762385375962184        0.9821115355165274
...
```
