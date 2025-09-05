This directory contains example data for testing purposes.  
All files are provided in simplified form and are not intended for biological interpretation.  

Files included:  

1. **sample1.bam**, **sample2.bam**, **sample3.bam**  
   - BAM files containing reads aligned to a limited set of genes:  
     *PRICKLE1, NUP107, SUPT6H, PROCA1, NUDT21, and CDC6*.  
   - These subsets are provided solely to facilitate testing and demonstration, and do not represent complete transcriptome data.  

2. **RefSeq_UTR_final.bed**  
   - BED file containing the 3′ untranslated regions (3′UTRs) annotation from RefSeq. Connected 3'UTRs are merged into one 3’UTR transcript, while unconnected ones are set as different 3’UTR transcripts. 3’UTR transcripts are named as “Gene%%Transcript number”.   

3. **RefSeq_UTR_lastexon_final.bed**  
   - BED file including the last exon in addition to the 3′UTRs (RefSeq_UTR_final.bed).  

