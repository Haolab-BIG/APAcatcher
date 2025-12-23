# Author: [ChengPeng]
import logging
import pandas as pd
from Bio import SeqIO
from typing import Dict, List, Tuple

def get_genomic_sequence(
    chromosome: str, 
    pos: int, 
    strand: str, 
    genome_dict: Dict, 
    upstream: int = 100, 
    downstream: int = 100
) -> str:
    
    if chromosome not in genome_dict:
        logging.warning(f"Chromosome {chromosome} not found in genome.")
        return ""

    seq_record = genome_dict[chromosome].seq
    seq_len = len(seq_record)

    if strand == '+':
        start_pos = pos - upstream - 1
        end_pos = pos + downstream
    else:
        start_pos = pos - downstream - 1
        end_pos = pos + upstream

    start_final = max(0, start_final := start_pos)
    end_final = min(seq_len, end_pos)

    extracted_seq = seq_record[start_final:end_final]

    if strand == '-':
        return str(extracted_seq.reverse_complement())
    return str(extracted_seq)


def extract_sequences_from_results(
    results: List[Tuple], 
    genome_file: str, 
    flanking: int = 100, 
    species: str = "Default"
) -> List[List]:
    logging.info(f"Loading genome index from {genome_file}...")
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    
    config = {
        "Arabidopsis": {"up": 201, "down": 100},
        "Default": {"up": flanking, "down": flanking}
    }
    
    spec_cfg = config.get(species, config["Default"])
    upstream_len = spec_cfg["up"]
    downstream_len = spec_cfg["down"]

    sequences = []
    logging.info(f"Extracting sequences (Upstream:{upstream_len}, Downstream:{downstream_len})...")
    for result in results:
        chrom, start, end, _, gene, strand = result
        
        extracted_seq = get_genomic_sequence(
            chrom, start, strand, genome, 
            upstream=upstream_len, 
            downstream=downstream_len
        )
        
        sequences.append([chrom, start, end, gene, strand, extracted_seq])
        
    return sequences

