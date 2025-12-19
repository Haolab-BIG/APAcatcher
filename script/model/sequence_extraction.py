# Author: [ChengPeng]
import pandas as pd
from Bio import SeqIO


def get_sequence(chromosome, start, strand, genome, flanking=100):
    seq = genome[chromosome].seq
    start_pos = max(0, start - flanking)  # Ensure not to go out of bounds
    end_pos = start + flanking  # Extract sequence from start position

    extracted_seq = seq[start_pos:end_pos + 1]
    if strand == '-':
        return str(extracted_seq.reverse_complement())
    else:
        return str(extracted_seq)




def get_sequence_A(chromosome, start, strand, genome):

    if strand == '+':
        adj_start = start - 201
        adj_end   = start + 100
    else:  # strand == '-'
        adj_start = start - 100
        adj_end   = start + 201

    # Boundary protection
    seq_record = genome[chromosome].seq
    adj_start = max(0, adj_start)
    adj_end = min(len(seq_record), adj_end)

    # Extract DNA sequence
    extracted_seq = seq_record[adj_start:adj_end]

    # Reverse complement for negative strand
    if strand == '-':
        return str(extracted_seq.reverse_complement())
    else:
        return str(extracted_seq)


def extract_sequences_from_results(results, genome_file, species, flanking=100):
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    sequences = []
    if species == "Arabidopsis":
        for result in results:
            chromosome, start, end, _, gene, strand = result
            extracted_seq = get_sequence_A(chromosome, start, strand, genome)
            sequences.append([chromosome, start, end, gene, strand, extracted_seq])
    else:
        for result in results:
            chromosome, start, end, _, gene, strand = result
            extracted_seq = get_sequence(chromosome, start, strand, genome, flanking)
            sequences.append([chromosome, start, end, gene, strand, extracted_seq])
    return sequences
    
