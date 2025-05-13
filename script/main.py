"""
Author: [ChengPeng]
"""
import argparse
from multiprocessing import Pool, TimeoutError
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import logging
from model.model import PAS_CNN
from model.data_processing import group_and_filter_by_gene, process_gene
from model.sequence_extraction import extract_sequences_from_results
from model.model_inference import filter_sequences_with_model
import torch
import pandas as pd
# Set up logging
import os
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def combine_and_write_final_output(filtered_sequences, retained_positions, output_file):
    logging.info("Combining and writing final output.")
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        if hasattr(filtered_sequences, 'to_csv'):
            filtered_sequences.to_csv(file, sep='\t', header=False, index=False)
        else:
            logging.warning("filtered_sequences is not a DataFrame。")
        for pos in retained_positions:
            writer.writerow(pos)
    logging.info("Final output written to file.")


def main(input_file, genome_file, output_file, tpm_threshold, length_threshold, penalty, min_size, num_processes,
         flanking_bp):
    logging.info("Starting the main process.")

    # Step 1: Obtain results and retained_positions
    logging.info("Grouping and filtering genes.")
    filtered_gene_groups = group_and_filter_by_gene(input_file, tpm_threshold)
    num_genes = len(filtered_gene_groups)
    logging.info(f"Processing {num_genes} genes with {num_processes} processes.")
    all_genes = set(filtered_gene_groups.keys())
    results = []
    retained_positions = []
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = {executor.submit(process_gene, gene, data, penalty, min_size, length_threshold): gene for
                   gene, data in
                   filtered_gene_groups.items()}
        for future in tqdm(as_completed(futures), total=len(filtered_gene_groups), desc="Processing Genes"):
            try:
                gene = futures[future]
                res, ret = future.result()
                all_genes.discard(gene)
                results.extend(res)
                retained_positions.extend(ret)
                if len(all_genes) <= 5:
                    logging.info(f"Remaining genes: {all_genes}")
            except TimeoutError:
                logging.warning("A process timed out and was skipped.")

    # Step 2: Extract sequences
    logging.info("Extracting sequences.")
    sequences = extract_sequences_from_results(results, genome_file, flanking=flanking_bp)

    # Step 3: Filter sequences with model
    logging.info("Filtering sequences with the model.")
    model = PAS_CNN()
    model.load_state_dict(torch.load('model.pth', map_location=torch.device('cpu')))
    filtered_sequences = filter_sequences_with_model(sequences, model, max_len=201)
    # Step 4: Combine and write final output
    output_file_path = os.path.join(output_file, f"{os.path.basename(input_file).replace('.txt', '_output.bed')}")
    combine_and_write_final_output(filtered_sequences, retained_positions, output_file_path)
    logging.info("Main process completed.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process all .txt files in a directory and filter APA sites.')

    # Input and output file paths
    parser.add_argument('--input_folder', required=True, help='Path to the input folder containing .txt files.')
    parser.add_argument('--genome_file', required=True, help='Path to the genome fasta file.')
    parser.add_argument('--output_folder', required=True, help='Path to the output folder where results will be saved.')

    # Filtering and processing parameters
    parser.add_argument('--tpm_threshold', type=int, default=1, help='Threshold for the tpm.')
    parser.add_argument('--length_threshold', type=int, default=100, required=True, help='Threshold for the length of 3'UTR.')
    parser.add_argument('--penalty', type=float, default=50, help='Penalty value for change point detection.')
    parser.add_argument('--min_size', type=int, default=30, help='Minimum size for change point detection.')
    parser.add_argument('--num_processes', type=int, default=4, help='Number of parallel processes to use.')
    parser.add_argument('--flanking_bp', type=int, default=100, help='Number of flanking base pairs for sequence extraction.')

    args = parser.parse_args()

    os.makedirs(args.output_folder, exist_ok=True)

    for file_name in os.listdir(args.input_folder):
        if file_name.endswith(".txt"):
            input_file_path = os.path.join(args.input_folder, file_name)
            main(input_file_path, args.genome_file, args.output_folder, args.tpm_threshold, args.length_threshold,
                 args.penalty, args.min_size, args.num_processes, args.flanking_bp)

