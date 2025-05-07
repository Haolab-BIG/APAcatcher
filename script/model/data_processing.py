#Author: [ChengPeng]
import csv
import numpy as np
import pandas as pd
import ruptures as rpt
from ruptures.base import BaseCost
from collections import defaultdict
import argparse
from multiprocessing import Pool
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from ruptures.base import BaseCost
from collections import defaultdict




class PoissonLoss(BaseCost):
    model = ""

    def __init__(self, min_size=10):
        self.min_size = min_size

    def fit(self, signal):
        self.signal = np.array(signal)
        return self

    def error(self, start, end):
        segment = self.signal[start:end]
        mean = np.mean(segment)
        return -np.sum(segment * np.log(mean + 1e-10) - mean + 1e-10)


def group_and_filter_by_gene(input_file, tmp_threshold):
    gene_dict = defaultdict(lambda: {"counts": [], "positions": [], "chromosome": "", "strand": ""})
    total_sum = 0

    with open(input_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            chromosome = row[0]
            position = int(row[1])
            count_value = int(row[2])
            gene_name = row[3]
            strand = row[4]

            gene_dict[gene_name]["counts"].append(count_value)
            gene_dict[gene_name]["positions"].append(position)
            gene_dict[gene_name]["chromosome"] = chromosome
            gene_dict[gene_name]["strand"] = strand
            total_sum += count_value

    threshold = total_sum * tmp_threshold // 10**6
    filtered_gene_dict = {gene: data for gene, data in gene_dict.items()
                          if sum(data["counts"]) > threshold}
    sorted_filtered_gene_dict = dict(sorted(filtered_gene_dict.items(),
                                            key=lambda item: len(item[1]["counts"]),
                                            reverse=True))
    return sorted_filtered_gene_dict



def filter_change_points_recursively(cp, coverage, strand):
    while True:
        initial_cp_len = len(cp)
        cp, coverage = filter_change_points_around_well(cp, coverage, strand)
        if len(cp) == initial_cp_len:
            break
        elif len(cp) <=1:
            break
    return cp, coverage


def filter_change_points_around_well(cp, coverage, strand):
    M = [sum(segment) / (len(segment) + 1) for segment in coverage]

    def remove_change_point(index):
        if 0 <= index < len(cp):
            cp.pop(index)
            coverage.pop(index)

    if strand == '+':
        i = 1
        while i < len(M) - 1:
            if M[i - 1] > M[i] < M[i + 1]:
                if M[i - 1] == M[i + 1]:
                    remove_change_point(i)
                    remove_change_point(i)
                elif M[i - 1] > M[i + 1]:
                    remove_change_point(i + 1)
                else:
                    remove_change_point(i - 1)
            i += 1
    else:
        i = 1
        while i < len(M) - 1:
            if M[i - 1] > M[i] < M[i + 1]:
                if M[i - 1] == M[i + 1]:
                    remove_change_point(i)
                    remove_change_point(i - 1)
                elif M[i - 1] < M[i + 1]:
                    remove_change_point(i - 1)
                else:
                    remove_change_point(i)
            i += 1

    return cp, coverage


def filter_redundant_change_points(cp, coverage, strand):
    i = 0
    while i < len(cp) - 1:
        start_index = max(0, cp[i] - 10)  
        end_index = min(len(coverage), cp[i] + 10)  
        
        if strand == '+':
            if ((sum(coverage[start_index:cp[i]]) // (cp[i] - start_index)) > 
                (sum(coverage[cp[i]:end_index]) // (end_index - cp[i]))):
                cp.pop(i)
            else:
                i += 1
        else:
            if ((sum(coverage[start_index:cp[i]]) // (cp[i] - start_index)) < 
                (sum(coverage[cp[i]:end_index]) // (end_index - cp[i]))):
                cp.pop(i)
            else:
                i += 1
    return cp

def process_gene(gene, data, penalty, min_size, threshold=10):
    try:
        # Extract relevant data
        signal = np.array(data["counts"])
        positions = data["positions"]
        chromosome = data["chromosome"]
        strand = data["strand"]
        utr_all_length = len(signal)
        if utr_all_length < threshold:
            # Handle short signals directly
            results = []
            retained_positions = []
            if strand == '+':
                last_position = positions[-1]
                retained_positions.append([chromosome, last_position, last_position + 1, 1, gene, strand])
            elif strand == '-':
                first_position = positions[0]
                retained_positions.append([chromosome, first_position, first_position + 1, 1, gene, strand])
            return results, retained_positions

        # Run PELT algorithm
        algo = rpt.Pelt(custom_cost=PoissonLoss(min_size=min_size)).fit(signal)
        best_bkps = algo.predict(pen=penalty)

        # Filter change points
        
        best_bkps_filtered_1, _ = filter_change_points_recursively(
            best_bkps,
            [signal[i:j] for i, j in zip([0] + best_bkps[:-1], best_bkps)] + [signal[best_bkps[-1]:]],
            strand
        )
        
        #best_bkps_filtered_2 = filter_redundant_change_points(best_bkps_filtered_1, signal, strand)
        # Generate results and retained positions
        
        results = []
        retained_positions = []
        for bkp in best_bkps_filtered_1:
            if bkp < len(positions):
                position = positions[bkp - 1]
                if strand == '+':
                    results.append([chromosome, position, position + 1, 1, gene, strand])
                elif strand == '-':
                    results.append([chromosome, position, position + 1, 1, gene, strand])

        # Ensure the correct position is retained based on the strand
        if strand == '+':
            last_position = positions[-1]
            retained_positions.append([chromosome, last_position, last_position + 1, 1, gene, strand])
        elif strand == '-':
            first_position = positions[0]
            retained_positions.append([chromosome, first_position, first_position + 1, 1, gene, strand])

        return results, retained_positions

    except Exception as e:
        print(f"Error processing gene {gene}: {e}")
        return [], []
