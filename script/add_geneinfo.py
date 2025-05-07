#Author: [ChengPeng]
import os
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import logging
import pandas as pd
import tempfile
import shutil

# Configure logging with timestamp and severity level
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_gene_info(gene_file):
    """
    Load gene information from BED-formatted file and create a position-to-gene mapping dictionary.
    """
    # Read BED file columns: chrom, start, end, gene_name, strand
    gene_df = pd.read_csv(gene_file, sep='\t', header=None, usecols=[0, 1, 2, 3, 5])
    gene_df.columns = ['chrom', 'start', 'end', 'gene_name', 'strand']

    # Create dictionary mapping (chromosome, position) to (gene_name, strand)
    gene_info = {}
    for _, row in gene_df.iterrows():
        for pos in range(row['start'], row['end'] + 1):
            gene_info[(row['chrom'], pos)] = (row['gene_name'], row['strand'])

    return gene_info


def safe_write_csv(df, original_file):
    """
    Safely write DataFrame to file using temporary file strategy to prevent data corruption.
    """
    try:
        # Get system default temporary directory
        temp_dir = tempfile.gettempdir()

        # Create temporary file in system temporary directory
        with tempfile.NamedTemporaryFile('w', dir=temp_dir, delete=False, suffix='.tmp') as tmpfile:
            tmp_file_name = tmpfile.name
            df.to_csv(tmp_file_name, sep='\t', header=False, index=False)

        # Use secure file copy operation
        shutil.copy(tmp_file_name, original_file)

        # Clean up temporary file
        os.remove(tmp_file_name)

        # Remove original file after successful copy
        os.remove(original_file)
    except Exception as e:
        logging.error(f"Failed to replace file: {original_file}, Error: {e}")
        # Attempt to clean up temporary file if operation failed
        if 'tmp_file_name' in locals() and os.path.exists(tmp_file_name):
            os.remove(tmp_file_name)
        raise


def add_gene_info(data_file, gene_info):
    """
    Annotate data file with gene information based on genomic coordinates.
    """
    try:
        # Read data file with three required columns: chrom, pos, count
        data_df = pd.read_csv(data_file, sep='\t', header=None, names=['chrom', 'pos', 'count'])
    except Exception as e:
        logging.error(f"Failed to read file: {data_file}, Error: {e}")
        return f"Error reading {data_file}: {e}"

    # Create list of (chrom, pos) tuples for lookup
    keys = list(zip(data_df['chrom'], data_df['pos']))
    
    # Get gene information for each position
    gene_names, strands = zip(*[gene_info.get(key, ("N/A", "N/A")) for key in keys])

    # Add gene annotation columns
    data_df['gene_name'] = gene_names
    data_df['strand'] = strands

    try:
        # Write updated data back to original file
        safe_write_csv(data_df, data_file)
    except Exception as e:
        return f"Error writing {data_file}: {e}"

    return data_file


def process_file(data_file, gene_info):
    """
    Process single file with gene annotation workflow.
    """
    try:
        return add_gene_info(data_file, gene_info)
    except Exception as e:
        return f"Error processing {data_file}: {e}"


def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(description='Annotate genomic data files with gene information')
    parser.add_argument('-g', '--gene_file', required=True, 
                        help='Path to gene annotation file in BED format')
    parser.add_argument('-d', '--data_dir', required=True,
                        help='Directory containing coverage data files to process')
    parser.add_argument('-p', '--processes', type=int, default=4,
                        help='Number of parallel threads to use (default: 4)')

    args = parser.parse_args()

    # Identify target data files with specific naming pattern
    data_files = [os.path.join(args.data_dir, filename) 
                 for filename in os.listdir(args.data_dir) 
                 if filename.endswith('.txt') and '_read_coverage' in filename]
    
    logging.info("Starting gene information dictionary construction.")

    # Build gene position lookup table
    gene_info = load_gene_info(args.gene_file)
    logging.info("Gene information dictionary construction completed.")

    # Process files in parallel with progress tracking
    with ThreadPoolExecutor(max_workers=args.processes) as executor:
        futures = {executor.submit(process_file, f, gene_info): f for f in data_files}

        # Monitor progress using progress bar
        for future in tqdm(as_completed(futures), total=len(futures)):
            data_file = futures[future]
            try:
                result = future.result()
                print(f"Processed {result}")
            except Exception as e:
                print(f"Error processing {data_file}: {e}")
    
    logging.info("All data files have been successfully processed.")


if __name__ == "__main__":
    main()
