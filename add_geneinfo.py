import os
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import logging
import pandas as pd
import tempfile
import shutil
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')   

def load_gene_info(gene_file):
    gene_df = pd.read_csv(gene_file, sep='\t', header=None, usecols=[0, 1, 2, 3, 5])
    gene_df.columns = ['chrom', 'start', 'end', 'gene_name', 'strand']

    gene_info = {}
    for _, row in gene_df.iterrows():
        for pos in range(row['start'], row['end'] + 1):
            gene_info[(row['chrom'], pos)] = (row['gene_name'], row['strand'])

    return gene_info
"""
def safe_write_csv(df, original_file):
    dir_name = os.path.dirname(original_file)
    base_name = os.path.basename(original_file)
    
    # 创建临时文件
    with tempfile.NamedTemporaryFile('w', dir=dir_name, delete=False, suffix='.tmp') as tmpfile:
        tmp_file_name = tmpfile.name
        df.to_csv(tmp_file_name, sep='\t', header=False, index=False)
    
    try:
        # 原子性地替换原始文件
        os.replace(tmp_file_name, original_file)
    except Exception as e:
        logging.error(f"替换文件失败: {original_file}，错误: {e}")
        # 如果替换失败，删除临时文件
        os.remove(tmp_file_name)
        raise
"""

def safe_write_csv(df, original_file):
    try:
        dir_name = os.path.dirname(original_file)
        # 使用系统推荐的临时目录
        temp_dir = tempfile.gettempdir()

        with tempfile.NamedTemporaryFile('w', dir=temp_dir, delete=False, suffix='.tmp') as tmpfile:
            tmp_file_name = tmpfile.name
            df.to_csv(tmp_file_name, sep='\t', header=False, index=False)

        # 使用 shutil.copy() 代替 shutil.move()
        shutil.copy(tmp_file_name, original_file)

        # 删除临时文件
        os.remove(tmp_file_name)

        # 如果需要，删除原始文件
        os.remove(original_file)
    except Exception as e:
        logging.error(f"替换文件失败: {original_file}，错误: {e}")
        # 如果替换失败，尝试删除临时文件
        if 'tmp_file_name' in locals() and os.path.exists(tmp_file_name):
            os.remove(tmp_file_name)
        raise

def add_gene_info(data_file, gene_info):
    try:
        data_df = pd.read_csv(data_file, sep='\t', header=None, names=['chrom', 'pos', 'count'])
    except Exception as e:
        logging.error(f"读取文件失败: {data_file}，错误: {e}")
        return f"Error reading {data_file}: {e}"

    keys = list(zip(data_df['chrom'], data_df['pos']))
    gene_names, strands = zip(*[gene_info.get(key, ("N/A", "N/A")) for key in keys])

    data_df['gene_name'] = gene_names
    data_df['strand'] = strands

    try:
        safe_write_csv(data_df, data_file)
    except Exception as e:
        return f"Error writing {data_file}: {e}"

    return data_file

def process_file(data_file, gene_info):
    try:
        return add_gene_info(data_file, gene_info)
    except Exception as e:
        return f"Error processing {data_file}: {e}"

def main():
    parser = argparse.ArgumentParser(description='Add gene information to data based on chromosome position.')
    parser.add_argument('-g', '--gene_file', required=True, help='Path to the gene information file (BED format)')
    parser.add_argument('-d', '--data_dir', required=True,
                        help='Path to the directory containing data files to be processed')
    parser.add_argument('-p', '--processes', type=int, default=4, help='Number of threads to use (default: 4)')

    args = parser.parse_args()
   # data_files = [os.path.join(args.data_dir, filename) for filename in os.listdir(args.data_dir) if filename.endswith('10M_read_coverage.txt')]
   
    data_files = [os.path.join(args.data_dir, filename) for filename in os.listdir(args.data_dir) if filename.endswith('.txt') and '_read_coverage' in filename]
    logging.info("Starting to build the gene information dictionary.")

    gene_info = load_gene_info(args.gene_file)
    logging.info("Finished building the gene information dictionary.")

    with ThreadPoolExecutor(max_workers=args.processes) as executor:
        futures = {executor.submit(process_file, f, gene_info): f for f in data_files}

        for future in tqdm(as_completed(futures), total=len(futures)):
            data_file = futures[future]
            try:
                result = future.result()
                print(f"Processed {result}")
            except Exception as e:
                print(f"Error processing {data_file}: {e}")
    logging.info("All files have been processed.")

if __name__ == "__main__":
    main()

