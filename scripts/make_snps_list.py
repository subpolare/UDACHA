#!/usr/bin/env python3
import sys
import multiprocessing
from collections import defaultdict
from pyfaidx import Fasta
import argparse

def init_worker(genome_path):
    global genome_global
    genome_global = Fasta(genome_path)

def process_chromosome(args):
    chrom, records = args
    chrom_seq = genome_global[chrom][:].seq
    results = []
    for rec in records:
        _, end, variant_id, ref, alt = rec
        left_seq = chrom_seq[end - 30 : end - 1]
        right_seq = chrom_seq[end:end + 30]
        result_line = f'{variant_id}\t{left_seq}[{ref}/{alt}]{right_seq}'
        results.append(result_line)
    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome', required = True, help = 'Path to genome file')
    parser.add_argument('--threads', type = int, default = 1, help = 'Number of threads to use')
    parser.add_argument('--input', required = True, help = 'Input file with variants')
    args = parser.parse_args()
    records_by_chrom = defaultdict(list)
    with open(args.input, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 6:
                continue
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            variant_id = fields[3]
            ref = fields[4]
            alt = fields[5]
            records_by_chrom[chrom].append((start, end, variant_id, ref, alt))
    tasks = list(records_by_chrom.items())
    pool = multiprocessing.Pool(processes = args.threads, initializer = init_worker, initargs = (args.genome,))
    results = pool.map(process_chromosome, tasks)
    pool.close()
    pool.join()
    for chrom_results in results:
        for line in chrom_results:
            print(line)

if __name__ == '__main__':
    main()  
