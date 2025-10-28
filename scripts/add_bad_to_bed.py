#!/usr/bin/env python3
import argparse
import csv
import sys

def main():
    parser = argparse.ArgumentParser(
        description = 'Merge BED and BAD files into one BED file with additional columns.'
    )
    parser.add_argument('--bed', required = True, help = 'Input BED file')
    parser.add_argument('--bad', required = True, help = 'Input BAD file with BAD calculations')
    parser.add_argument('-o', '--output', required = True, help = 'Output file to write the merged table')
    args = parser.parse_args()

    bad_intervals = {}
    try:
        with open(args.bad, 'r') as bad_file:
            reader = csv.reader(bad_file, delimiter = '\t')
            header = next(reader)
            for row in reader:
                if not row or row[0].startswith('#'):
                    continue
                chrom = row[0]
                try:
                    interval_start = int(row[1])
                    interval_end = int(row[2])
                except ValueError:
                    raise Exception(f'ValueError: Could not convert BAD interval coordinates to integer in row: {row}')
                bad_value = row[3]
                snp_count = row[4]
                sum_cover = row[6]
                if chrom not in bad_intervals:
                    bad_intervals[chrom] = []
                bad_intervals[chrom].append((interval_start, interval_end, bad_value, snp_count, sum_cover))
    except FileNotFoundError:
        raise Exception(f'FileError: BAD file {args.bad} not found.')

    for chrom in bad_intervals:
        bad_intervals[chrom].sort(key = lambda x : x[0])

    bed_rows = []
    try:
        with open(args.bed, 'r') as bed_file:
            reader = csv.reader(bed_file, delimiter = '\t')
            header = next(reader)
            for row in reader:
                if not row or row[0].startswith('#'):
                    continue
                bed_rows.append(row)
    except FileNotFoundError:
        raise Exception(f'FileError: BED file {args.bed} not found.')

    output_rows = []
    for row in bed_rows:
        if row[3] == '.':
            continue
        chrom = row[0]
        try:
            snp_start = int(row[1])
        except ValueError:
            raise Exception(f'ValueError: Could not convert BED start coordinate to integer in row: {row}')
        matching_interval = None
        if chrom in bad_intervals:
            for interval in bad_intervals[chrom]:
                interval_start, interval_end, bad_val, snp_ct, sum_cov = interval
                if interval_start <= snp_start < interval_end:
                    matching_interval = interval
                    break
        if matching_interval is None:
            continue
        bad_val = matching_interval[2]
        snp_ct = matching_interval[3]
        try:
            total_cover = int(row[6]) + int(row[7])
        except ValueError:
            raise Exception(f'ValueError: Could not convert ref_count or alt_count to integer in row: {row}')
        new_row = []
        new_row.extend(row[0:8])
        new_row.extend([bad_val, snp_ct, total_cover])
        new_row.append(row[8])
        output_rows.append(new_row)

    try:
        with open(args.output, 'w', newline = '') as out_file:
            writer = csv.writer(out_file, delimiter = '\t')
            writer.writerow([
                '#chr', 'start', 'end', 'id', 'ref', 'alt',
                'ref_count', 'alt_count', 'bad', 'SNP_per_segment', 'total_cover', 'sample_id'
            ])
            writer.writerows(output_rows)
    except Exception as e:
        raise Exception(f'WriteError: Could not write to output file {args.output}. {str(e)}')

if __name__ == '__main__':
    main()
