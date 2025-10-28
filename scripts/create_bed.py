#!/usr/bin/env python3

import sys
import gzip
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description = 'Convert VCF to BED with specific columns.')
    parser.add_argument('-i', '--input', required = True, help = 'Input VCF file (can be gzipped).')
    parser.add_argument('-o', '--output', required = True, help = 'Output BED file.')
    return parser.parse_args()

def open_vcf(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

args = parse_args()
with open_vcf(args.input) as vcf, open(args.output, 'w') as bed:
    bed_fields = [
                '#chr',
                'start',
                'end',
                'id',
                'ref',
                'alt',
                'ref_count',
                'alt_count',
                'sample_id'
            ]
    bed.write('\t'.join(bed_fields) + '\n')
    sample_ids = []
    for line in vcf:
        line = line.strip()
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            headers = line.split('\t')
            sample_ids = headers[9:]
            continue
        if not sample_ids:
            print('Error: No sample columns found in VCF.', file=sys.stderr)
            sys.exit(1)
        
        fields = line.split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        var_id = fields[2]
        ref = fields[3]
        alt = fields[4]
        format_field = fields[8]
        sample_fields = fields[9:]

        start = pos - 1
        end = pos

        format_keys = format_field.split(':')
        try:
            ad_index = format_keys.index('AD')

        except ValueError:
            ad_index = None

        for sample_id, sample in zip(sample_ids, sample_fields):
            if ad_index is not None:
                sample_values = sample.split(':')
                if ad_index < len(sample_values):
                    ad = sample_values[ad_index]
                    if ad != '.' and ad != './.':
                        ad_counts = ad.split(',')
                        if len(ad_counts) >= 2:
                            ref_count = ad_counts[0]
                            alt_count = ad_counts[1]
                        elif len(ad_counts) == 1:
                            ref_count = ad_counts[0]
                            alt_count = '0'
                        else:
                            ref_count = alt_count = '0'
                    else:
                        ref_count = alt_count = '0'
                else:
                    ref_count = alt_count = '0'
            else:
                ref_count = alt_count = '0'
            
            bed_id = var_id if var_id != '.' else '.'

            bed_fields = [
                chrom,
                str(start),
                str(end),
                bed_id,
                ref,
                alt,
                ref_count,
                alt_count,
                sample_id
            ]
            
            bed.write('\t'.join(bed_fields) + '\n')
