import sys, os, gzip, argparse, csv
from collections import defaultdict
from tqdm.auto import tqdm

BED_HEADER = [
    '#chr',
    'start',
    'end',
    'id',
    'ref',
    'alt',
    'ref_count',
    'alt_count',
    'sample_id',
]

def parse_args():
    parser = argparse.ArgumentParser(
        description='Convert clusters of VCF files (from metadata.clustered.tsv) to one BED file per indiv_id.'
    )
    parser.add_argument(
        '-m', '--metadata',
        required=True,
        help='TSV metadata file (e.g. clustering/metadata.clustered.tsv).'
    )
    parser.add_argument(
        '--work',
        required=True,
        help='Working directory containing VCFs, BEDs, tmp.'
    )
    parser.add_argument(
        '--indiv-column',
        default='indiv_id',
        help='Column name with cluster id (default: indiv_id).'
    )
    parser.add_argument(
        '--path-column',
        default='path',
        help='Column name with VCF path (default: path).'
    )
    return parser.parse_args()


def open_vcf(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    return open(file_path, 'r')


def make_vcf_path(path_from_metadata, work_dir):
    name = os.path.basename(path_from_metadata)
    if name.endswith('.without_MAF.vcf.gz'):
        name = name.replace('.without_MAF.vcf.gz', '.vcf.gz')
    elif '.without_MAF' in name:
        name = name.replace('.without_MAF', '')
    return os.path.join(work_dir, 'VCFs', name)


def extract_variants_from_vcf(vcf_path):
    records = []
    try:
        vcf = open_vcf(vcf_path)
    except OSError as e:
        print(f'Warning: cannot open VCF {vcf_path}: {e}', file=sys.stderr)
        return records

    with vcf:
        sample_ids = []
        for line in vcf:
            line = line.strip()
            if not line:
                continue
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                headers = line.split('\t')
                sample_ids = headers[9:]
                continue
            if not sample_ids:
                print(f'Error: No sample columns found in VCF {vcf_path}.', file=sys.stderr)
                return records

            fields = line.split('\t')
            if len(fields) < 8:
                continue

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
                if bed_id == '.':
                    continue

                bed_fields = [
                    chrom,
                    str(start),
                    str(end),
                    bed_id,
                    ref,
                    alt,
                    ref_count,
                    alt_count,
                    sample_id,
                ]
                records.append(bed_fields)
    return records


def load_clusters(metadata_path, indiv_col, path_col, work_dir):
    clusters = defaultdict(list)

    with open(metadata_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if indiv_col not in reader.fieldnames or path_col not in reader.fieldnames:
            print(
                'Error: metadata file must contain columns '
                f'"{indiv_col}" and "{path_col}". '
                f'Found columns: {reader.fieldnames}',
                file=sys.stderr
            )
            sys.exit(1)

        for row in reader:
            indiv_id = row[indiv_col]
            path_from_meta = row[path_col]

            if not indiv_id or not path_from_meta or path_from_meta == 'nan':
                continue

            vcf_path = make_vcf_path(path_from_meta, work_dir)

            if not os.path.exists(vcf_path):
                print(f'Warning: VCF file does not exist, skipping: {vcf_path}', file=sys.stderr)
                continue

            clusters[indiv_id].append(vcf_path)

    return clusters


def chrom_sort_key(chrom):
    c = chrom
    if c.lower().startswith('chr'):
        c = c[3:]
    try:
        n = int(c)
        if 1 <= n <= 22:
            return (0, n)
        else:
            return (1, n)
    except ValueError:
        return (2, c)


def main():
    args = parse_args()

    work_dir = os.path.abspath(args.work)
    outdir = os.path.join(work_dir, 'BEDs')
    os.makedirs(outdir, exist_ok=True)

    clusters = load_clusters(
        metadata_path=args.metadata,
        indiv_col=args.indiv_column,
        path_col=args.path_column,
        work_dir=work_dir,
    )

    if not clusters:
        print('No clusters found (check metadata / paths).', file=sys.stderr)
        sys.exit(1)

    for indiv_id, vcf_paths in tqdm(clusters.items(), desc='Clusters'):
        if not vcf_paths:
            continue

        out_bed_path = os.path.join(outdir, f'{indiv_id}.bed')

        all_records = []
        for vcf_path in vcf_paths:
            all_records.extend(extract_variants_from_vcf(vcf_path))

        all_records.sort(key=lambda row: (*chrom_sort_key(row[0]), int(row[1])))

        with open(out_bed_path, 'w') as bed:
            bed.write('\t'.join(BED_HEADER) + '\n')
            for rec in all_records:
                bed.write('\t'.join(rec) + '\n')


if __name__ == '__main__':
    main()
