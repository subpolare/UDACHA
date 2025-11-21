import argparse, csv, sys
from typing import Dict, List, Tuple


def read_bad_file(path: str) -> Tuple[Dict[str, List[Tuple[int, int, str, str, str]]], bool]:
    bad_intervals: Dict[str, List[Tuple[int, int, str, str, str]]] = {}
    has_bad_data = False

    try:
        with open(path, 'r') as bad_file:
            reader = csv.reader(bad_file, delimiter = '\t')
            header = next(reader, None)

            for row in reader:
                if not row or row[0].startswith('#'):
                    continue

                chrom = row[0]
                try:
                    interval_start = int(row[1])
                    interval_end = int(row[2])
                except ValueError:
                    raise Exception(
                        f'ValueError: Could not convert BAD interval coordinates to integer in row: {row}'
                    )

                bad_value = row[3]
                snp_count = row[4]
                sum_cover = row[6]

                if chrom not in bad_intervals:
                    bad_intervals[chrom] = []
                bad_intervals[chrom].append(
                    (interval_start, interval_end, bad_value, snp_count, sum_cover)
                )
                has_bad_data = True
    except FileNotFoundError:
        raise Exception(f'FileError: BAD file {path} not found.')

    for chrom in bad_intervals:
        bad_intervals[chrom].sort(key = lambda x: x[0])

    return bad_intervals, has_bad_data


def read_bed_file(path: str) -> Tuple[List[str], List[List[str]]]:
    bed_rows: List[List[str]] = []

    try:
        with open(path, 'r') as bed_file:
            reader = csv.reader(bed_file, delimiter = '\t')
            header = next(reader)

            for row in reader:
                if not row or row[0].startswith('#'):
                    continue
                bed_rows.append(row)
    except FileNotFoundError:
        raise Exception(f'FileError: BED file {path} not found.')

    return header, bed_rows


def merge_bed_and_bad(
    bed_header: List[str],
    bed_rows: List[List[str]],
    bad_intervals: Dict[str, List[Tuple[int, int, str, str, str]]],
    has_bad_data: bool,
) -> List[List[str]]:
    output_rows: List[List[str]] = []

    for row in bed_rows:
        if len(row) < 9:
            raise Exception(f'FormatError: BED row has fewer than 9 columns: {row}')

        if row[3] == '.':
            continue

        chrom = row[0]
        try:
            snp_start = int(row[1])
        except ValueError:
            raise Exception(
                f'ValueError: Could not convert BED start coordinate to integer in row: {row}'
            )

        bad_val = '1'
        snp_ct = '1'

        if has_bad_data:
            intervals = bad_intervals.get(chrom)
            if intervals:
                matching_interval = None
                for interval_start, interval_end, interval_bad, interval_snp_ct, sum_cov in intervals:
                    if interval_start <= snp_start < interval_end:
                        matching_interval = (interval_bad, interval_snp_ct)
                        break

                if matching_interval is not None:
                    bad_val, snp_ct = matching_interval

        try:
            ref_count = int(row[6])
            alt_count = int(row[7])
        except ValueError:
            raise Exception(
                f'ValueError: Could not convert ref_count or alt_count to integer in row: {row}'
            )

        total_cover = ref_count + alt_count

        new_row: List[str] = []
        new_row.extend(row[0:8])
        new_row.extend([bad_val, snp_ct, str(total_cover)])
        new_row.append(row[8])

        output_rows.append(new_row)

    return output_rows


def main():
    parser = argparse.ArgumentParser(
        description = 'Merge BED and BAD files into one BED file with additional columns.'
    )
    parser.add_argument('--bed', required = True, help = 'Input BED file')
    parser.add_argument('--bad', required = True, help = 'Input BAD file with BAD calculations')
    parser.add_argument(
        '-o', '--output', required = True, help = 'Output file to write the merged table'
    )
    args = parser.parse_args()

    bad_intervals, has_bad_data = read_bad_file(args.bad)
    bed_header, bed_rows = read_bed_file(args.bed)
    output_rows = merge_bed_and_bad(bed_header, bed_rows, bad_intervals, has_bad_data)

    try:
        with open(args.output, 'w', newline = '') as out_file:
            writer = csv.writer(out_file, delimiter = '\t')
            writer.writerow(
                [
                    '#chr',
                    'start',
                    'end',
                    'id',
                    'ref',
                    'alt',
                    'ref_count',
                    'alt_count',
                    'bad',
                    'SNP_per_segment',
                    'total_cover',
                    'sample_id',
                ]
            )
            writer.writerows(output_rows)
    except Exception as e:
        raise Exception(f'WriteError: Could not write to output file {args.output}. {str(e)}')

    total_bed = len(bed_rows)
    total_out = len(output_rows)
    if total_out == 0:
        print(
            f'[WARN] Output file {args.output} has only header: no SNPs passed filters or merging.',
            file = sys.stderr,
        )
    else:
        print(
            f'[INFO] {args.bed}: {total_bed} input SNPs, {total_out} written to {args.output} '
            f'(BAD {"present" if has_bad_data else "absent, treated as BAD=1"}).',
            file = sys.stderr,
        )


if __name__ == '__main__':
    main()

