#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Add sequence length to sequence record descriptions.')
parser.add_argument(
    '-i', metavar='in_format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-o', metavar='out_format', type=str, help="Output format (fastq).", default='fastq')
parser.add_argument('input_fastx', nargs='?', help='Input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastx', nargs='?', help='Output file (default: stdout).',
                    type=argparse.FileType('w'), default=sys.stdout)


def _record_annotate_length(input_iter):
    """ Add sequence length to record description.
    """
    for record in input_iter:
        record.description = record.description + " seq_length={}".format(len(record.seq))
        yield record


if __name__ == '__main__':
    args = parser.parse_args()

    if args.i == 'fasta' and args.o == 'fastq':
        sys.stderr.write(
            "Cannot produce fastq output from fasta! Use fasta_to_mock_fastq.py instead.\n")
        sys.exit(1)

    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=args.i)

    output_iterator = _record_annotate_length(input_iterator)

    seq_util.write_seq_records(output_iterator, args.output_fastx, format=args.o)
