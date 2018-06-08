#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Generate a tab separated file with the sequence lengths in the input file.')
parser.add_argument(
    '-i', metavar='in_format', type=str, help="Input format (fasta).", default='fasta')
parser.add_argument('input_fastx', nargs='?', help='Input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_tsv', nargs='?', help='Output file (default: stdout).',
                    type=argparse.FileType('w'), default=sys.stdout)

if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=args.i)

    args.output_tsv.write("Reference\tLength\n")
    for rec in input_iterator:
        args.output_tsv.write("{}\t{}\n".format(rec.id, len(rec.seq)))

    args.output_tsv.flush()
    args.output_tsv.close()
