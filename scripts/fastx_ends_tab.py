#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import argparse
import sys
from os import path

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Generate a tab separated file with the first and last -n bases of the sequences.')
parser.add_argument(
    '-i', metavar='in_format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-n', metavar='nr_bases', type=int, help=".", default=100)
parser.add_argument('input_fastx', nargs='?', help='Input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_tsv', nargs='?', help='Output file (default: stdout).',
                    type=argparse.FileType('w'), default=sys.stdout)


if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=args.i)

    for rec in input_iterator:
        args.output_tsv.write("{}\t{}\t{}\n".format(rec.id, rec.seq[0:args.n], rec.seq[-args.n:]))

    args.output_tsv.flush()
    args.output_tsv.close()
