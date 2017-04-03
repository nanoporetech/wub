#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Filter out sequences present in the first file from the second file.')
parser.add_argument(
    '-i', metavar='in_format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-o', metavar='out_format', type=str, help="Output format (fastq).", default='fastq')
parser.add_argument('input_fastx_bait', nargs='?', help='First input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('input_fastx_target', nargs='?', help='Second input file.',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastx', nargs='?', help='Output file (default: stdout).',
                    type=argparse.FileType('w'), default=sys.stdout)


def _record_filter(input_iter_bait, input_iter_target):
    """ Filter out SeqRecord objects present in the first iterator. """
    bait_ids = [read.id for read in input_iter_bait]
    for record in input_iter_target:
        if record.id not in bait_ids:
            yield record


if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator_bait = seq_util.read_seq_records(
        args.input_fastx_bait, format=args.i)

    input_iterator_target = seq_util.read_seq_records(
        args.input_fastx_target, format=args.i)

    output_iterator = _record_filter(input_iterator_bait, input_iterator_target)

    seq_util.write_seq_records(output_iterator, args.output_fastx, format=args.o)
