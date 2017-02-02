#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Filter sequences by length and mean quality value.')
parser.add_argument(
    '-i', metavar='in_format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-o', metavar='out_format', type=str, help="Output format (fastq).", default='fastq')
parser.add_argument(
    '-q', metavar='min_qual', type=float, help="Minimum mean quality value (0.0).", default=0.0)
parser.add_argument(
    '-l', metavar='min_length', type=int, help="Minimum length (0).", default=0)
parser.add_argument(
    '-u', metavar='max_length', type=int, help="Maximum length (None).", default=None)
parser.add_argument('input_fastx', nargs='?', help='Input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastx', nargs='?', help='Output file (default: stdout).',
                    type=argparse.FileType('w'), default=sys.stdout)


def record_filter(input_iter, in_format, min_qual, min_len, max_len):
    """ Filter SeqRecord objects by length and mean quality.

    :param input_iter: Iterator of SeqRecord objects.
    :param in_format: Input format.
    :param min_qual: Minimum mean quality.
    :param min_len: Minimum length.
    :param max_len: Maximum length.
    :returns: SeqRecord object.
    :rtype: generator
    """
    for record in input_iter:
        # Quality filtering:
        if in_format == 'fastq':
            mean_quality = seq_util.mean_qscore(record.letter_annotations["phred_quality"])
            if mean_quality < min_qual:
                continue
        # Length filtering:
        if len(record) < min_len:
            continue
        if max_len is not None and len(record) > max_len:
            continue
        yield record

if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=args.i)

    output_iterator = record_filter(input_iterator, args.i, args.q, args.l, args.u)

    seq_util.write_seq_records(output_iterator, args.output_fastx, format=args.o)
