#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Filter sequence files by read name.')
parser.add_argument(
    '-i', metavar='in_format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-o', metavar='out_format', type=str, help="Output format (fastq).", default='fastq')
parser.add_argument(
    '-n', metavar='read_names', type=str, help="Comma separated list of read names to select.", default="")
parser.add_argument('input_fastx', nargs='?', help='Input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastx', nargs='?', help='Output file (default: stdout).',
                    type=argparse.FileType('w'), default=sys.stdout)


def record_filter(input_iter, in_format, read_names):
    """ Filter SeqRecord objects by length and mean quality.

    :param input_iter: Iterator of SeqRecord objects.
    :param in_format: Input format.
    :param to_alphabet: Convert to this alphabet.
    :returns: SeqRecord object.
    :rtype: generator
    """
    for record in input_iter:
        if record.id in read_names:
            yield record


if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=args.i)

    names = args.n.split(',')

    output_iterator = record_filter(input_iterator, args.i, names)

    seq_util.write_seq_records(output_iterator, args.output_fastx, format=args.o)
