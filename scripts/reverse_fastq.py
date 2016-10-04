#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Reverse (but not complement!) sequences and qualities in fastq file.')
parser.add_argument('input_fastq', nargs='?', help='Input fastq (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastq', nargs='?', help='Output fastq (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)


def reverse_seq_records(input_iterator):
    """Reverse SeqRecord objects.

    :param input_iterator: Iterator of SeqRecord objects.
    :returns: Generator of reversed SeqRecord objects.
    :rtype: generator
    """
    for record in input_iterator:
        yield record[::-1]

if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastq, format='fastq')
    output_iterator = reverse_seq_records(input_iterator)
    seq_util.write_seq_records(
        output_iterator, args.output_fastq, format='fastq')
