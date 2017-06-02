#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Convert between DNA and RNA alphabets.')
parser.add_argument(
    '-i', metavar='in_format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-o', metavar='out_format', type=str, help="Output format (fastq).", default='fastq')
parser.add_argument(
    '-D', action='store_true', help="RNA->DNA alphabet conversion.", default=False)
parser.add_argument(
    '-R', action='store_true', help="DNA->RNA alphabet conversion.", default=False)
parser.add_argument('input_fastx', nargs='?', help='Input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastx', nargs='?', help='Output file (default: stdout).',
                    type=argparse.FileType('w'), default=sys.stdout)


def record_filter(input_iter, in_format, to_alphabet):
    """ Filter SeqRecord objects by length and mean quality.

    :param input_iter: Iterator of SeqRecord objects.
    :param in_format: Input format.
    :param to_alphabet: Convert to this alphabet.
    :returns: SeqRecord object.
    :rtype: generator
    """
    for record in input_iter:
        if to_alphabet == 'DNA':
            yield seq_util.rna_record_to_dna(record)
        elif to_alphabet == 'RNA':
            yield seq_util.dna_record_to_rna(record)
        else:
            raise Exception('Invalid alphabet type')


if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=args.i)

    to_alphabet = None
    if args.D and args.R:
        sys.stderr.write("-D and -R are mutually exclusive!\n")
        sys.exit(1)
    elif not args.D and not args.R:
        sys.stderr.write("Either -D or -R must be specified!\n")
        sys.exit(1)
    elif args.D:
        to_alphabet = 'DNA'
    elif args.R:
        to_alphabet = 'RNA'

    output_iterator = record_filter(input_iterator, args.i, to_alphabet)

    seq_util.write_seq_records(output_iterator, args.output_fastx, format=args.o)
