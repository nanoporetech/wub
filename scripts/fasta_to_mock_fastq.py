#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq
from Bio import SeqIO

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Convert fasta file to fastq with mock qualities.')
parser.add_argument(
    '-q', metavar='mock_quals', type=int, help="Mock quality value (40).", default=40)
parser.add_argument('input_fasta', nargs='?', help='Input fasta (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastq', nargs='?', help='Output fastq (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)


if __name__ == '__main__':
    args = parser.parse_args()

    mock_qual = args.q

    input_iterator = SeqIO.parse(args.input_fasta, 'fasta')
    output_iterator = (seq.mock_qualities(record, mock_qual) for record in input_iterator)
    SeqIO.write(output_iterator, args.output_fastq, 'fastq')
