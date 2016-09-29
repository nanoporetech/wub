#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq
from Bio import SeqIO

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Convert fasta file to fastq with mock qualities. Reads from standard input and writes to standard output.')
parser.add_argument(
    '-q', metavar='mock_quals', type=int, help="Mock quality value (40).", default=40)

if __name__ == '__main__':
    args = parser.parse_args()

    mock_qual = args.q

    input_iterator = SeqIO.parse(sys.stdin, 'fasta')
    output_iterator = (seq.mock_qualities(record, mock_qual) for record in input_iterator)
    SeqIO.write(output_iterator, sys.stdout, 'fastq')
