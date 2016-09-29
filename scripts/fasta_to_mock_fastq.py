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

def add_qualities(input_iter, mock_qual):
    for record in input_iter:
        yield seq.mock_qualities(record, mock_qual)

if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = SeqIO.parse(sys.stdin, 'fasta')
    SeqIO.write(add_qualities(input_iterator, args.q), sys.stdout, 'fastq')
