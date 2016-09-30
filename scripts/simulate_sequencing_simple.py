#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from Bio import SeqIO
import numpy as np

from wub.simulate import seq as sim_seq
from wub.util import parse as parse_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Sample fragments from the input genome and simulate sequencing errors.')
parser.add_argument(
    '-n', metavar='nr_reads', type=int, help="Number of simulated reads (1).", default=1)
parser.add_argument('-m', metavar='mean_length', type=int,
                    help="Mean read length (5000).", default=5000)
parser.add_argument(
    '-a', metavar='gamma_shape', type=float, help="Read length distribution: gamma shape parameter (1).", default=1.0)
parser.add_argument(
    '-l', metavar='low_trunc', type=int, help="Read length distribution: lower truncation point (None).", default=None)
parser.add_argument(
    '-e', metavar='error_rate', type=float, help="Total rate of substitutions insertions and deletions (0.1).", default=0.1)
parser.add_argument('-w', metavar='error_weights', type=str,
                    help="Relative frequency of substitutions,insertions,deletions (1,1,4).", default="1,1,4")
parser.add_argument('input_fasta', nargs='?', help='Input genome in fasta format (default: stdin).', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastq', nargs='?', help='Output fastq (default: stdout)', type=argparse.FileType('w'), default=sys.stdout)


if __name__ == '__main__':
    args = parser.parse_args()

    for line in args.input_fasta:
        print line
