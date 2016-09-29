#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Simulate genome sequence with the specified number of chromosomes, length distribution (truncated gamma) and base composition. Output is fasta to the standard output.')
parser.add_argument('-n', metavar='nr_chrom', type=int, help="Number of chromosomes (23).", default=23)
parser.add_argument('-m', metavar='mean_length', type=int, help="Mean length of chromosomes (50000000).", default=50000000)
parser.add_argument('-a', metavar='gamma_shape', type=float, help="Gamma shape parameter (1).", default=1.0)
parser.add_argument('-l', metavar='low_trunc', type=int, help="Lower truncation point (None).", default=None)
parser.add_argument('-u', metavar='high_trunc', type=int, help="Upper truncation point (None).", default=None)
parser.add_argument('-b', metavar='base_freqs', type=str, help="Relative base frequencies in A,C,G,T order (1,1,1,1).", default="1,1,1,1")


if __name__ == '__main__':
    args = parser.parse_args()
