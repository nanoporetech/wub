#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from Bio import SeqIO
import numpy as np

from wub.simulate import genome as sim_genome
from wub.util import parse as parse_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Simulate genome sequence with the specified number of chromosomes, length distribution (truncated gamma) and base composition. Output is fasta to the standard output.')
parser.add_argument(
    '-n', metavar='nr_chrom', type=int, help="Number of chromosomes (23).", default=23)
parser.add_argument('-m', metavar='mean_length', type=int,
                    help="Mean length of chromosomes (5000000).", default=5000000)
parser.add_argument(
    '-a', metavar='gamma_shape', type=float, help="Gamma shape parameter (1).", default=1.0)
parser.add_argument(
    '-l', metavar='low_trunc', type=int, help="Lower truncation point (None).", default=None)
parser.add_argument(
    '-u', metavar='high_trunc', type=int, help="Upper truncation point (None).", default=None)
parser.add_argument('-b', metavar='base_freqs', type=str,
                    help="Relative base frequencies in A,C,G,T order (1,1,1,1).", default="1,1,1,1")


if __name__ == '__main__':
    args = parser.parse_args()

    base_frequencies = np.array(parse_util.separated_list_to_floats(args.b))
    # Normalise relative base frequencies to probabilities:
    base_frequencies = base_frequencies / np.sum(base_frequencies)

    simulation_iterator = sim_genome.simulate_genome(
        args.n, args.m, args.a, args.l, args.u, base_frequencies)
    SeqIO.write(simulation_iterator, sys.stdout, 'fasta')
