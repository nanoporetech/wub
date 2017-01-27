#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import numpy as np

from wub.simulate import genome as sim_genome
from wub.util import parse as parse_util
from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Simulate genome sequence with the specified number of chromosomes,
    length distribution (truncated gamma) and base composition.""")
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
                    help="Relative base frequencies in A,C,G,T order (1,1,1,1) or \"random\".", default="1,1,1,1")
parser.add_argument('-z', metavar='random_seed', type=int,
                    help="Random seed (None).", default=None)
parser.add_argument('output_fasta', nargs='?', help='Output fasta (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)


if __name__ == '__main__':
    args = parser.parse_args()

    # Set random seed:
    if args.z is not None:
        np.random.seed(args.z)

    if args.b == "random":
        base_frequencies = np.random.uniform(size=4)
        base_frequencies = base_frequencies / np.sum(base_frequencies)
    else:
        base_frequencies = np.array(parse_util.separated_list_to_floats(args.b))
    # Normalise relative base frequencies to probabilities:
    base_frequencies = parse_util.normalise_array(base_frequencies)

    simulation_iterator = sim_genome.simulate_genome(
        args.n, args.m, args.a, args.l, args.u, base_frequencies)
    seq_util.write_seq_records(simulation_iterator, args.output_fasta, format='fasta')
