#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import numpy as np

from wub.simulate import seq as sim_seq
from wub.util import parse as parse_util
from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Simulate sequences of fixed length and specified base composition.""")
parser.add_argument(
    '-n', metavar='nr_seq', type=int, help="Number of sequences (1).", default=1)
parser.add_argument('-m', metavar='length', type=int,
                    help="Length of simulated sequences (3000).", default=3000)
parser.add_argument('-b', metavar='base_freqs', type=str,
                    help="Relative base frequencies in A,C,G,T order (1,1,1,1).", default="1,1,1,1")
parser.add_argument('-z', metavar='random_seed', type=int,
                    help="Random seed (None).", default=None)
parser.add_argument('output_fasta', nargs='?', help='Output fasta (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)


if __name__ == '__main__':
    args = parser.parse_args()

    # Set random seed:
    if args.z is not None:
        np.random.seed(args.z)

    base_frequencies = np.array(parse_util.separated_list_to_floats(args.b))
    # Normalise relative base frequencies to probabilities:
    base_frequencies = parse_util.normalise_array(base_frequencies)

    simulation_iterator = (seq_util.new_dna_record(sim_seq.simulate_sequence(args.m, base_frequencies), "seq_{}".format(i)) for i in range(args.n))

    seq_util.write_seq_records(simulation_iterator, args.output_fasta, format='fasta')
