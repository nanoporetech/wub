#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import numpy as np
from Bio.Seq import Seq

from wub.simulate import seq as sim_seq
from wub.simulate import genome as sim_genome
from wub.util import parse as parse_util
from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Simulate sequencing errors for each input sequence.
    """)
parser.add_argument('-e', metavar='error_rate', type=float,
                    help="Total rate of substitutions insertions and deletions (0.1).", default=0.1)
parser.add_argument('-w', metavar='error_weights', type=str,
                    help="Relative frequency of substitutions,insertions,deletions (1,1,4).", default="1,1,4")
parser.add_argument('input_fasta', nargs='?', help='Input fasta (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fasta', nargs='?', help='Output fasta (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)


def simulate_errors(input_iter, error_rate, error_weights):
    """Simulate sequencing errors for each SeqRecord object in the input iterator.

    :param input_iter: Iterator of SeqRecord objects.
    :para error_rate: Total error rate of substitutions, insertions and deletions.
    :param error_weights: Relative frequency of substitutions,insertions,deletions.
    :returns: Generator of SeqRecord objects.
    :rtype: generator
    """
    for record in input_iter:
        mutated_seq = sim_seq.simulate_sequencing_errors(record.seq, error_rate, error_weights).seq
        record.seq = Seq(mutated_seq)
        yield record

if __name__ == '__main__':
    args = parser.parse_args()

    # Process error weights:
    error_weights = np.array(parse_util.separated_list_to_floats(args.w))
    # Normalise error weights to probabilities:
    error_weights = parse_util.normalise_array(error_weights)
    error_weights = dict(
        zip(['substitution', 'insertion', 'deletion'], error_weights))

    input_iterator = seq_util.read_seq_records(args.input_fasta, format='fasta')

    simulation_iterator = simulate_errors(input_iterator, args.e, error_weights)

    seq_util.write_seq_records(
        simulation_iterator, args.output_fasta, format='fasta')
