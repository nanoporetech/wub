#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from Bio.Seq import Seq

from wub.simulate import seq as sim_seq
from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Add a specified number of errors to random sites for each input sequence.""")
parser.add_argument('-n', metavar='nr_errors', type=int,
                    help="Number of errors to introduce (0).", default=0)
parser.add_argument('-t', metavar='error_type', type=str,
                    help="Error type: substitution, insertion or deletion.", choices=['substitution', 'insertion', 'deletion'], default='substitution')
parser.add_argument('input_fasta', nargs='?', help='Input fasta (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fasta', nargs='?', help='Output fasta (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)


def add_fixed_errors(input_iter, nr_errors, error_type):
    """Simulate sequencing errors for each SeqRecord object in the input iterator.

    :param input_iter: Iterator of SeqRecord objects.
    :para nr_errors: Number of errors to introduce.
    :param error_type: Error type: substitution, insertion or deletion.
    :returns: Generator of SeqRecord objects.
    :rtype: generator
    """
    for record in input_iter:
        mutated_seq = sim_seq.add_errors(record.seq, nr_errors, error_type)
        record.seq = Seq(mutated_seq)
        yield record

if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(args.input_fasta, format='fasta')

    simulation_iterator = add_fixed_errors(input_iterator, args.n, args.t)

    seq_util.write_seq_records(
        simulation_iterator, args.output_fasta, format='fasta')
